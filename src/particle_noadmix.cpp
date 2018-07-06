
#include "particle_noadmix.h"
#include "misc.h"
#include "probability.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// default constructor for particle class under no-admixture model
particle_noadmix::particle_noadmix() {}

//------------------------------------------------
// constructor for particle class under no-admixture model
particle_noadmix::particle_noadmix(vector<vector<int>> &data, vector<int> &Jl_, vector<int> &ploidy_, Rcpp::List &args_model, double beta_raised_) {
  
  // extract data and parameters
  data_ptr = &data;
  Jl = Jl_;
  ploidy = ploidy_;
  n = ploidy.size();
  L = Jl.size();
  
  beta_raised = beta_raised_;
  K = rcpp_to_int(args_model["K"]);
  lambda = rcpp_to_double(args_model["lambda"]);
  
  // initialise group allocation
  group = vector<int>(n);
  
  // allele frequencies
  log_allele_freqs = vector<vector<vector<double>>>(K, vector<vector<double>>(L));
  allele_counts = vector<vector<vector<int>>>(K, vector<vector<int>>(L));
  for (int k=0; k<K; k++) {
    for (int l=0; l<L; l++) {
      log_allele_freqs[k][l] = vector<double>(Jl[l]);
      allele_counts[k][l] = vector<int>(Jl[l]);
    }
  }
  
  // probabilities and likelihoods
  loglike = 0;
  
  // initialise Q matrix
  qmatrix_ind = vector<vector<double>>(n, vector<double>(K));
  log_qmatrix_ind = vector<vector<double>>(n, vector<double>(K));
  
  // initialise ordering of labels
  label_order = seq_int(0,K-1);
  label_order_new = vector<int>(K);
  
  // objects for solving label switching problem
  cost_mat = vector<vector<double>>(K, vector<double>(K));
  best_perm = vector<int>(K);
  best_perm_order = vector<int>(K);
  edges_left = vector<int>(K);
  edges_right = vector<int>(K);
  blocked_left = vector<int>(K);
  blocked_right = vector<int>(K);
  
}

//------------------------------------------------
// update group allocations
void particle_noadmix::update_group() {
  
  // loop through individuals
  int this_first_row = 0;
  for (int i=0; i<n; i++) {
    int this_ploidy = ploidy[i];
    
    // calculate log probability of being allocated to each group
    // TODO - reorder to only check missing once
    for (int k=0; k<K; k++) {
      log_qmatrix_ind[i][k] = 0;
      for (int l=0; l<L; l++) {
        for (int j=0; j<this_ploidy; j++) {
          int this_data = (*data_ptr)[this_first_row+j][l];
          if (this_data!=0) {
            log_qmatrix_ind[i][k] += beta_raised*log_allele_freqs[k][l][this_data-1];
          }
        }
      }
    }
    
    // exponentiate without underflow
    double max_log_qmatrix_ind = max(log_qmatrix_ind[i]);
    double qmatrix_ind_sum = 0;
    for (int k=0; k<K; k++) {
      log_qmatrix_ind[i][k] -= max_log_qmatrix_ind;
      qmatrix_ind[i][k] = exp(log_qmatrix_ind[i][k]);
      qmatrix_ind_sum += qmatrix_ind[i][k];
    }
    double log_inv_qmatrix_sum = -log(qmatrix_ind_sum);
    double inv_qmatrix_ind_sum = 1/qmatrix_ind_sum;
    for (int k=0; k<K; k++) {
      log_qmatrix_ind[i][k] += log_inv_qmatrix_sum;
      qmatrix_ind[i][k] *= inv_qmatrix_ind_sum;
    }
    
    // draw new group
    group[i] = sample1(qmatrix_ind[i]) - 1;
    
    this_first_row += this_ploidy;
  } // end loop over i
  
}

//------------------------------------------------
// update allele frequencies
void particle_noadmix::update_allele_freqs() {
  
  // reset allele counts
  for (int k=0; k<K; k++) {
    for (int l=0; l<L; l++) {
      fill(allele_counts[k][l].begin(), allele_counts[k][l].end(), 0);
    }
  }
  
  // recalculate allele counts
  int this_first_row = 0;
  for (int i=0; i<n; i++) {
    int this_ploidy = ploidy[i];
    int this_group = group[i];
    
    for (int l=0; l<L; l++) {
      for (int j=0; j<this_ploidy; j++) {
        int this_data = (*data_ptr)[this_first_row+j][l];
        if (this_data!=0) {
          allele_counts[this_group][l][this_data-1]++;
        }
      }
    }
    
    this_first_row += this_ploidy;
  } // end loop over i
  
  // draw allele freqs
  for (int k=0; k<K; k++) {
    for (int l=0; l<L; l++) {
      rdirichlet2(log_allele_freqs[k][l], allele_counts[k][l], beta_raised, lambda);
    }
  }
  
}

//------------------------------------------------
// fix label switching problem
void particle_noadmix::solve_label_switching(const vector<vector<double>> &log_qmatrix_ind_running) {
  
  foo();
  
  // fill in cost matrix
  for (int k1=0; k1<K; k1++) {
    fill(cost_mat[k1].begin(), cost_mat[k1].end(), 0);
    for (int k2=0; k2<K; k2++) {
      for (int i=0; i<n; i++) {
        //cost_mat[k1][k2] += qmatrix_ind[i][label_order[k1]]*(log_qmatrix_ind[i][label_order[k1]] - log_qmatrix_ind_running[i][k2]);
        cost_mat[k1][k2] += qmatrix_ind[i][k1]*(log_qmatrix_ind[i][k1] - log_qmatrix_ind_running[i][k2]);
      }
    }
  }
  
  // find best permutation of current labels using Hungarian algorithm
  best_perm = hungarian(cost_mat, edges_left, edges_right, blocked_left, blocked_right);
  
  // define best_perm_order
  for (int k=0; k<K; k++) {
    best_perm_order[best_perm[k]] = k;
  }
  
  // replace old label order with new
  for (int k=0; k<K; k++) {
    label_order_new[k] = label_order[best_perm_order[k]];
  }
  label_order = label_order_new;
  
}

//------------------------------------------------
// calculate overall log-likelihood
void particle_noadmix::calculate_loglike() {
  
  // Multinomial-Dirichlet likelihood
  loglike = 0;
  for (int k=0; k<K; k++) {
    for (int l=0; l<L; l++) {
      int allele_counts_total = 0;
      for (int j=0; j<Jl[l]; j++) {
        loglike += lgamma(lambda + allele_counts[k][l][j]) - lgamma(lambda);
        allele_counts_total += allele_counts[k][l][j];
      }
      loglike += lgamma(Jl[l]*lambda) - lgamma(Jl[l]*lambda + allele_counts_total);
    }
  }
  
}
