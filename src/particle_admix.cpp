
#include "particle_admix.h"
#include "misc.h"
#include "probability.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// default constructor for particle class under admixture model
particle_admix::particle_admix() {}

//------------------------------------------------
// constructor for particle class under no-admixture model
particle_admix::particle_admix(vector<vector<int>> &data, vector<int> &Jl_, vector<int> &ploidy_, Rcpp::List &args_model, double beta_raised_) {
  
  // extract data and parameters
  data_ptr = &data;
  Jl = Jl_;
  ploidy = ploidy_;
  n = ploidy.size();
  L = Jl.size();
  
  beta_raised = beta_raised_;
  K = rcpp_to_int(args_model["K"]);
  lambda = rcpp_to_double(args_model["lambda"]);
  alpha = rcpp_to_double(args_model["alpha"]);
  
  // initialise group allocation
  group = vector<vector<int>>(sum(ploidy), vector<int>(L));
  
  // allele frequencies
  log_allele_freqs = vector<vector<vector<double>>>(K, vector<vector<double>>(L));
  allele_counts = vector<vector<vector<int>>>(K, vector<vector<int>>(L));
  for (int k=0; k<K; k++) {
    for (int l=0; l<L; l++) {
      log_allele_freqs[k][l] = vector<double>(Jl[l]);
      allele_counts[k][l] = vector<int>(Jl[l]);
    }
  }
  
  // admixture proportions
  log_admix_freqs = vector<vector<double>>(n, vector<double>(K));
  admix_counts = vector<vector<int>>(n, vector<int>(K));
  admix_counts_totals = vector<int>(n);
  int this_first_row = 0;
  for (int i=0; i<n; i++) {
    int this_ploidy = ploidy[i];
    for (int l=0; l<L; l++) {
      for (int j=0; j<this_ploidy; j++) {
        int this_data = (*data_ptr)[this_first_row+j][l];
        if (this_data!=0) {
          admix_counts_totals[i]++;
        }
      }
    }
  }
  
  // individual-level group update objects
  log_admix_freqs_prop = vector<double>(K);
  log_qmatrix_prop = vector<vector<vector<double>>>(max(ploidy), vector<vector<double>>(L, vector<double>(K)));
  qmatrix_prop = log_qmatrix_prop;
  group_prop = vector<vector<int>>(max(ploidy), vector<int>(L));
  
  // probabilities and likelihoods
  loglike = 0;
  
  // initialise Q matrices
  log_qmatrix_gene = vector<vector<vector<double>>>(sum(ploidy), vector<vector<double>>(L, vector<double>(K)));
  qmatrix_gene = log_qmatrix_gene;
  
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
void particle_admix::update_group() {
  
  // loop through individuals
  int this_first_row = 0;
  for (int i=0; i<n; i++) {
    int this_ploidy = ploidy[i];
    
    // allocate gene copies to groups by drawing from conditional posterior. 
    // Simultaneously calculate the log likelihood and log proposal probability 
    // of this group allocation.
    double loglike_old = 0;
    double log_prop_old = 0;
    for (int l=0; l<L; l++) {
      for (int j=0; j<this_ploidy; j++) {
        
        // calculate log probability of being allocated to each group
        int this_data = (*data_ptr)[this_first_row+j][l];
        if (this_data==0) {
          continue;
        }
        for (int k=0; k<K; k++) {
          log_qmatrix_gene[this_first_row+j][l][k] = log_admix_freqs[i][k] + beta_raised*log_allele_freqs[k][l][this_data-1];
        }
        
        // exponentiate without underflow
        double max_log_qmatrix = max(log_qmatrix_gene[this_first_row+j][l]);
        double qmatrix_sum = 0;
        for (int k=0; k<K; k++) {
          log_qmatrix_gene[this_first_row+j][l][k] -= max_log_qmatrix;
          qmatrix_gene[this_first_row+j][l][k] = exp(log_qmatrix_gene[this_first_row+j][l][k]);
          qmatrix_sum += qmatrix_gene[this_first_row+j][l][k];
        }
        double log_inv_qmatrix_sum = -log(qmatrix_sum);
        double inv_qmatrix_sum = 1/qmatrix_sum;
        for (int k=0; k<K; k++) {
          log_qmatrix_gene[this_first_row+j][l][k] += log_inv_qmatrix_sum;
          qmatrix_gene[this_first_row+j][l][k] *= inv_qmatrix_sum;
        }
        
        // draw new group
        int new_group = sample1(qmatrix_gene[this_first_row+j][l]) - 1;
        group[this_first_row+j][l] = new_group;
        
        // calculate log likelihood and proposal probability
        log_prop_old += log_qmatrix_gene[this_first_row+j][l][new_group];
        loglike_old += log_admix_freqs[i][new_group] + log_allele_freqs[new_group][l][this_data-1];
        
      } // end j loop
    } // end l loop
    
    // propose new admixture frequencies by drawing from prior
    rdirichlet2(log_admix_freqs_prop, admix_counts[i], 0, alpha);
    
    // propose new grouping from new admixture frequencies
    double loglike_new = 0;
    double log_prop_new = 0;
    for (int l=0; l<L; l++) {
      for (int j=0; j<this_ploidy; j++) {
        
        // calculate log probability of being allocated to each group
        int this_data = (*data_ptr)[this_first_row+j][l];
        if (this_data==0) {
          continue;
        }
        for (int k=0; k<K; k++) {
          log_qmatrix_prop[j][l][k] = log_admix_freqs_prop[k] + beta_raised*log_allele_freqs[k][l][this_data-1];
        }
        
        // exponentiate without underflow
        double max_log_qmatrix = max(log_qmatrix_prop[j][l]);
        double qmatrix_sum = 0;
        for (int k=0; k<K; k++) {
          log_qmatrix_prop[j][l][k] -= max_log_qmatrix;
          qmatrix_prop[j][l][k] = exp(log_qmatrix_prop[j][l][k]);
          qmatrix_sum += qmatrix_prop[j][l][k];
        }
        
        double log_inv_qmatrix_sum = -log(qmatrix_sum);
        double inv_qmatrix_sum = 1/qmatrix_sum;
        for (int k=0; k<K; k++) {
          log_qmatrix_prop[j][l][k] += log_inv_qmatrix_sum;
          qmatrix_prop[j][l][k] *= inv_qmatrix_sum;
        }
        
        // draw new group
        int new_group = sample1(qmatrix_prop[j][l]) - 1;
        group_prop[j][l] = new_group;
        log_prop_new += log_qmatrix_prop[j][l][new_group];
        loglike_new += log_admix_freqs_prop[new_group] + log_allele_freqs[new_group][l][this_data-1];
        
      } // end j loop
    } // end l loop
    
    // Metropolis-Hastings step
    double MH = (loglike_new - loglike_old) - (log_prop_new - log_prop_old);
    if (log(runif_0_1()) < MH) {
      for (int l=0; l<L; l++) {
        for (int j=0; j<this_ploidy; j++) {
          group[this_first_row+j][l] = group_prop[j][l];
          log_qmatrix_gene[this_first_row+j][l] = log_qmatrix_prop[j][l];
          qmatrix_gene[this_first_row+j][l] = qmatrix_prop[j][l];
        }
      }
      log_admix_freqs[i] = log_admix_freqs_prop;
    }
    
    this_first_row += this_ploidy;
  } // end i loop
  
}

//------------------------------------------------
// update allele and admixture frequencies
void particle_admix::update_allele_admix_freqs() {
  
  // reset allele counts
  for (int k=0; k<K; k++) {
    for (int l=0; l<L; l++) {
      fill(allele_counts[k][l].begin(), allele_counts[k][l].end(), 0);
    }
  }
  
  // reset admix counts
  for (int i=0; i<n; i++) {
    fill(admix_counts[i].begin(), admix_counts[i].end(), 0);
  }
  
  // recalculate allele and admix counts
  int this_first_row = 0;
  for (int i=0; i<n; i++) {
    int this_ploidy = ploidy[i];
    
    for (int l=0; l<L; l++) {
      for (int j=0; j<this_ploidy; j++) {
        int this_group = group[this_first_row+j][l];
        int this_data = (*data_ptr)[this_first_row+j][l];
        if (this_data!=0) {
          allele_counts[this_group][l][this_data-1]++;
          admix_counts[i][this_group]++;
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
  
  // draw admix freqs
  for (int i=0; i<n; i++) {
    rdirichlet2(log_admix_freqs[i], admix_counts[i], 1.0, alpha);
  }
  
}

//------------------------------------------------
// update admixture parameter alpha
void particle_admix::update_alpha2() {
  
  double alpha_prop = rnorm1(alpha, 0.5);
  if (alpha_prop<0) {
    alpha_prop = -alpha_prop;
  }
  
  double ll_old = 0;
  double ll_new = 0;
  for (int i=0; i<n; i++) {
    ll_old += lgamma(K*alpha);
    ll_new += lgamma(K*alpha_prop);
    for (int k=0; k<K; k++) {
      ll_old += (alpha-1)*log_admix_freqs[i][k] - lgamma(alpha);
      ll_new += (alpha_prop-1)*log_admix_freqs[i][k] - lgamma(alpha_prop);
    }
  }
  
  if (log(runif_0_1()) < (ll_new-ll_old)) {
    alpha = alpha_prop;
  }
  
}

//------------------------------------------------
// update admixture parameter alpha
void particle_admix::update_alpha() {
  
  // prior parameters on alpha (shape and rate of gamma prior)
  double alpha_shape = 1.0;
  double alpha_rate = 0.2;
  
  // draw latent variables
  int z = 0;
  double log_phi = 0;
  for (int i=0; i<n; i++) {
    log_phi += log(rbeta1(K*alpha+1, admix_counts_totals[i]));
    for (int k=0; k<K; k++) {
      z += rCRPgroups(admix_counts[i][k], alpha);
    }
  }
  
  // draw final gamma parameters
  double w1 = alpha_shape + z - n;
  double w2 = alpha_rate - K*log_phi;
  double p;
  for (int i=0; i<n; i++) {
    p = admix_counts_totals[i]/(admix_counts_totals[i] + K*w1/w2);
    w1 += rbernoulli1(1-p);
  }
  
  // draw new alpha
  alpha = rgamma1(w1, w2);
  
}

//------------------------------------------------
// fix label switching problem
void particle_admix::solve_label_switching(const vector<vector<vector<double>>> &log_qmatrix_gene_running) {
  
  //log_qmatrix_gene = vector<vector<vector<double>>>(sum(ploidy), vector<vector<double>>(L, vector<double>(K)));
  
  // fill in cost matrix
  for (int k1=0; k1<K; k1++) {
    fill(cost_mat[k1].begin(), cost_mat[k1].end(), 0);
    for (int k2=0; k2<K; k2++) {
      int this_first_row = 0;
      for (int i=0; i<n; i++) {
        int this_ploidy = ploidy[i];
        for (int l=0; l<L; l++) {
          for (int j=0; j<this_ploidy; j++) {
            cost_mat[k1][k2] += qmatrix_gene[this_first_row+j][l][label_order[k1]]*(log_qmatrix_gene[this_first_row+j][l][label_order[k1]] - log_qmatrix_gene_running[this_first_row+j][l][k2]);
          }
        }
        this_first_row += this_ploidy;
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
void particle_admix::calculate_loglike() {
  
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
