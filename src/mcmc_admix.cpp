
#include "mcmc_admix.h"
#include "misc.h"
#include "probability.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// constructor for MCMC class under the admixture model
mcmc_admix::mcmc_admix(Rcpp::List &args_data, Rcpp::List &args_model) {
  
  // extract data and parameters
  data = rcpp_to_mat_int(args_data["data"]);
  Jl = rcpp_to_vector_int(args_data["Jl"]);
  ploidy = rcpp_to_vector_int(args_data["ploidy"]);
  n = ploidy.size();
  L = Jl.size();
  
  K = rcpp_to_int(args_model["K"]);
  burnin = rcpp_to_int(args_model["burnin"]);
  samples = rcpp_to_int(args_model["samples"]);
  rungs = rcpp_to_int(args_model["rungs"]);
  GTI_pow = rcpp_to_double(args_model["GTI_pow"]);
  auto_converge = rcpp_to_bool(args_model["auto_converge"]);
  solve_label_switching_on = rcpp_to_bool(args_model["solve_label_switching_on"]);
  coupling_on = rcpp_to_bool(args_model["coupling_on"]);
  splitmerge_on = rcpp_to_bool(args_model["splitmerge_on"]);
  pb_markdown = rcpp_to_bool(args_model["pb_markdown"]);
  silent = rcpp_to_bool(args_model["silent"]);
  estimate_alpha = rcpp_to_bool(args_model["estimate_alpha"]);
  
  // thermodynamic parameters. The object beta_raised stores values of beta (the
  // thermodynamic power), raised to the power GTI_pow
  beta_raised = vector<double>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    beta_raised[rung] = pow((rung+1)/double(rungs), GTI_pow);
  }
  rung_order = seq_int(0,rungs-1);
  cold_rung = rung_order[rungs-1];
  
  // vector of particles
  particle_vec = vector<particle_admix>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    particle_vec[rung] = particle_admix(data, Jl, ploidy, args_model, beta_raised[rung]);
  }
  
  // Q-matrices
  log_qmatrix_gene_running = vector<vector<vector<double>>>(sum(ploidy), vector<vector<double>>(L, vector<double>(K)));
  qmatrix_gene = log_qmatrix_gene_running;
  qmatrix_ind = vector<vector<double>>(n, vector<double>(K));
  
  // initialise ordering of labels
  label_order = seq_int(0,K-1);
  
  // objects for solving label switching problem
  cost_mat = vector<vector<double>>(K, vector<double>(K));
  best_perm = vector<int>(K);
  best_perm_order = vector<int>(K);
  edges_left = vector<int>(K);
  edges_right = vector<int>(K);
  blocked_left = vector<int>(K);
  blocked_right = vector<int>(K);
  
  // objects for storing results
  loglike_burnin = vector<vector<double>>(rungs, vector<double>(burnin));
  loglike_sampling = vector<vector<double>>(rungs, vector<double>(samples));
  alpha_store = vector<double>(samples);
  
  // objects for storing acceptance rates
  coupling_accept = vector<int>(rungs-1);
  
}

//------------------------------------------------
// loop through burn-in iterations
void mcmc_admix::burnin_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress) {
  
  // skip if K==1
  if (K==1) {
    if (!silent) {
      print("Calculating exact solution for K = 1");
    }
    return;
  }
  
  // print header
  if (!silent) {
    print("Running MCMC for K =", K);
    print("Burn-in phase");
  }
  
  // read in R functions
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // define points at which convergence checked
  vector<int> convergence_checkpoint(10);
  for (int i=0; i<10; i++) {
    convergence_checkpoint[i] = (i+1)*double(burnin)/10;
  }
  int checkpoint_i = 0;
  
  // loop through burnin iterations
  bool convergence_reached = false;
  for (int rep=0; rep<burnin; rep++) {
    
    // update particles
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      
      // update group allocation of all individuals
      particle_vec[rung].update_group();
      
      // update allele and admixture frequencies
      particle_vec[rung].update_allele_admix_freqs();
      
      // update alpha
      if (estimate_alpha) {
        particle_vec[rung].update_alpha();
      }
      
      // calculate group log-likelihood
      particle_vec[rung].calculate_loglike();
      
    } // end loop over rungs
    
    // Metropolis-coupling
    if (coupling_on) {
      metropolis_coupling();
    }
    
    // focus on cold rung
    cold_rung = rung_order[rungs-1];
    
    // fix labels
    if (solve_label_switching_on) {
      //particle_vec[cold_rung].solve_label_switching(log_qmatrix_gene_running);
      solve_label_switching();
    }
    
    // add particle log_qmatrix to log_qmatrix_running
    update_log_qmatrix_gene_running();
    
    // store loglikelihood
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      loglike_burnin[r][rep] = particle_vec[rung].loglike;
    }
    
    // update progress bars
    if (!silent) {
      if ((rep+1)==burnin) {
        update_progress(args_progress, "pb_burnin", rep+1, burnin);
      } else {
        int remainder = rep % int(ceil(double(burnin)/100));
        if (remainder==0 && !pb_markdown) {
          update_progress(args_progress, "pb_burnin", rep+1, burnin);
        }
      }
    }
    
    // check for convergence
    if (auto_converge && (rep+1)==convergence_checkpoint[checkpoint_i]) {
      
      // break if convergence reached
      convergence_reached = true;
      for (int r=0; r<rungs; r++) {
        bool this_converged = rcpp_to_bool(test_convergence(loglike_burnin[r], rep+1));
        if (!this_converged) {
          convergence_reached = false;
          break;
        }
      }
      if (convergence_reached) {
        if (!silent) {
          update_progress(args_progress, "pb_burnin", burnin, burnin);
          print("   converged within", rep+1, "iterations");
        }
        for (int r=0; r<rungs; r++) {
          loglike_burnin[r].resize(rep+1);
        }
        break;
      }
      checkpoint_i++;
    }
    
  } // end burn-in iterations
  
  //print_array(log_qmatrix_gene_running);
  
  // warning if still not converged
  if (!convergence_reached && !silent) {
    print("   Warning: convergence still not reached within", burnin, "iterations");
  }
  
}

//------------------------------------------------
// loop through sampling iterations
void mcmc_admix::sampling_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress) {
  
  // skip if K==1
  if (K==1) {
    return;
  }
  
  // print header
  if (!silent) {
    print("Sampling phase");
  }
  
  // read in R functions
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // reset acceptance rates
  coupling_accept = vector<int>(rungs-1);
  
  // loop through sampling iterations
  for (int rep=0; rep<samples; rep++) {
    
    // update particles
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      
      // update group allocation of all gene copies
      particle_vec[rung].update_group();
      
      // update allele and admixture frequencies
      particle_vec[rung].update_allele_admix_freqs();
      
      // update alpha
      if (estimate_alpha) {
        particle_vec[rung].update_alpha();
      }
      
      // calculate group log-likelihood
      particle_vec[rung].calculate_loglike();
      
    } // end loop over rungs
    
    // Metropolis-coupling
    if (coupling_on) {
      metropolis_coupling();
    }
    
    // focus on cold rung
    cold_rung = rung_order[rungs-1];
    
    // fix labels
    if (solve_label_switching_on) {
      //particle_vec[cold_rung].solve_label_switching(log_qmatrix_gene_running);
      solve_label_switching();
    }
    
    // add particle log_qmatrix_gene to log_qmatrix_gene_running
    update_log_qmatrix_gene_running();
    
    // add particle qmatrix to qmatrix_gene
    update_qmatrix_gene();
    
    // store loglikelihood
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      loglike_sampling[r][rep] = particle_vec[rung].loglike;
    }
    
    // store alpha
    alpha_store[rep] = particle_vec[cold_rung].alpha;
    
    // update progress bars
    if (!silent) {
      if ((rep+1)==samples) {
        update_progress(args_progress, "pb_samples", rep+1, samples);
      } else {
        int remainder = rep % int(ceil(double(samples)/100));
        if (remainder==0 && !pb_markdown) {
          update_progress(args_progress, "pb_samples", rep+1, samples);
        }
      }
    }
    
  } // end sampling iterations
  
  // finalise qmatrix_gene and qmatrix_ind
  double inv_samples = 1.0/double(samples);
  int this_first_row = 0;
  for (int i=0; i<n; i++) {
    int this_ploidy = ploidy[i];
    for (int l=0; l<L; l++) {
      for (int j=0; j<this_ploidy; j++) {
        for (int k=0; k<K; k++) {
          qmatrix_gene[this_first_row+j][l][k] *= inv_samples;
          qmatrix_ind[i][k] += qmatrix_gene[this_first_row+j][l][k];
        }
      }
    }
    for (int k=0; k<K; k++) {
      qmatrix_ind[i][k] /= double(this_ploidy*L);
    }
    this_first_row += this_ploidy;
  }
  
}

//------------------------------------------------
// update running estimate of qmatrix
void mcmc_admix::update_log_qmatrix_gene_running() {
  
  int this_first_row = 0;
  for (int i=0; i<n; i++) {
    int this_ploidy = ploidy[i];
    for (int l=0; l<L; l++) {
      for (int j=0; j<this_ploidy; j++) {
        for (int k=0; k<K; k++) {
          //int tmp = particle_vec[cold_rung].log_qmatrix_gene[this_first_row+j][l][k];
          int this_k = label_order[k];
          log_qmatrix_gene_running[this_first_row+j][l][k] = log_sum(log_qmatrix_gene_running[this_first_row+j][l][k], particle_vec[cold_rung].log_qmatrix_gene[this_first_row+j][l][this_k]);
        }
      }
    }
    this_first_row += this_ploidy;
  }
  
}

//------------------------------------------------
// add qmatrix to qmatrix_gene
void mcmc_admix::update_qmatrix_gene() {
  
  int this_first_row = 0;
  for (int i=0; i<n; i++) {
    int this_ploidy = ploidy[i];
    for (int l=0; l<L; l++) {
      for (int j=0; j<this_ploidy; j++) {
        for (int k=0; k<K; k++) {
          int this_k = label_order[k];
          qmatrix_gene[this_first_row+j][l][k] += particle_vec[cold_rung].qmatrix_gene[this_first_row+j][l][this_k];
        }
      }
    }
    this_first_row += this_ploidy;
  }
  
}

//------------------------------------------------
// Metropolis-coupling to propose swaps between temperature rungs
void mcmc_admix::metropolis_coupling() {
  
  // loop over rungs, starting with the hottest chain and moving to the cold
  // chain. Each time propose a swap with the next rung up
  for (int i=0; i<(rungs-1); i++) {
    
    // define rungs of interest
    int rung1 = rung_order[i];
    int rung2 = rung_order[i+1];
    
    // get log-likelihoods and beta values of two chains in the comparison
    double loglike1 = particle_vec[rung1].loglike;
    double loglike2 = particle_vec[rung2].loglike;
    
    double beta_raised1 = particle_vec[rung1].beta_raised;
    double beta_raised2 = particle_vec[rung2].beta_raised;
    
    // calculate acceptance ratio (still in log space)
    double acceptance = (loglike2*beta_raised1 + loglike1*beta_raised2) - (loglike1*beta_raised1 + loglike2*beta_raised2);
    
    // accept or reject move
    if (log(runif1())<acceptance) {
      
      // swap beta values
      particle_vec[rung1].beta_raised = beta_raised2;
      particle_vec[rung2].beta_raised = beta_raised1;
      
      // swap rung order
      rung_order[i] = rung2;
      rung_order[i+1] = rung1;
      
      // update acceptance rates
      coupling_accept[i]++;
    }
  }
  
}

//------------------------------------------------
// fix label switching problem
void mcmc_admix::solve_label_switching() {
  
  // fill in cost matrix
  for (int k1=0; k1<K; k1++) {
    fill(cost_mat[k1].begin(), cost_mat[k1].end(), 0);
    for (int k2=0; k2<K; k2++) {
      int this_first_row = 0;
      for (int i=0; i<n; i++) {
        int this_ploidy = ploidy[i];
        for (int l=0; l<L; l++) {
          for (int j=0; j<this_ploidy; j++) {
            cost_mat[k1][k2] += particle_vec[cold_rung].qmatrix_gene[this_first_row+j][l][k1]*(particle_vec[cold_rung].log_qmatrix_gene[this_first_row+j][l][k1] - log_qmatrix_gene_running[this_first_row+j][l][k2]);
            //cost_mat[k1][k2] += particle_vec[cold_rung].qmatrix_gene[this_first_row+j][l][label_order[k1]]*(particle_vec[cold_rung].log_qmatrix_gene[this_first_row+j][l][label_order[k1]] - log_qmatrix_gene_running[this_first_row+j][l][k2]);
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
    label_order[best_perm[k]] = k;
  }
  
}
