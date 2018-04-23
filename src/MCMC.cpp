
#include "MCMC.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// constructor for MCMC class
MCMC::MCMC(vector<double> &x, Rcpp::List &args) {
  
  // extract data and parameters
  x_ptr = &x;
  n = x.size();
  K = rcpp_to_int(args["K"]);
  burnin = rcpp_to_int(args["burnin"]);
  samples = rcpp_to_int(args["samples"]);
  iterations = burnin + samples;
  rungs = rcpp_to_int(args["rungs"]);
  solve_label_switching_on = rcpp_to_bool(args["solve_label_switching_on"]);
  coupling_on = rcpp_to_bool(args["coupling_on"]);
  scaf_on = rcpp_to_bool(args["scaffold_on"]);
  scaf_n = rcpp_to_int(args["scaf_n"]);
  splitmerge_on = rcpp_to_bool(args["splitmerge_on"]);
  parallel_on = rcpp_to_bool(args["parallel_on"]);
  print_console = !parallel_on;
  
  // thermodynamic parameters
  beta_vec = vector<double>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    beta_vec[rung] = (rungs==1) ? 1 : rung/double(rungs-1);
  }
  rung_order = seq_int(0,rungs-1);
  cold_rung = rung_order[rungs-1];
  
  // vector of particles
  particle_vec = std::vector<particle>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    particle_vec[rung] = particle(x, args, beta_vec[rung]);
  }
  
  // scaffold objects
  scaf_group = vector<vector<vector<int>>>(rungs, vector<vector<int>>(scaf_n, vector<int>(n)));
  scaf_counts = vector<vector<vector<int>>>(rungs, vector<vector<int>>(scaf_n, vector<int>(K)));
  scaf_x_sum = vector<vector<vector<double>>>(rungs, vector<vector<double>>(scaf_n, vector<double>(K)));
  
  // objects for storing results
  mu_store = vector<vector<double>>(iterations, vector<double>(K));
  loglike_store = vector<vector<double>>(rungs, vector<double>(iterations));
  qmatrix_running = vector<vector<double>>(n, vector<double>(K));
  log_qmatrix_running = vector<vector<double>>(n, vector<double>(K));
  qmatrix_final = vector<vector<double>>(n, vector<double>(K));
  
  // objects for storing acceptance rates
  mc_accept = vector<int>(rungs-1);
  scaf_accept = vector<int>(rungs);
  splitmerge_accept = vector<int>(rungs);
}

//------------------------------------------------
// generate scaffold groupings
void MCMC::scaffold_mcmc(Rcpp::List &args) {
  
  // print header
  if (print_console) {
    print("Generating", scaf_n, "scaffolds");
  }
  
  // extract R functions
  Rcpp::Function test_convergence = args["test_convergence"];
  Rcpp::Function update_progress = args["update_progress"];
  
  // initialise objects needed for generating scaffolds
  int batch_size = 10;
  int max_batches = 10;
  vector<double> batch_vec(batch_size);
  int dummy_accept = 0;
  
  // generate multiple scaffolds
  for (int scaf_rep=0; scaf_rep<scaf_n; scaf_rep++) {
    
    // reset particles
    for (int r=0; r<rungs; r++) {
      particle_vec[r].reset();
      particle_vec[r].beta = beta_vec[r];
    }
    rung_order = seq_int(0,rungs-1);
    
    // store log-likelihoods
    vector<vector<double>> scaf_log_like(rungs, vector<double>(batch_size));
    
    // loop through iterations in batches
    for (int batch=0; batch<max_batches; batch++) {
      
      // iterations of this batch
      for (int rep=0; rep<batch_size; rep++) {
        
        // update particles
        for (int r=0; r<rungs; r++) {
          int rung = rung_order[r];
          
          // update parameters
          particle_vec[rung].update_mu();
          particle_vec[rung].update_group();
          
          // calculate log-likelihood
          particle_vec[rung].calc_log_like();
          
          // split-merge step
          if (splitmerge_on) {
            particle_vec[rung].splitmerge_propose(dummy_accept);
          }
        }
        
        // Metropolis-coupling
        if (coupling_on) {
          metropolis_coupling();
        }
        
        // store log-likelihoods
        for (int r=0; r<rungs; r++) {
          int rung = rung_order[r];
          scaf_log_like[rung][batch*batch_size + rep] = particle_vec[rung].loglike;
        }
        
      } // end iterations of this batch
      
      // break if converged
      bool all_converged = true;
      for (int r=0; r<rungs; r++) {
        int rung = rung_order[r];
        bool rung_converged = test_convergence(scaf_log_like[rung]);
        if (!rung_converged) {
          all_converged = false;
          break;
        }
      }
      if (all_converged) {
        break;
      }
      
      // if not converged expand scaf_log_like
      for (int r=0; r<rungs; r++) {
        push_back_multiple(scaf_log_like[r], batch_vec);
      }
      
      // throw warning if still not converged at end of batches
      if (batch == (max_batches-1)) {
        print("warning: scaffold did not converge");
      }
      
    } // loop over batches
    
    // loop over rungs
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      
      // re-order group allocation to be always-increasing
      particle_vec[rung].group_increasing();
      
      // store this particle
      scaf_group[rung][scaf_rep] = particle_vec[rung].group;
      for (int i=0; i<n; i++) {
        scaf_counts[rung][scaf_rep][scaf_group[rung][scaf_rep][i]] ++;
        scaf_x_sum[rung][scaf_rep][scaf_group[rung][scaf_rep][i]] += (*x_ptr)[i];
      }
    }
    
    // update scaffold progress bar
    if (print_console) {
      update_progress(args, 1, scaf_rep+1, scaf_n);
    }
    
  } // loop over scaf_n
  
  // load these scaffolds back into particles
  for (int r=0; r<rungs; r++) {
    particle_vec[r].scaf_group = scaf_group[r];
    particle_vec[r].scaf_counts = scaf_counts[r];
    particle_vec[r].scaf_x_sum = scaf_x_sum[r];
  }
  
}

//------------------------------------------------
// run MCMC
void MCMC::main_mcmc(Rcpp::List &args) {
  
  // extract R functions
  Rcpp::Function update_progress = args["update_progress"];
  
  // loop through burn-in iterations
  for (int rep=0; rep<iterations; rep++) {
    
    // print headers
    if (rep==0 && print_console) {
      print("MCMC burn-in phase");
    }
    if (rep==burnin && print_console) {
      print("MCMC sampling phase");
    }
    
    // update particles
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      
      // update parameters
      particle_vec[rung].update_mu();
      particle_vec[rung].update_group();
      
      // calculate log-likelihood
      particle_vec[rung].calc_log_like();
      
      // propose swap with scaffolds
      if (scaf_on) {
        particle_vec[rung].scaf_propose(scaf_accept[r]);
      }
      
      // split-merge step
      if (splitmerge_on) {
        particle_vec[rung].splitmerge_propose(splitmerge_accept[r]);
      }
    }
    
    // Metropolis-coupling
    if (coupling_on) {
      metropolis_coupling();
    }
    
    // focus on cold rung
    cold_rung = rung_order[rungs-1];
    
    // fix labels
    if (solve_label_switching_on) {
      particle_vec[cold_rung].solve_label_switching(log_qmatrix_running);
    }
    
    // add qmatrix to qmatrix_running and qmatrix_final
    update_qmatrix_running();
    if (rep >= burnin) {
      update_qmatrix_final();
    }
    
    // store results of this iteration
    for (int k=0; k<K; k++) {
      mu_store[rep][k] = particle_vec[cold_rung].mu[particle_vec[cold_rung].label_order[k]];
    }
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      loglike_store[r][rep] = particle_vec[rung].loglike;
    }
    
    // update progress bars
    if (print_console) {
      if (rep<burnin) {
        int remainder = rep % int(ceil(burnin/100));
        if (remainder==0 || (rep+1)==burnin) {
          update_progress(args, 2, rep+1, burnin);
        }
      } else {
        int remainder = (rep-burnin) % int(ceil(samples/100));
        if (remainder==0 || (rep-burnin+1)==samples) {
          update_progress(args, 3, rep-burnin+1, samples);
        }
      }
    }
    
  } // end MCMC loop
  
  // finalise Q-matrix
  for (int i=0; i<n; i++) {
    for (int j=0; j<K; j++) {
      qmatrix_final[i][j] /= samples;
    }
  }
  
}

//------------------------------------------------
// add qmatrix to qmatrix_running
void MCMC::update_qmatrix_running() {
  
  for (int i=0; i<n; i++) {
    for (int k=0; k<K; k++) {
      qmatrix_running[i][k] += particle_vec[cold_rung].qmatrix[i][particle_vec[cold_rung].label_order[k]];
      double lq = log(qmatrix_running[i][k]);
      if (lq < -OVERFLO) {
        log_qmatrix_running[i][k] = -OVERFLO;
      } else {
        log_qmatrix_running[i][k] = lq;
      }
    }
  }
  
}

//------------------------------------------------
// add qmatrix to qmatrix_final
void MCMC::update_qmatrix_final() {
  
  for (int i=0; i<n; i++) {
    for (int k=0; k<K; k++) {
      qmatrix_final[i][k] += particle_vec[cold_rung].qmatrix[i][particle_vec[cold_rung].label_order[k]];
    }
  }
  
}

//------------------------------------------------
// Metropolis-coupling to propose swaps between temperature rungs
void MCMC::metropolis_coupling() {
  
  // loop over rungs, starting with the hottest chain and moving to the cold
  // chain. Each time propose a swap with the next rung up
  for (int i=0; i<(rungs-1); i++) {
    
    // define rungs of interest
    int rung1 = rung_order[i];
    int rung2 = rung_order[i+1];
    
    // get log-likelihoods and beta values of two chains in the comparison
    double log_like1 = particle_vec[rung1].loglike;
    double log_like2 = particle_vec[rung2].loglike;
    
    double beta1 = particle_vec[rung1].beta;
    double beta2 = particle_vec[rung2].beta;
    
    // calculate acceptance ratio (still in log space)
    double acceptance = (log_like2*beta1 + log_like1*beta2) - (log_like1*beta1 + log_like2*beta2);
    
    // accept or reject move
    double rand1 = runif1();
    if (log(rand1)<acceptance) {
      
      // swap beta values
      particle_vec[rung1].beta = beta2;
      particle_vec[rung2].beta = beta1;
      
      // swap rung order
      rung_order[i] = rung2;
      rung_order[i+1] = rung1;
      
      // update acceptance rates
      mc_accept[i] ++;
    }
  }
}
