
#include <chrono>
#include "main.h"
#include "misc.h"
#include "probability.h"
#include "mcmc_noadmix.h"
#include "mcmc_admix.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// run main MCMC
// [[Rcpp::export]]
Rcpp::List run_mcmc_cpp(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_data = args["args_data"];
  Rcpp::List args_model = args["args_model"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  
  // extract arguments
  bool silent = rcpp_to_bool(args_model["silent"]);
  bool admix_on = rcpp_to_bool(args_model["admix_on"]);
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // define objects for storing results
  vector<vector<double>> qmatrix_ind;
  vector<vector<double>> loglike_burnin;
  vector<vector<double>> loglike_sampling;
  vector<double> alpha_store;
  vector<int> coupling_accept;
  vector<bool> rung_converged;
  
  // run MCMC and store results
  if (admix_on) { // admixture model
    
    // initialise MCMC object
    mcmc_admix m(args_data, args_model);
    
    // find starting group
    m.starting_group();
    
    // run MCMC
    m.burnin_mcmc(args_functions, args_progress);
    m.sampling_mcmc(args_functions, args_progress);
    
    // store results
    qmatrix_ind = get_qmatrix_ind(m);
    loglike_burnin = get_loglike_burnin(m);
    loglike_sampling = get_loglike_sampling(m);
    alpha_store = get_alpha_store(m);
    coupling_accept = get_coupling_accept(m);
    rung_converged = get_rung_converged(m);
    
  } else { // no-admixture model
    
    // initialise MCMC object
    mcmc_noadmix m(args_data, args_model);
    
    // find starting group
    m.starting_group();
    
    // run MCMC
    m.burnin_mcmc(args_functions, args_progress);
    m.sampling_mcmc(args_functions, args_progress);
    
    // store results
    qmatrix_ind = get_qmatrix_ind(m);
    loglike_burnin = get_loglike_burnin(m);
    loglike_sampling = get_loglike_sampling(m);
    coupling_accept = get_coupling_accept(m);
    rung_converged = get_rung_converged(m);
  }
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  if (!silent) {
    print("   completed in", time_span.count(), "seconds\n");
  }
  
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( loglike_burnin ));
  ret.push_back(Rcpp::wrap( loglike_sampling ));
  ret.push_back(Rcpp::wrap( qmatrix_ind ));
  ret.push_back(Rcpp::wrap( alpha_store ));
  ret.push_back(Rcpp::wrap( coupling_accept ));
  ret.push_back(Rcpp::wrap( rung_converged ));
  ret.push_back(Rcpp::wrap( time_span.count() ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("loglike_burnin");
  ret_names.push_back("loglike_sampling");
  ret_names.push_back("qmatrix_ind");
  ret_names.push_back("alpha_store");
  ret_names.push_back("coupling_accept");
  ret_names.push_back("rung_converged");
  ret_names.push_back("run_time");
  
  ret.names() = ret_names;
  return ret;
}

//------------------------------------------------
// estimate quantiles of posterior probability of K by simulation
// [[Rcpp::export]]
Rcpp::List GTI_posterior_K_sim_cpp(Rcpp::List args) {
  
  // extract arguments
  vector<double> m = rcpp_to_vector_double(args["mean"]);
  vector<double> s = rcpp_to_vector_double(args["SE"]);
  int reps = rcpp_to_int(args["reps"]);
  int K = m.size();
  
  // obtain normalised draws
  vector<double> y(K);
  double y_max = 0;
  vector<double> y_trans(K);
  double y_trans_sum = 0;
  double y_trans_inv_sum = 0;
  vector<vector<double>> ret(K, vector<double>(reps));
  for (int i=0; i<reps; i++) {
    for (int k=0; k<K; k++) {
      y[k] = rnorm1(m[k], s[k]);
      if (k==0 || y[k]>y_max) {
        y_max = y[k];
      }
    }
    y_trans_sum = 0;
    for (int k=0; k<K; k++) {
      y_trans[k] = exp(y[k]-y_max);
      y_trans_sum += y_trans[k];
    }
    y_trans_inv_sum = 1/y_trans_sum;
    for (int k=0; k<K; k++) {
      ret[k][i] = y_trans[k]*y_trans_inv_sum;
    }
  }
  
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("ret")=ret);
}

//------------------------------------------------
// integrate log-evidence over K by simulation
// [[Rcpp::export]]
Rcpp::List GTI_integrated_K_sim_cpp(Rcpp::List args) {
  
  // extract arguments
  vector<double> m = rcpp_to_vector_double(args["mean"]);
  vector<double> s = rcpp_to_vector_double(args["SE"]);
  int reps = rcpp_to_int(args["reps"]);
  int K = m.size();
  
  // obtain integrated draws
  vector<double> y(K);
  double y_max = 0;
  double y_trans_sum = 0;
  double x = 0;
  double x_sum = 0;
  double x_sum_squared = 0;
  for (int i=0; i<reps; i++) {
    for (int k=0; k<K; k++) {
      y[k] = rnorm1(m[k], s[k]);
      if (k==0 || y[k]>y_max) {
        y_max = y[k];
      }
    }
    y_trans_sum = 0;
    for (int k=0; k<K; k++) {
      y_trans_sum += exp(y[k]-y_max);
    }
    x = y_max + log(y_trans_sum) - log(K);
    x_sum += x;
    x_sum_squared += x*x;
  }
  double x_mean = x_sum/double(reps);
  double x_var = x_sum_squared/double(reps) - x_mean*x_mean;
  if (x_var<0) {
    x_var = 0;
  }
  
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("mean")=x_mean,
                            Rcpp::Named("SE")=x_var);
}