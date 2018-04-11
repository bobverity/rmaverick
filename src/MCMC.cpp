
#include "MCMC.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// constructor for MCMC class
MCMC::MCMC(Rcpp::NumericVector x_, Rcpp::List args_params) {
  
  // extract data and parameters
  x = x_;
  n = x.length();
  K = rcpp_to_int(args_params["K"]);
  mu_prior_mean = rcpp_to_double(args_params["mu_prior_mean"]);
  mu_prior_var = rcpp_to_double(args_params["mu_prior_var"]);
  burnin = rcpp_to_int(args_params["burnin"]);
  samples = rcpp_to_int(args_params["samples"]);
  iterations = burnin + samples;
  rungs = rcpp_to_int(args_params["rungs"]);
  mc_interval = rcpp_to_int(args_params["mc_interval"]);
  num_threads = rcpp_to_int(args_params["num_threads"]);
  
  // thermodynamic powers
  beta_vec = Rcpp::NumericVector(rungs);
  for (int i=0; i<rungs; i++) {
    beta_vec(i) = i/double(rungs-1);
  }
  
  // initialize component means
  mu = Rcpp::NumericMatrix(rungs, K);
  
  // initialize group allocation
  group = Rcpp::IntegerMatrix(rungs, n);
  for (int i=0; i<rungs; i++) {
    for (int j=0; j<n; j++) {
      group(i, j) = sample2(0, K-1);
    }
  }
  
  // define objects for storing results
  mu_store = Rcpp::NumericMatrix(rungs, K*iterations);
  
}

//------------------------------------------------
// run serial version of MCMC
void MCMC::serial_mcmc() {
  
  print("serial MCMC");
  
}

//------------------------------------------------
// functor to run multiple MCMC chains
struct functor_run_chains : public RcppParallel::Worker {
  
  // input objects
  RcppParallel::RVector<int> start_it;
  RcppParallel::RVector<int> end_it;
  RcppParallel::RMatrix<int> group;
  RcppParallel::RMatrix<double> mu;
  const double mu_prior_mean;
  const double mu_prior_var;
  const int K;
  RcppParallel::RVector<double> beta;
  
  // output objects
  RcppParallel::RMatrix<double> mu_store;
  
  // initialize with input and output
  functor_run_chains(Rcpp::IntegerVector start_it, Rcpp::IntegerVector end_it, const Rcpp::IntegerMatrix group, const Rcpp::NumericMatrix mu, const double mu_prior_mean, const double mu_prior_var, const int K, const Rcpp::NumericVector beta, Rcpp::NumericMatrix mu_store)
    : start_it(start_it), end_it(end_it), group(group), mu(mu), mu_prior_mean(mu_prior_mean), mu_prior_var(mu_prior_var), K(K), beta(beta), mu_store(mu_store) {}
  
  // operator to be called by parallel functions
  void operator()(size_t begin, size_t end) {
    
    // loop over MCMC chains
    for (size_t rung=begin; rung<end; rung++) {
      
      // MCMC for this chunk
      for (int rep=start_it[rung]; rep<end_it[rung]; rep++) {
        
        // update component means
        update_mu(mu_prior_mean, rung*mu_prior_var, mu, rung, K);
        
        // store results
        for (int k=0; k<K; k++) {
          mu_store(rung, rep*K + k) = mu(rung, k);
        }
      }
      
      // move chunk forward
      start_it[rung] = end_it[rung];
      
    }
    
  }
};

//------------------------------------------------
// Run parallel version of MCMC
void MCMC::parallel_mcmc() {
  
  // print header
  print("Running", num_threads, "threads in parallel");
  
  // store start and end iterations for each parallel chunk
  Rcpp::IntegerVector start_it(rungs);
  Rcpp::IntegerVector end_it(rungs);
  fill(end_it.begin(), end_it.end(), mc_interval);
  
  // pass input and output to functor
  functor_run_chains run_chains(start_it, end_it, group, mu, mu_prior_mean, mu_prior_var, K, beta_vec, mu_store);
  
  // carry out MCMC in chunks
  for (int i=0; i<(iterations/mc_interval); i++) {
    
    //print("chunk", i);
    
    // run this chunk in parallel
    parallelFor(0, rungs, run_chains);
    
    // move to next chunk
    fill(start_it.begin(), start_it.end(), end_it[0]);
    fill(end_it.begin(), end_it.end(), end_it[0] + mc_interval);
    
    rcpp_print_vector(start_it);
  }
  
}

//------------------------------------------------
// update mu
void update_mu(const double mu_prior_mean, const double mu_prior_var, RcppParallel::RMatrix<double> &mu, const int rung, const int K) {
  for (int k=0; k<K; k++) {
    mu(rung, k) = rnorm1(mu_prior_mean, mu_prior_var);
  }
}

