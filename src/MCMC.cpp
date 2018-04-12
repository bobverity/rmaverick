
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
  sigma = rcpp_to_double(args_params["sigma"]);
  sigma_sq = sigma * sigma;
  burnin = rcpp_to_int(args_params["burnin"]);
  samples = rcpp_to_int(args_params["samples"]);
  iterations = burnin + samples;
  rungs = rcpp_to_int(args_params["rungs"]);
  mc_interval = rcpp_to_int(args_params["mc_interval"]);
  num_threads = rcpp_to_int(args_params["num_threads"]);
  
  // thermodynamic powers
  beta = Rcpp::NumericVector(rungs);
  for (int i=0; i<rungs; i++) {
    beta(i) = (rungs==1) ? 1 : i/double(rungs-1);
  }
  
  // initialize component means
  mu = Rcpp::NumericMatrix(rungs, K);
  
  // initialize group allocation and counts
  group = Rcpp::IntegerMatrix(rungs, n);
  x_sum = Rcpp::NumericMatrix(rungs, K);
  counts = Rcpp::IntegerMatrix(rungs, K);
  for (int i=0; i<rungs; i++) {
    for (int j=0; j<n; j++) {
      int this_group = sample2(0, K-1);
      group(i, j) = this_group;
      x_sum(i, this_group) += x[j];
      counts(i, this_group) ++;
    }
  }
  
  // Q-matrix
  q_matrix = Rcpp::NumericMatrix(rungs, K);
  
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
  RcppParallel::RVector<double> x;
  RcppParallel::RMatrix<int> group;
  RcppParallel::RMatrix<double> x_sum;
  RcppParallel::RMatrix<int> counts;
  RcppParallel::RMatrix<double> q_matrix;
  RcppParallel::RMatrix<double> mu;
  const double mu_prior_mean;
  const double mu_prior_var;
  const double sigma_sq;
  const int K;
  const RcppParallel::RVector<double> beta;
  
  // output objects
  RcppParallel::RMatrix<double> mu_store;
  
  // initialize with input and output
  functor_run_chains(Rcpp::IntegerVector start_it, Rcpp::IntegerVector end_it, const Rcpp::NumericVector x, Rcpp::IntegerMatrix group, Rcpp::NumericMatrix x_sum, Rcpp::IntegerMatrix counts, Rcpp::NumericMatrix q_matrix, Rcpp::NumericMatrix mu, const double mu_prior_mean, const double mu_prior_var, const double sigma_sq, const int K, const Rcpp::NumericVector beta, Rcpp::NumericMatrix mu_store)
    : start_it(start_it), end_it(end_it), x(x), group(group), x_sum(x_sum), counts(counts), q_matrix(q_matrix), mu(mu), mu_prior_mean(mu_prior_mean), mu_prior_var(mu_prior_var), sigma_sq(sigma_sq), K(K), beta(beta), mu_store(mu_store) {}
  
  // operator to be called by parallel functions
  void operator()(size_t begin, size_t end) {
    
    // loop over MCMC chains
    for (size_t rung=begin; rung<end; rung++) {
      
      // MCMC for this chunk
      for (int rep=start_it[rung]; rep<end_it[rung]; rep++) {
        
        // update component means and grouping
        update_mu_group(x, group, x_sum, counts, q_matrix, mu, mu_prior_mean, mu_prior_var, sigma_sq, beta, rung);
        
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
  functor_run_chains run_chains(start_it, end_it, x, group, x_sum, counts, q_matrix, mu, mu_prior_mean, mu_prior_var, sigma_sq, K, beta, mu_store);
  
  // carry out MCMC in chunks
  for (int i=0; i<(iterations/mc_interval); i++) {
    
    //print("chunk", i);
    
    // run this chunk in parallel
    parallelFor(0, rungs, run_chains);
    
    // move to next chunk
    fill(start_it.begin(), start_it.end(), end_it[0]);
    fill(end_it.begin(), end_it.end(), end_it[0] + mc_interval);
    
    //rcpp_print_vector(start_it);
  }
  
}

//------------------------------------------------
// update mu and group
void update_mu_group(RcppParallel::RVector<double> &x, RcppParallel::RMatrix<int> &group, RcppParallel::RMatrix<double> &x_sum, RcppParallel::RMatrix<int> &counts, RcppParallel::RMatrix<double> &q_matrix, RcppParallel::RMatrix<double> &mu, const double mu_prior_mean, const double mu_prior_var, const double sigma, const RcppParallel::RVector<double> &beta, const int rung) {
  
  int n = group.ncol();
  int K = counts.ncol();
  
  /// update mu
  for (int k=0; k<K; k++) {
    double mu_post_var = 1/(counts(rung,k)/(sigma*sigma) + 1/mu_prior_var);
    double mu_post_mean = mu_post_var * x_sum(rung,k)/(sigma*sigma) + mu_prior_mean/mu_prior_var;
    mu(rung, k) = rnorm1(mu_post_mean, mu_post_var);
  }
  
  // update group
  for (int i=0; i<n; i++) {
    double prob_vec_sum = 0;
    for (int k=0; k<K; k++) {
      q_matrix(rung, k) = dnorm2(x, mu, sigma, i, rung, k, true);
      prob_vec_sum += q_matrix(rung, k);
    }
    //print(prob_vec_sum);
  }
  
}

