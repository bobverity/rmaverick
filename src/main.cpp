
#include "main.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// Run MCMC
// [[Rcpp::export]]
Rcpp::List example_mcmc_cpp(Rcpp::NumericVector x, Rcpp::List args_params) {
  
  // extract parameter values
  double mu_prior_mean = rcpp_to_double(args_params["mu_prior_mean"]);
  double mu_prior_var = rcpp_to_double(args_params["mu_prior_var"]);
  int burnin = rcpp_to_int(args_params["burnin"]);
  int samples = rcpp_to_int(args_params["samples"]);
  int rungs = rcpp_to_int(args_params["rungs"]);
  int num_threads = rcpp_to_int(args_params["num_threads"]);
  
  int iterations = burnin + samples;
  
  // define objects for storing results
  Rcpp::NumericMatrix mu_store(rungs, iterations);
  
  //rcpp_print_vector(mu_store( Rcpp::_, 0 ));
  //print_stars();
  
  //---------------------
  
  // carry out MCMC in serial or in parallel
  //linear_mcmc(x);
  parallel_mcmc(x, mu_store, num_threads);
  
  //---------------------
  
  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( mu_store ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("mu");
  
  ret.names() = ret_names;
  return ret;
}

//------------------------------------------------
// Run serial version of MCMC
void serial_mcmc(Rcpp::NumericVector &x) {
  
  print("serial MCMC");
  
}

//------------------------------------------------
// Functor to run multiple MCMC chains
struct functor_run_chains : public RcppParallel::Worker {
  
  // source vector
  const RcppParallel::RVector<double> input;
  
  // destination vector
  RcppParallel::RMatrix<double> output;
  
  // initialize with input and output
  functor_run_chains(const Rcpp::NumericVector input, Rcpp::NumericMatrix output)
    : input(input), output(output) {}
  
  // operator to be called by parallel functions
  void operator()(size_t begin, size_t end) {
    
    // loop over MCMC chains
    for (size_t i=begin; i<end; i++) {
      
      // MCMC
      for (int j=0; j<output.ncol(); j++) {
        output(i,j) = rnorm1(0, 1);
      }
      
    }
    
  }
};

//------------------------------------------------
// Run parallel version of MCMC
void parallel_mcmc(Rcpp::NumericVector &x, Rcpp::NumericMatrix &mu_store, int num_threads) {
  
  // print header
  print("Running", num_threads, "threads in parallel");
  
  // pass input and output to functor
  functor_run_chains run_chains(x, mu_store);
  
  // call parallelFor
  parallelFor(0, mu_store.nrow(), run_chains);
  
}
