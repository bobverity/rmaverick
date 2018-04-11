
#include "main.h"
#include "misc.h"
#include "probability.h"
#include "MCMC.h"

using namespace std;

//------------------------------------------------
// Run MCMC
// [[Rcpp::export]]
Rcpp::List example_mcmc_cpp(Rcpp::NumericVector x, Rcpp::List args_params) {
  
  // create MCMC object
  MCMC mainMCMC(x, args_params);
  
  // run MCMC
  mainMCMC.parallel_mcmc();
  
  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( mainMCMC.mu_store ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("mu");
  
  ret.names() = ret_names;
  return ret;
}

