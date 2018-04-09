
#pragma once

#include <Rcpp.h>
#include <RcppParallel.h>

//------------------------------------------------
// Run MCMC
Rcpp::List example_mcmc_cpp(Rcpp::NumericVector x, Rcpp::List args_params);

//------------------------------------------------
// Run serial version of MCMC
void serial_mcmc(Rcpp::NumericVector &x);

//------------------------------------------------
// Run parallel version of MCMC
void parallel_mcmc(Rcpp::NumericVector &x, Rcpp::NumericMatrix &mu_store, int num_threads);