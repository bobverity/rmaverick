
#pragma once

#include <Rcpp.h>
#include <RcppParallel.h>

//------------------------------------------------
// Run MCMC
Rcpp::List example_mcmc_cpp(Rcpp::NumericVector x, Rcpp::List args_params);
