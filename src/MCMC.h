
#pragma once

#include <Rcpp.h>
#include <RcppParallel.h>
#include "particle.h"
//#include "probability.h"

//------------------------------------------------
// class defining MCMC
class MCMC {
  
public:
  
  // PUBLIC OBJECTS
  
  // data and basic data properties
  Rcpp::NumericVector x;
  int n;
  int K;
  double mu_prior_mean;
  double mu_prior_var;
  int burnin;
  int samples;
  int iterations;
  int rungs;
  int mc_interval;
  int num_threads;
  
  // thermodynamic powers
  Rcpp::NumericVector beta_vec;
  
  // component means and group allocation
  Rcpp::NumericMatrix mu;
  Rcpp::IntegerMatrix group;
  
  // define objects for storing results
  Rcpp::NumericMatrix mu_store;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  MCMC(Rcpp::NumericVector x, Rcpp::List args_params);
  
  // other functions
  void serial_mcmc();
  void parallel_mcmc();
  
};

//------------------------------------------------
// update mu
void update_mu(const double mu_prior_mean, const double mu_prior_var, RcppParallel::RMatrix<double> &mu, const int rung, const int k);

