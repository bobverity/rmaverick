
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
  double sigma;
  double sigma_sq;
  int burnin;
  int samples;
  int iterations;
  int rungs;
  int mc_interval;
  int num_threads;
  
  // thermodynamic powers
  Rcpp::NumericVector beta;
  
  // component means and group allocation
  Rcpp::NumericMatrix mu;
  Rcpp::IntegerMatrix group;
  Rcpp::NumericMatrix x_sum;
  Rcpp::IntegerMatrix counts;
  
  // Q-matrix
  Rcpp::NumericMatrix q_matrix;
  
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
void update_mu_group(RcppParallel::RVector<double> &x, RcppParallel::RMatrix<int> &group, RcppParallel::RMatrix<double> &x_sum, RcppParallel::RMatrix<int> &counts, RcppParallel::RMatrix<double> &q_matrix, RcppParallel::RMatrix<double> &mu, const double mu_prior_mean, const double mu_prior_var, const double sigma, const RcppParallel::RVector<double> &beta, const int rung);

