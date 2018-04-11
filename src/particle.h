
#pragma once

#include <Rcpp.h>
#include <RcppParallel.h>

//------------------------------------------------
// class defining particle
class particle {
  
public:
  
  // PUBLIC OBJECTS
  
  // data and basic data properties
  // extract parameter values
  double mu_prior_mean;
  double mu_prior_var;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  particle();
  particle(Rcpp::NumericVector x, Rcpp::List args_params);
  
  // other functions
  void dummy1();
  
};