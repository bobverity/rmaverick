
#include "particle.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// default constructor for particle class
particle::particle() {}

//------------------------------------------------
// constructor for particle class
particle::particle(Rcpp::NumericVector x, Rcpp::List args_params) {
  
  // extract parameter values
  mu_prior_mean = rcpp_to_double(args_params["mu_prior_mean"]);
  mu_prior_var = rcpp_to_double(args_params["mu_prior_var"]);
  
}

//------------------------------------------------
// test function for particle class
void particle::dummy1() {
  bar();
}