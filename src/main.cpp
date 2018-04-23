
#include "main.h"
#include "misc.h"
#include "probability.h"
#include "MCMC.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// Run MCMC
// [[Rcpp::export]]
Rcpp::List example_mcmc_cpp(Rcpp::List args) {
  
  // convert data to base C++ format
  vector<double> x = rcpp_to_vector_double(args["x"]);
  bool scaf_on = rcpp_to_bool(args["scaffold_on"]);
  
  // create MCMC object
  MCMC m(x, args);
  
  // generate scaffold groupings
  if (scaf_on) {
    m.scaffold_mcmc(args);
  }
  
  // run MCMC
  m.main_mcmc(args);
  
  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( m.loglike_store ));
  ret.push_back(Rcpp::wrap( m.mu_store ));
  ret.push_back(Rcpp::wrap( m.qmatrix_final ));
  ret.push_back(Rcpp::wrap( m.mc_accept ));
  ret.push_back(Rcpp::wrap( m.scaf_accept ));
  ret.push_back(Rcpp::wrap( m.splitmerge_accept ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("loglike");
  ret_names.push_back("mu");
  ret_names.push_back("qmatrix");
  ret_names.push_back("mc_accept");
  ret_names.push_back("scaf_accept");
  ret_names.push_back("splitmerge_accept");
  
  ret.names() = ret_names;
  return ret;
}

