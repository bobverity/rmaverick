
#pragma once

#include <Rcpp.h>
#include "misc.h"

//------------------------------------------------
// run main MCMC
Rcpp::List run_mcmc_cpp(Rcpp::List args);

//------------------------------------------------
// get qmatrix_ind
template<class TYPE>
std::vector<std::vector<double>> get_qmatrix_ind(TYPE x) {
  return x.qmatrix_ind;
}

//------------------------------------------------
// get loglike_burnin
template<class TYPE>
std::vector<std::vector<double>> get_loglike_burnin(TYPE x) {
  return x.loglike_burnin;
}

//------------------------------------------------
// get loglike_sampling
template<class TYPE>
std::vector<std::vector<double>> get_loglike_sampling(TYPE x) {
  return x.loglike_sampling;
}

//------------------------------------------------
// get coupling_accept
template<class TYPE>
std::vector<int> get_coupling_accept(TYPE x) {
  return x.coupling_accept;
}

//------------------------------------------------
// get alpha_store
template<class TYPE>
std::vector<double> get_alpha_store(TYPE x) {
  return x.alpha_store;
}

//------------------------------------------------
// estimate evidence quantiles by simulation
Rcpp::List GTI_evidence_sim_cpp(Rcpp::List args);

// integrate log-evidence over K by simulation
Rcpp::List GTI_integrated_K_sim_cpp(Rcpp::List args);
