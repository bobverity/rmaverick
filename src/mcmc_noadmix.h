
#pragma once

#include <Rcpp.h>
#include "particle_noadmix.h"

//------------------------------------------------
// class defining MCMC under the no-admixture model
class mcmc_noadmix {
  
public:
  
  // PUBLIC OBJECTS
  
  // extract data and parameters
  std::vector<std::vector<int>> data;
  int n;
  std::vector<int> Jl;
  std::vector<int> ploidy;
  
  int K;
  int burnin;
  int samples;
  int rungs;
  double GTI_pow;
  bool auto_converge;
  bool solve_label_switching_on;
  bool coupling_on;
  bool splitmerge_on;
  bool silent;
  
  // thermodynamic parameters
  std::vector<double> beta_raised;
  std::vector<int> rung_order;
  int cold_rung;
  
  // vector of particles
  std::vector<particle_noadmix> particle_vec;
  
  // Q-matrices
  std::vector<std::vector<double>> log_qmatrix_ind_running;
  std::vector<std::vector<double>> qmatrix_ind;
  
  // objects for storing results
  std::vector<std::vector<double>> loglike_burnin;
  std::vector<std::vector<double>> loglike_sampling;
  
  // objects for storing acceptance rates
  std::vector<int> coupling_accept;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  mcmc_noadmix(Rcpp::List &args_data, Rcpp::List &args_model);
  
  // other functions
  void burnin_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress);
  void sampling_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress);
  
  void update_log_qmatrix_ind_running();
  void update_qmatrix_ind();
  void metropolis_coupling();
  
};