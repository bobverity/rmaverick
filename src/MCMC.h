
#pragma once

#include <Rcpp.h>
#include "particle.h"

//------------------------------------------------
// class defining MCMC
class MCMC {
  
public:
  
  // PUBLIC OBJECTS
  
  // extract data and parameters
  std::vector<double>* x_ptr;
  int n;
  int K;
  int burnin;
  int samples;
  int iterations;
  int rungs;
  bool solve_label_switching_on;
  bool coupling_on;
  bool scaf_on;
  int scaf_n;
  bool splitmerge_on;
  bool parallel_on;
  bool print_console;
  
  // thermodynamic parameters
  std::vector<double> beta_vec;
  std::vector<int> rung_order;
  int cold_rung;
  
  // vector of particles
  std::vector<particle> particle_vec;
  
  // scaffold objects
  std::vector<std::vector<std::vector<int>>> scaf_group;
  std::vector<std::vector<std::vector<int>>> scaf_counts;
  std::vector<std::vector<std::vector<double>>> scaf_x_sum;
  
  // objects for storing results
  std::vector<std::vector<double>> mu_store;
  std::vector<std::vector<double>> loglike_store;
  std::vector<std::vector<double>> qmatrix_running;
  std::vector<std::vector<double>> log_qmatrix_running;
  std::vector<std::vector<double>> qmatrix_final;
  
  // objects for storing acceptance rates
  std::vector<int> mc_accept;
  std::vector<int> scaf_accept;
  std::vector<int> splitmerge_accept;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  MCMC(std::vector<double> &x, Rcpp::List &args);
  
  // other functions
  void scaffold_mcmc(Rcpp::List &args);
  void main_mcmc(Rcpp::List &args);
  void update_qmatrix_running();
  void update_qmatrix_final();
  void metropolis_coupling();
  
};
