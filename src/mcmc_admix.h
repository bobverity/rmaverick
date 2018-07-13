
#pragma once

#include <Rcpp.h>
#include "particle_admix.h"

//------------------------------------------------
// class defining MCMC under the admixture model
class mcmc_admix {
  
public:
  
  // PUBLIC OBJECTS
  
  // extract data and parameters
  std::vector<std::vector<int>> data;
  std::vector<int> Jl;
  std::vector<int> ploidy;
  int n;
  int L;
  
  int K;
  int burnin;
  int samples;
  int rungs;
  double GTI_pow;
  bool auto_converge;
  int converge_test;
  bool solve_label_switching_on;
  bool coupling_on;
  bool pb_markdown;
  bool silent;
  bool estimate_alpha;
  
  // thermodynamic parameters
  std::vector<double> beta_raised;
  std::vector<int> rung_order;
  int cold_rung;
  
  // vector of particles
  std::vector<particle_admix> particle_vec;
  
  // Q-matrices
  std::vector<std::vector<std::vector<double>>> log_qmatrix_gene_running;
  std::vector<std::vector<std::vector<double>>> qmatrix_gene;
  std::vector<std::vector<double>> qmatrix_ind;
  
  // ordering of labels
  std::vector<int> label_order;
  
  // objects for solving label switching problem
  std::vector<std::vector<double>> cost_mat;
  std::vector<int> best_perm;
  std::vector<int> best_perm_order;
  std::vector<int> edges_left;
  std::vector<int> edges_right;
  std::vector<int> blocked_left;
  std::vector<int> blocked_right;
  
  // objects for storing results
  std::vector<std::vector<double>> loglike_burnin;
  std::vector<std::vector<double>> loglike_sampling;
  std::vector<double> alpha_store;
  
  // objects for storing acceptance rates
  std::vector<int> coupling_accept;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  mcmc_admix(Rcpp::List &args_data, Rcpp::List &args_model);
  
  // other functions
  void burnin_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress);
  void sampling_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress);
  
  void update_log_qmatrix_gene_running();
  void update_qmatrix_gene();
  void metropolis_coupling();
  void solve_label_switching();
  
};
