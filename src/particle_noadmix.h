
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class defining particle under no-admixture model
class particle_noadmix {
  
public:
  
  // PUBLIC OBJECTS
  
  // extract data and parameters
  std::vector<std::vector<int>> *data_ptr;
  std::vector<int> Jl;
  std::vector<int> ploidy;
  int n;
  int L;
  
  double beta_raised;
  int K;
  double lambda;
  
  // group allocation
  std::vector<int> group;
  
  // allele frequencies
  std::vector<std::vector<std::vector<double>>> log_allele_freqs;
  std::vector<std::vector<std::vector<int>>> allele_counts;
  
  // probabilities and likelihoods
  double loglike;
  
  // Q matrix
  std::vector<std::vector<double>> qmatrix_ind;
  std::vector<std::vector<double>> log_qmatrix_ind;
  
  // initialise ordering of labels
  std::vector<int> label_order;
  std::vector<int> label_order_new;
  
  // objects for solving label switching problem
  std::vector<std::vector<double>> cost_mat;
  std::vector<int> best_perm;
  std::vector<int> best_perm_order;
  std::vector<int> edges_left;
  std::vector<int> edges_right;
  std::vector<int> blocked_left;
  std::vector<int> blocked_right;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  particle_noadmix();
  particle_noadmix(std::vector<std::vector<int>> &data, std::vector<int> &Jl_, std::vector<int> &ploidy_, Rcpp::List &args_model, double beta_raised_);
  
  // other functions
  void reset();
  void update_group();
  void update_allele_freqs();
  void solve_label_switching(const std::vector<std::vector<double>> &log_qmatrix_running);
  void calculate_loglike();
  
};