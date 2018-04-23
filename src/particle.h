
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class defining particle
class particle {
  
public:
  
  // PUBLIC OBJECTS
  
  // extract data and parameters
  std::vector<double>* x_ptr;
  int n;
  int K;
  double mu_prior_mean;
  double mu_prior_var;
  double sigma;
  double sigma_sq;
  double beta;
  double scaf_n;
  
  // probabilities and likelihoods
  double loglike;
  std::vector<std::vector<double>> log_prob_mat;
  
  // component means
  std::vector<double> mu;
  
  // group allocation
  std::vector<int> group;
  std::vector<int> counts;
  std::vector<double> x_sum;
  
  // Q matrix
  std::vector<std::vector<double>> qmatrix;
  std::vector<std::vector<double>> log_qmatrix;
  
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
  
  // scaffold objects
  std::vector<double> scaf_mu;
  std::vector<std::vector<int>> scaf_group;
  std::vector<std::vector<int>> scaf_counts;
  std::vector<std::vector<double>> scaf_x_sum;
  
  // objects for split-merge
  std::vector<int> splitmerge_targets;
  std::vector<double> splitmerge_mu;
  std::vector<int> splitmerge_group;
  std::vector<int> splitmerge_counts;
  std::vector<double> splitmerge_x_sum;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  particle();
  particle(std::vector<double> &x, Rcpp::List &args_params, double beta_);
  
  // other functions
  void reset();
  void update_mu();
  void update_group();
  void group_increasing();
  double scaf_prop_logprob(const std::vector<int> &prop_group);
  void scaf_propose(int &scaf_accept);
  void splitmerge_propose(int &splitmerge_accept);
  void solve_label_switching(const std::vector<std::vector<double>> &log_qmatrix_running);
  void calc_log_like();
  
};