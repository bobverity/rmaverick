
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class defining particle under admixture model
class particle_admix {
  
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
  double alpha;
  
  // group allocation
  std::vector<std::vector<int>> group;
  
  // allele frequencies
  std::vector<std::vector<std::vector<double>>> log_allele_freqs;
  std::vector<std::vector<std::vector<int>>> allele_counts;
  
  // admixture proportions
  std::vector<std::vector<double>> log_admix_freqs;
  std::vector<std::vector<int>> admix_counts;
  std::vector<int> admix_counts_totals;
  
  // individual-level group update objects
  std::vector<double> log_admix_freqs_prop;
  std::vector<std::vector<std::vector<double>>> log_qmatrix_prop;
  std::vector<std::vector<std::vector<double>>> qmatrix_prop;
  std::vector<std::vector<int>> group_prop;
  
  // probabilities and likelihoods
  double loglike;
  
  // Q matrices
  std::vector<std::vector<std::vector<double>>> log_qmatrix_gene;
  std::vector<std::vector<std::vector<double>>> qmatrix_gene;
  
  // ordering of labels
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
  particle_admix();
  particle_admix(std::vector<std::vector<int>> &data, std::vector<int> &Jl_, std::vector<int> &ploidy_, Rcpp::List &args_model, double beta_raised_);
  
  // other functions
  void reset();
  void update_group();
  void update_allele_admix_freqs();
  void update_alpha();
  void solve_label_switching(const std::vector<std::vector<std::vector<double>>> &log_qmatrix_gene_running);
  void calculate_loglike();
  
};