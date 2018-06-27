
#include "misc.h"
#include "Hungarian.h"
#include "test_wrappers.h"

using namespace std;

//------------------------------------------------
// call Hungarian algorithm for binding best matching in a linear sum assigment problem
// [[Rcpp::export]]
Rcpp::List call_hungarian_cpp(Rcpp::List args) {
  
  // objects for calling Hungarian algorithm
  vector<vector<double>> cost_mat = rcpp_to_mat_double(args["cost_mat"]);
  int n = cost_mat.size();
  vector<int> edges_left(n);
  vector<int> edges_right(n);
  vector<int> blocked_left(n);
  vector<int> blocked_right(n);
  vector<int> best_perm = hungarian(cost_mat, edges_left, edges_right, blocked_left, blocked_right);
  
  // return
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("best_matching") = best_perm);
  return ret;
}
