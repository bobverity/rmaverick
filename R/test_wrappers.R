
#------------------------------------------------
# calls C++ implementation of the Hungarian algorithm for binding best matching
# in a linear sum assigment problem. This is function is used in testing.
# (not exported)
#' @noRd
call_hungarian <- function(x) {
  args <- list(cost_mat = mat_to_rcpp(x))
  call_hungarian_cpp(args)
}
