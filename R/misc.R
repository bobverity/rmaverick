
#------------------------------------------------
# replace NULL value with default
# (not exported)

define_default <- function(x, default_value) {
  if (is.null(x)) {
    x <- default_value
  }
  return(x)
}

# -----------------------------------
# mat_to_Rcpp
# takes matrix as input, converts to list format for use within Rcpp code
# (not exported)

mat_to_rcpp <- function(x) {
  return(split(x, f=1:nrow(x)))
}

# -----------------------------------
# Rcpp_to_mat
# Takes list format returned from Rcpp and converts to matrix.
# (not exported)

rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  return(ret)
}
