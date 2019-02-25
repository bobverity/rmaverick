
#------------------------------------------------
# assert x is class mavproject
assert_mavproject <- function(x, message = "%s must be of class 'mavproject'", name = deparse(substitute(x))) {
  if (!is.mavproject(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# assert x is class cluster
assert_cluster <- function(x, message = "%s must be of class 'cluster'", name = deparse(substitute(x))) {
  if (!is.cluster(x)) {
    stop(sprintf(message, name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
#' @title Import file
#'
#' @description Import file from the inst/extdata folder of this package
#' 
#' @param name name of file
#'
#' @export

maverick_file <- function(name) {
  
  # load file from inst/extdata folder
  name_full <- system.file("extdata/", name, package = 'maverick', mustWork = TRUE)
  ret <- readRDS(name_full)
  
  # return
  return(ret)
}

#------------------------------------------------
# replace NULL value with default
#' @noRd
define_default <- function(x, default_value) {
  if (is.null(x)) {
    x <- default_value
  }
  return(x)
}

#------------------------------------------------
# force scalar to vector by repeating value 
#' @noRd
force_vector <- function(x, l) {
  if (length(x)==1) {
    x <- rep(x, l)
  }
  return(x)
}

# -----------------------------------
# ask user a yes/no question. Return TRUE/FALSE.
#' @noRd
user_yes_no <- function(x = "continue? (Y/N): ") {
  
  user_choice <- NA
  while (!user_choice %in% c("Y", "y" ,"N", "n")) {
    user_choice <- readline(x)
  }
  return(user_choice %in% c("Y", "y"))
}

# -----------------------------------
# draw from dirichlet distribution with vector of parameter inputs
#' @noRd
rdirichlet <- function(alpha_vec) {
  z <- rgamma(length(alpha_vec), shape = alpha_vec, scale = 1)
  ret <- z/sum(z)
  return(ret)
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
mat_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix.
#' @noRd
rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  return(ret)
}

#------------------------------------------------
# return 95% quantile
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs=c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#------------------------------------------------
# sum logged values without underflow, i.e. do log(sum(exp(x)))
#' @noRd
log_sum <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  x_max <- max(x, na.rm = TRUE)
  ret <- x_max + log(sum(exp(x-x_max)))
  return(ret)
}

#------------------------------------------------
# geweke_pvalue
# return p-value of Geweke's diagnostic convergence statistic, estimated from package coda
#' @noRd
geweke_pvalue <- function(x) {
  ret <- 2*pnorm(abs(geweke.diag(x)$z), lower.tail=FALSE)
  return(ret)
}

#------------------------------------------------
# check convergence on values x[1:n]
#' @noRd
test_convergence <- function(x, n) {
  # fail if n = 1
  if (n == 1) {
    return(FALSE)
  }
  
  # fail if ESS too small
  ESS <- coda::effectiveSize(x[1:n])
  if (ESS < 10) {
    return(FALSE)
  }
  
  # fail if geweke p-value < threshold
  g <- geweke_pvalue(mcmc(x[1:n]))
  ret <- (g > 0.01)
  if (is.na(ret)) {
    ret <- FALSE;
  }
  return(ret)
}

#------------------------------------------------
# update progress bar
#' @noRd
update_progress <- function(pb_list, name, i, max_i) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i==max_i) {
    close(pb_list[[name]])
  }
}

#------------------------------------------------
# function for determining if object is of class cluster
#' @noRd
is.cluster <- function(x) {
  inherits(x, "cluster")
}

################################################################
# MISC CLASSES

#------------------------------------------------
# Determine if object is of class maverick_GTI_path
#' @noRd
is.maverick_GTI_path <- function(x) {
  inherits(x, "maverick_GTI_path")
}

#------------------------------------------------
# Overload print() function for class maverick_GTI_path
#' @noRd
print.maverick_GTI_path <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
# Determine if object is of class maverick_qmatrix_ind
#' @noRd
is.maverick_qmatrix_ind <- function(x) {
  inherits(x, "maverick_qmatrix_ind")
}

#------------------------------------------------
# Overload print() function for class maverick_qmatrix_ind
#' @noRd
print.maverick_qmatrix_ind <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
# Determine if object is of class maverick_loglike_quantiles
#' @noRd
is.maverick_loglike_quantiles <- function(x) {
  inherits(x, "maverick_loglike_quantiles")
}

#------------------------------------------------
# Overload print() function for class maverick_loglike_quantiles
#' @noRd
print.maverick_loglike_quantiles <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
# Determine if object is of class maverick_GTI_logevidence
#' @noRd
is.maverick_GTI_logevidence <- function(x) {
  inherits(x, "maverick_GTI_logevidence")
}

#------------------------------------------------
# Overload print() function for class maverick_GTI_logevidence
#' @noRd
print.maverick_GTI_logevidence <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
# Determine if object is of class maverick_GTI_posterior
#' @noRd
is.maverick_GTI_posterior <- function(x) {
  inherits(x, "maverick_GTI_logevidence")
}

#------------------------------------------------
# Overload print() function for class maverick_GTI_posterior
#' @noRd
print.maverick_GTI_posterior <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
# Determine if object is of class maverick_GTI_logevidence_model
#' @noRd
is.maverick_GTI_logevidence_model <- function(x) {
  inherits(x, "maverick_GTI_logevidence")
}

#------------------------------------------------
# Overload print() function for class maverick_GTI_logevidence_model
#' @noRd
print.maverick_GTI_logevidence_model <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
# Determine if object is of class maverick_GTI_posterior_model
#' @noRd
is.maverick_GTI_posterior_model <- function(x) {
  inherits(x, "maverick_GTI_logevidence")
}

#------------------------------------------------
# Overload print() function for class maverick_GTI_posterior_model
#' @noRd
print.maverick_GTI_posterior_model <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

