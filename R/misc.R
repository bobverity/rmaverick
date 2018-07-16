
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
# (not exported)
#' @noRd
define_default <- function(x, default_value) {
  if (is.null(x)) {
    x <- default_value
  }
  return(x)
}

#------------------------------------------------
# force scalar to vector by repeating value 
# (not exported)
#' @noRd
force_vector <- function(x, l) {
  if (length(x)==1) {
    x <- rep(x, l)
  }
  return(x)
}

# -----------------------------------
# ask user a yes/no question. Return TRUE/FALSE.
# (not exported)
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
# (not exported)
#' @noRd
rdirichlet <- function(alpha_vec) {
  z <- rgamma(length(alpha_vec), shape = alpha_vec, scale = 1)
  ret <- z/sum(z)
  return(ret)
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
# (not exported)
#' @noRd
mat_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix.
# (not exported)
#' @noRd
rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  return(ret)
}

#------------------------------------------------
# return 95% quantile
# (not exported)
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs=c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#------------------------------------------------
# sum logged values without underflow, i.e. do log(sum(exp(x)))
# (not exported)
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
# (not exported)
#' @noRd
geweke_pvalue <- function(x) {
  ret <- 2*pnorm(abs(geweke.diag(x)$z), lower.tail=FALSE)
  return(ret)
}

#------------------------------------------------
# check that geweke p-value non-significant on values x[1:n]
# (not exported)
#' @noRd
test_convergence <- function(x, n) {
  if (n==1) {
    return(FALSE)
  }
  g <- geweke_pvalue(mcmc(x[1:n]))
  ret <- (g>0.01)
  if (is.na(ret)) {
    ret <- TRUE;
  }
  return(ret)
}

#------------------------------------------------
# update progress bar
# (not exported)
#' @noRd
update_progress <- function(pb_list, name, i, max_i) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i==max_i) {
    close(pb_list[[name]])
  }
}

#------------------------------------------------
# function for determining if object is of class cluster
# (not exported)
#' @noRd
is.cluster <- function(x) {
  inherits(x, "cluster")
}

################################################################
# MISC CLASSES

#------------------------------------------------
#' @title Determine if object is of class maverick_GTI_path
#'
#' @description Determine if object is of class maverick_GTI_path
#' 
#' @param x any R object
#'
#' @export

is.maverick_GTI_path <- function(x) {
  inherits(x, "maverick_GTI_path")
}

#------------------------------------------------
#' @title Overload print() function for class maverick_GTI_path
#'
#' @description Overload print() function for class maverick_GTI_path
#' 
#' @param x an object of class \code{maverick_GTI_path}
#' @param ... other arguments passed to \code{print()} function
#'
#' @export

print.maverick_GTI_path <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
#' @title Determine if object is of class maverick_qmatrix_ind
#'
#' @description Determine if object is of class maverick_qmatrix_ind
#' 
#' @param x any R object
#'
#' @export

is.maverick_qmatrix_ind <- function(x) {
  inherits(x, "maverick_qmatrix_ind")
}

#------------------------------------------------
#' @title Overload print() function for class maverick_qmatrix_ind
#'
#' @description Overload print() function for class maverick_qmatrix_ind
#' 
#' @param x an object of class \code{maverick_qmatrix_ind}
#' @param ... other arguments passed to \code{print()} function
#'
#' @export

print.maverick_qmatrix_ind <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
#' @title Determine if object is of class maverick_loglike_quantiles
#'
#' @description Determine if object is of class maverick_loglike_quantiles
#' 
#' @param x any R object
#'
#' @export

is.maverick_loglike_quantiles <- function(x) {
  inherits(x, "maverick_loglike_quantiles")
}

#------------------------------------------------
#' @title Overload print() function for class maverick_loglike_quantiles
#'
#' @description Overload print() function for class maverick_loglike_quantiles
#' 
#' @param x an object of class \code{maverick_loglike_quantiles}
#' @param ... other arguments passed to \code{print()} function
#'
#' @export

print.maverick_loglike_quantiles <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
#' @title Determine if object is of class maverick_GTI_logevidence
#'
#' @description Determine if object is of class maverick_GTI_logevidence
#' 
#' @param x any R object
#'
#' @export

is.maverick_GTI_logevidence <- function(x) {
  inherits(x, "maverick_GTI_logevidence")
}

#------------------------------------------------
#' @title Overload print() function for class maverick_GTI_logevidence
#'
#' @description Overload print() function for class maverick_GTI_logevidence
#' 
#' @param x an object of class \code{maverick_GTI_logevidence}
#' @param ... other arguments passed to \code{print()} function
#'
#' @export

print.maverick_GTI_logevidence <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
#' @title Determine if object is of class maverick_GTI_posterior
#'
#' @description Determine if object is of class maverick_GTI_posterior
#' 
#' @param x any R object
#'
#' @export

is.maverick_GTI_posterior <- function(x) {
  inherits(x, "maverick_GTI_logevidence")
}

#------------------------------------------------
#' @title Overload print() function for class maverick_GTI_posterior
#'
#' @description Overload print() function for class maverick_GTI_posterior
#' 
#' @param x an object of class \code{maverick_GTI_posterior}
#' @param ... other arguments passed to \code{print()} function
#'
#' @export

print.maverick_GTI_posterior <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
#' @title Determine if object is of class maverick_GTI_logevidence_model
#'
#' @description Determine if object is of class maverick_GTI_logevidence_model
#' 
#' @param x any R object
#'
#' @export

is.maverick_GTI_logevidence_model <- function(x) {
  inherits(x, "maverick_GTI_logevidence")
}

#------------------------------------------------
#' @title Overload print() function for class maverick_GTI_logevidence_model
#'
#' @description Overload print() function for class maverick_GTI_logevidence_model
#' 
#' @param x an object of class \code{maverick_GTI_logevidence_model}
#' @param ... other arguments passed to \code{print()} function
#'
#' @export

print.maverick_GTI_logevidence_model <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
#' @title Determine if object is of class maverick_GTI_posterior_model
#'
#' @description Determine if object is of class maverick_GTI_posterior_model
#' 
#' @param x any R object
#'
#' @export

is.maverick_GTI_posterior_model <- function(x) {
  inherits(x, "maverick_GTI_logevidence")
}

#------------------------------------------------
#' @title Overload print() function for class maverick_GTI_posterior_model
#'
#' @description Overload print() function for class maverick_GTI_posterior_model
#' 
#' @param x an object of class \code{maverick_GTI_posterior_model}
#' @param ... other arguments passed to \code{print()} function
#'
#' @export

print.maverick_GTI_posterior_model <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

