
#------------------------------------------------
#' @title Import file
#'
#' @description Import file from the inst/extdata folder of this package
#'
#' @details TODO
#' 
#' @param name name of file
#'
#' @export
#' @examples
#' # TODO

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
# define default colours for barplots
# (not exported)
#' @noRd
default_colours <- function(K) {
  
  # generate palette and colours
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  my_palette <- colorRampPalette(raw_cols)
  ret <- my_palette(max(K,6))
  
  # if fewer than 6 colours then choose manually
  if (K==5) {
    ret = ret[c(1,2,3,5,6)]
  } else if (K==4) {
    ret = ret[c(1,2,3,5)]
  } else if (K==3) {
    ret = ret[c(1,3,5)]
  } else if (K==2) {
    ret = ret[c(1,5)]
  } else if (K==1) {
    ret = ret[1]
  }
  
  return(ret)
}

#------------------------------------------------
# return 95% quantile
# (not exported)
#' @noRd
quantile_95 <- function(x) {
  quantile(x, probs=c(0.025, 0.5, 0.975))
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
    ret <- FALSE;
  }
  return(ret)
}

#------------------------------------------------
# update progress bar
# (not exported)
#' @noRd
update_progress <- function(args, type, i, max_i) {
  
  # split by type
  if (type==1) { # scaffold progress bar
    setTxtProgressBar(args$pb_scaf, i)
    if (i==max_i) {
      close(args$pb_scaf)
    }
  } else if (type==2) { # burn-in iterations progress bar
    setTxtProgressBar(args$pb_burnin, i)
    if (i==max_i) {
      close(args$pb_burnin)
    }
  } else if (type==3) { # sampling iterations progress bar
    setTxtProgressBar(args$pb_samples, i)
    if (i==max_i) {
      close(args$pb_samples)
    }
  }
}
