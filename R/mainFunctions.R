#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib rmaverick
#' @import assertthat
#' @import parallel
#' @import coda
#' @importFrom Rcpp evalCpp
#' @import stats
#' @import utils
NULL

#------------------------------------------------
#' Example MCMC
#'
#' Demonstrates good algorithmic and coding practices for a simple normal mixture model. Implements advanced MCMC methods, including Metropolis-coupling over temperature rungs, a split-merge proposal, and an additional form of proposal referred to here as a "scaffold" proposal. Uses Rcpp and the package parallel for speed and efficiency.
#'
#' @param x the raw data (vector)
#' @param K the number of mixture components
#' @param mu_prior_mean the mean of the (normal) prior on mixture component locations
#' @param mu_prior_var the variance of the (normal) prior on mixture component locations
#' @param sigma the standard deviation of mixture components
#' @param burnin the number of burn-in iterations
#' @param samples the number of sampling iterations
#' @param rungs the number of temperature rungs
#' @param solve_label_switching_on whether to implement the Stevens' solution to the label-switching problem
#' @param coupling_on whether to implement Metropolis-coupling over temperature rungs
#' @param scaffold_on whether to use scaffolds to improve mixing
#' @param scaffold_n the number of scaffolds to use
#' @param splitmerge_on whether to implement a split-merge proposal
#' @param parallel_on whether to run each value of K in parallel
#' @param num_cores number of cores to use in parallelisation
#'
#' @export
#' @examples
#' # run example MCMC
#' m <- example_mcmc(1:5, parallel_on = FALSE)

example_mcmc <- function(x, K = 3, mu_prior_mean = 0, mu_prior_var = 1e3, sigma = 1, burnin =1e2, samples = 1e3, rungs = 10, solve_label_switching_on = TRUE, coupling_on = TRUE, scaffold_on = TRUE, scaffold_n = 10, splitmerge_on = TRUE, parallel_on = TRUE, num_cores = NULL) {
  
  # set defaults
  num_cores <- define_default(num_cores, detectCores())
  
  # check inputs
  assert_that(length(x) > 1)
  assert_that(all(K > 1))
  assert_that(mu_prior_var > 0)
  assert_that(sigma > 0)
  assert_that(burnin > 0)
  assert_that(samples > 0)
  assert_that(rungs > 0)
  assert_that(scaffold_n > 0)
  assert_that(num_cores > 0)
  assert_that(num_cores <= detectCores(), msg = paste("num_cores must be less than or equal to", detectCores()))
  
  # define argument list
  parallel_args <- list()
  for (i in 1:length(K)) {
    
    # create progress bars
    pb_scaf <- txtProgressBar(min = 0, max = scaffold_n, initial = NA, style = 3)
    pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
    pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
    
    # create argument list
    parallel_args[[i]] <- list(x = x, K = K[i], mu_prior_mean = mu_prior_mean, mu_prior_var = mu_prior_var, sigma = sigma, burnin = burnin, samples = samples, rungs = rungs, solve_label_switching_on = solve_label_switching_on, coupling_on = coupling_on, scaffold_on = scaffold_on, scaf_n = scaffold_n, splitmerge_on = splitmerge_on, parallel_on = parallel_on, test_convergence = test_convergence, update_progress = update_progress, pb_scaf = pb_scaf, pb_burnin = pb_burnin, pb_samples = pb_samples)
  }
  
  #------------------------
  
  # run efficient Rcpp function
  if (parallel_on) {  # run in parallel
    
    cl <- makeCluster(num_cores, type="FORK")
    output_raw <- clusterApplyLB(cl = cl, parallel_args, example_mcmc_cpp)
    stopCluster(cl)
    
  } else {  # run in serial
    
    output_raw <- list()
    for (i in 1:length(K)) {
      output_raw[[i]] <- example_mcmc_cpp(parallel_args[[i]])
    }
    
  }
  
  #------------------------
  
  # begin processing results
  if (!parallel_on) {
    cat("Processing results\n")
  }
  
  # loop through K
  ret <- list()
  for (i in 1:length(K)) {
    
    # get loglikelihood
    loglike <- t(rcpp_to_mat(output_raw[[i]]$loglike))
    
    # split into burn-in vs sampling iterations
    loglike_burnin <- mcmc(loglike[1:burnin,])
    loglike_sampling <- mcmc(loglike[-(1:burnin),])
    
    # get quantiles over sampling loglikelihoods
    loglike_quantiles <- t(apply(loglike_sampling, 2, function(x){quantile(x, probs=c(0.025, 0.5, 0.975))}))
    
    # process mu draws
    mu <- rcpp_to_mat(output_raw[[i]]$mu)
    
    # process Q-matrix
    qmatrix <- rcpp_to_mat(output_raw[[i]]$qmatrix)
    
    # process acceptance rates
    mc_accept <- 100 * output_raw[[i]]$mc_accept/(burnin+samples)
    scaf_accept <- 100 * output_raw[[i]]$scaf_accept/(burnin+samples)
    splitmerge_accept <- 100 * output_raw[[i]]$splitmerge_accept/(burnin+samples)
    
    # store results of this K
    ret[[i]] <- list(loglike_burnin = loglike_burnin, loglike_sampling = loglike_sampling, loglike_quantiles = loglike_quantiles, mu = mu, qmatrix = qmatrix, mc_accept = mc_accept, scaf_accept = scaf_accept, splitmerge_accept = splitmerge_accept)
  }
  
  # return as list
  return(ret)
}

#------------------------------------------------
# geweke_pvalue
# return p-value of Geweke's diagnostic convergence statistic, estimated from package coda
# (not exported)

geweke_pvalue <- function(x) {
  ret <- 2*pnorm(abs(geweke.diag(x)$z), lower.tail=FALSE)
  return(ret)
}

#------------------------------------------------
# test convergence
# (not exported)

test_convergence <- function(x) {
  
  # get Geweke p-value
  g <- geweke_pvalue(mcmc(x))
  ret <- (g>0.05)
  
  return(ret)
}

#------------------------------------------------
# update progress bar
# (not exported)

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

#------------------------------------------------
# call_hungarian
# calls C++ implementation of the Hungarian algorithm for binding best matching
# in a linear sum assigment problem. This is function is used in testing.
# (not exported)

call_hungarian <- function(x) {
  args <- list(cost_mat = mat_to_rcpp(x))
  call_hungarian_cpp(args)
}

