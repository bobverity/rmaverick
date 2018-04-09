#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib rmaverick
#' @import assertthat
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs defaultNumThreads setThreadOptions
#' @import stats
NULL

#------------------------------------------------
#' Example MCMC
#'
#' Demonstrates the use of Rcpp and RcppParallel through a simple MCMC.
#'
#' @param x TODO
#' @param mu_prior_mean TODO
#' @param mu_prior_var TODO
#' @param burnin TODO
#' @param samples TODO
#' @param rungs TODO
#' @param num_threads TODO
#'
#' @export
#' @examples
#' example_mcmc(1:5)

example_mcmc <- function(x, mu_prior_mean=0, mu_prior_var=1e3, burnin=1e2, samples=1e3, rungs=10, num_threads=NULL) {
  
  # set defaults
  num_threads <- define_default(num_threads, defaultNumThreads())
  
  # check inputs
  assert_that(mu_prior_var > 0)
  assert_that(burnin > 0)
  assert_that(samples > 0)
  assert_that(rungs > 0)
  assert_that(num_threads <= defaultNumThreads(), msg = paste("num_threads is not less than or equal to", defaultNumThreads()))
  
  # set thread options
  setThreadOptions(numThreads = num_threads)
  
  # run efficient Rcpp function
  args_params <- list(mu_prior_mean = mu_prior_mean, mu_prior_var = mu_prior_var, burnin = burnin, samples = samples, rungs = rungs, num_threads = num_threads)
  output_raw <- example_mcmc_cpp(x, args_params)
  
  return(output_raw)
}
