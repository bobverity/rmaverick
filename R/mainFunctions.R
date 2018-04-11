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
#' @param mc_interval TODO
#' @param num_threads TODO
#'
#' @export
#' @examples
#' example_mcmc(1:5)

example_mcmc <- function(x, K = 3, mu_prior_mean = 0, mu_prior_var = 1e3, burnin =1e2, samples = 1e3, rungs = 10, mc_interval = 1e2, num_threads = NULL, parallel_on = TRUE) {
  
  # set defaults
  num_threads <- define_default(num_threads, defaultNumThreads())
  
  # check inputs
  assert_that(K > 0)
  assert_that(mu_prior_var > 0)
  assert_that(burnin > 0)
  assert_that(samples > 0)
  assert_that(rungs > 0)
  assert_that(mc_interval > 0)
  assert_that(num_threads <= defaultNumThreads(), msg = paste("num_threads is not less than or equal to", defaultNumThreads()))
  
  # define composite parameters
  iterations <- burnin + samples
  
  # check that mc_interval divides simply into both burnin and samples
  if (parallel_on && (burnin %% mc_interval)!=0) {
    burnin <- ceiling(burnin/mc_interval) * mc_interval
    cat(paste("Note: increasing burnin to", burnin, "to make simple product of mc_interval\n"))
  }
  if (parallel_on && (samples %% mc_interval)!=0) {
    samples <- ceiling(samples/mc_interval) * mc_interval
    cat(paste("Note: increasing samples to", samples, "to make simple product of mc_interval\n"))
  }
  
  # set thread options
  setThreadOptions(numThreads = num_threads)
  
  # run efficient Rcpp function
  args_params <- list(K = K, mu_prior_mean = mu_prior_mean, mu_prior_var = mu_prior_var, burnin = burnin, samples = samples, rungs = rungs, mc_interval = mc_interval, num_threads = num_threads, parallel_on = parallel_on)
  output_raw <- example_mcmc_cpp(x, args_params)
  
  return(output_raw)
  
  # process raw output
  mu <- list()
  for (rung in 1:rungs) {
    mu[[rung]] <- matrix(unlist(output_raw$mu[[rung]]), nrow=iterations)
  }
  
  # return as list
  ret <- list(mu=mu)
  return(ret)
}
