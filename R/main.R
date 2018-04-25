
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib rmaverick
#' @import assertthat
#' @import parallel
#' @import coda
#' @importFrom Rcpp evalCpp
#' @import stats
#' @import utils
#' @importFrom grDevices colorRampPalette
NULL

#------------------------------------------------
# load data

load_data <- function(project, df, ID_col = 1, pop_col = NULL, ploidy_col = NULL, data_cols = NULL, ploidy = 1, missing_data = -9, name = NULL, check_delete_output = TRUE) {
  
  # check before overwriting existing output
  if (project$active_set>0 & check_delete_output) {
    
    # ask before overwriting
    user_choice <- user_yes_no("All existing output and parameter sets for this project will be lost. Continue? (Y/N): ")
    
    # on abort, return original project
    if (!user_choice) {
      return(project)
    }
    
    # replace old project with fresh empty version
    project <- mavproject()
  }
  
  # process and perform checks on data
  dat_proc <- process_data(df, ID_col, pop_col, ploidy_col, data_cols, ploidy, missing_data)
  
  # add data to project
  project[["data"]] <- df
  project[["data_processed"]] <- dat_proc
  
  return(project)
}

#------------------------------------------------
# process data

process_data <- function(df, ID_col, pop_col, ploidy_col, data_cols, ploidy, missing_data) {
  
  # checks on input
  assert_that(is.data.frame(df))
  assert_that(is.int(ID_col))
  assert_that(!any(duplicated(c(ID_col, pop_col, ploidy_col, data_cols))))
  
  # get ploidy in final form
  if (is.null(ploidy_col)) {
    if (is.null(ploidy)) {
      ploidy <- 1
      cat("using default value of ploidy = 1")
    }
    assert_that(is.int(ploidy))
    assert_that(length(ploidy)==1)
    assert_that(ploidy>=1)
    assert_that(ploidy>=1)
    assert_that((nrow(df)%%ploidy)==0)
    ploidy <- rep(ploidy, nrow(df)/ploidy)
  } else {
    assert_that(is.int(ploidy_col))
    assert_that(length(ploidy_col)==1)
    assert_that(ploidy_col>=1)
    assert_that(ploidy_col<=ncol(df))
    ploidy_raw <- df[, ploidy_col]
    ploidy <- NULL
    i <- 1
    while (i<nrow(df)) {
      ploidy <- c(ploidy, ploidy_raw[i])
      i <- i + ploidy_raw[i]
    }
  }
  assert_that(nrow(df)==sum(ploidy))
  ind_first_row <- cumsum(ploidy) - ploidy + 1
  
  # get pop in final form
  if (!is.null(pop_col)) {
    assert_that(is.int(pop_col))
    assert_that(length(pop_col)==1)
    assert_that(pop_col>=1)
    assert_that(pop_col<=ncol(df))
    pop_raw <- df[,pop_col]
    pop <- pop_raw[ind_first_row]
    assert_that(all(is.int(pop)))
  } else {
    pop <- NA
  }
  
  # get data
  if (is.null(data_cols)) {
    data_cols <- setdiff(1:ncol(df), c(ID_col, pop_col, ploidy_col))
  }
  assert_that(all(is.int(data_cols)))
  assert_that(!any(duplicated(data_cols)))
  assert_that(all(data_cols>=1))
  assert_that(all(data_cols<=ncol(df)))
  dat <- as.matrix(df[,data_cols])
  assert_that(ncol(dat)>0)
  assert_that(nrow(dat)>0)
  
  # recode to remove redundancy
  for (j in 1:ncol(dat)) {
    dat[,j] <- match(dat[,j], unique(dat[,j][!is.na(dat[,j])]))
  }
  
  # return list
  ret <- list(dat = dat, ind_first_row = ind_first_row, pop = pop, ploidy = ploidy)
  return(ret)
}

#------------------------------------------------
# create new parameter set

new_set <- function(project, name = "(no name)", admix_on = FALSE, K = 1:3, rungs = 11) {
  
  # count current parameter sets and add one
  s <- length(project$parameter_sets) + 1
  
  # make new set active
  project[["active_set"]] <- s
  
  # create new parameter set
  project[["parameter_sets"]][[s]] <- list(name = name,
                                          admix_on = admix_on,
                                          K = K,
                                          rungs = rungs
                                          )
  
  # create new output corresponding to this set
  project[["output"]][[s]] <- list(name = name)
  
  # return
  return(project)
}

#------------------------------------------------
# delete parameter set

delete_set <- function(project, index = NULL, check_delete_output = TRUE) {
  
  # set index to activeSet by default
  index <- define_default(index, project$active_set)
  
  # check inputs
  assert_that(is.int(index))
  assert_that(length(index)==1)
  assert_that(index>=1)
  assert_that(index<=length(project$parameter_sets))
  
  # check before overwriting existing output
  if (project$active_set>0 & check_delete_output) {
    
    # ask before overwriting
    user_choice <- user_yes_no(sprintf("Output for set %s will be deleted. Continue? (Y/N): ", index))
    
    # on abort, return original project
    if (!user_choice) {
      return(project)
    }
  }
  
  # drop chosen parameter set
  project[["parameter_sets"]][[index]] <- NULL
  
  # drop chosen output
  project[["output"]][[index]] <- NULL
  
  # make new final set active
  project[["active_set"]] <- length(project$parameter_sets)
  
  # return
  return(project)
}

#------------------------------------------------
#' Generate scaffolds
#'
#' TODO
#'
#' @param project TODO
#' @param iterations TODO
#' @param n TODO
#' @param coupling_on TODO
#' @param splitmerge_on TODO
#' @param cluster TODO
#'
#' @export
#' @examples
#' # TODO

generate_scaffolds <- function(project, iterations = NULL, n = 10, coupling_on = TRUE, splitmerge_on = TRUE, cluster = NULL) {
  
  # set defaults
  num_cores <- define_default(num_cores, detectCores())
  
  
  
  # return
  return(project)
}

#------------------------------------------------
#' Run MCMC
#'
#' TODO
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
#' m <- run_mcmc(1:5, parallel_on = FALSE)

run_mcmc <- function(x, K = 3, mu_prior_mean = 0, mu_prior_var = 1e3, sigma = 1, burnin =1e2, samples = 1e3, rungs = 10, solve_label_switching_on = TRUE, coupling_on = TRUE, scaffold_on = TRUE, scaffold_n = 10, splitmerge_on = TRUE, parallel_on = TRUE, num_cores = NULL) {
  
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

