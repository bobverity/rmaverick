
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
#' @title Bind data to project
#'
#' @description Bind data to project
#'
#' @details TODO
#' 
#' @param project TODO
#' @param df TODO
#' @param ID_col TODO
#' @param pop_col TODO
#' @param ploidy_col TODO
#' @param data_cols TODO
#' @param ploidy TODO
#' @param missing_data TODO
#' @param name TODO
#' @param check_delete_output TODO
#'
#' @export
#' @examples
#' # TODO

bind_data <- function(project, df, ID_col = 1, pop_col = NULL, ploidy_col = NULL, data_cols = NULL, ploidy = 1, missing_data = -9, name = NULL, check_delete_output = TRUE) {
  
  # check before overwriting existing output
  if (project$active_set>0 & check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no("All existing output and parameter sets for this project will be lost. Continue? (Y/N): ")) {
      return(project)
    }
    
    # replace old project with fresh empty version
    project <- mavproject()
  }
  
  # process and perform checks on data
  dat_proc <- process_data(df, ID_col, pop_col, ploidy_col, data_cols, ploidy, missing_data)
  
  # add data to project
  project$data <- df
  project$data_processed <- dat_proc
  
  return(project)
}

#------------------------------------------------
# process data
# (not exported)
#' @noRd
process_data <- function(df, ID_col, pop_col, ploidy_col, data_cols, ploidy, missing_data) {
  
  # checks on input
  assert_dataframe(df)
  assert_scalar_pos_int(ID_col, zero_allowed = FALSE)
  assert_leq(ID_col, ncol(df))
  assert_noduplicates(c(ID_col, pop_col, ploidy_col, data_cols))
  
  # get ploidy in final form
  if (is.null(ploidy_col)) {
    if (is.null(ploidy)) {
      ploidy <- 1
      cat("using default value of ploidy = 1")
    }
    assert_scalar_pos_int(ploidy)
    assert_that((nrow(df)%%ploidy)==0)
    ploidy <- rep(ploidy, nrow(df)/ploidy)
  } else {
    assert_scalar_pos_int(ploidy_col, zero_allowed = FALSE)
    assert_leq(ploidy_col, ncol(df))
    ploidy_raw <- df[,ploidy_col]
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
    assert_scalar_pos_int(pop_col, zero_allowed = FALSE)
    assert_leq(pop_col, ncol(df))
    pop_raw <- df[,pop_col]
    pop <- pop_raw[ind_first_row]
    assert_pos_int(pop)
  } else {
    pop <- NA
  }
  
  # get genetic data in final form
  if (is.null(data_cols)) {
    data_cols <- setdiff(1:ncol(df), c(ID_col, pop_col, ploidy_col))
  }
  assert_pos_int(data_cols, zero_allowed = FALSE)
  assert_leq(data_cols, ncol(df))
  assert_noduplicates(data_cols)
  dat <- as.matrix(df[,data_cols,drop = FALSE])
  assert_that(ncol(dat)>0)
  assert_that(nrow(dat)>0)
  assert_that(all(apply(dat, 1, is.numeric)))
  dat[dat==missing_data] <- NA
  
  # recode to remove redundancy
  for (j in 1:ncol(dat)) {
    dat[,j] <- match(dat[,j], unique(dat[,j][!is.na(dat[,j])]))
  }
  
  # return list
  ret <- list(dat = dat,
              ind_first_row = ind_first_row,
              pop = pop,
              ploidy = ploidy)
  return(ret)
}

#------------------------------------------------
#' @title Create new parameter set
#'
#' @description Create new parameter set
#'
#' @details TODO
#' 
#' @param project TODO
#' @param name TODO
#' @param admix_on TODO
#' @param K TODO
#' @param rungs TODO
#' 
#' @export
#' @examples
#' # TODO

new_set <- function(project, name = "(no name)", admix_on = FALSE, K = 1:3, rungs = 11) {
  
  # count current parameter sets and add one
  s <- length(project$parameter_sets) + 1
  
  # make new set active
  project$active_set <- s
  
  # create new parameter set
  project$parameter_sets[[s]] <- list(name = name,
                                      admix_on = admix_on,
                                      K = K,
                                      rungs = rungs)
  
  # create new output corresponding to this set
  project$output[[s]] <- list(name = name)
  
  # return
  return(project)
}

#------------------------------------------------
#' @title Delete parameter set
#'
#' @description Delete parameter set
#'
#' @details TODO
#' 
#' @param project TODO
#' @param index TODO
#' @param check_delete_output TODO
#' 
#' @export
#' @examples
#' # TODO

delete_set <- function(project, index = NULL, check_delete_output = TRUE) {
  
  # check inputs
  assert_mavproject(project)
  
  # set index to activeSet by default
  index <- define_default(index, project$active_set)
  
  # further checks
  assert_scalar_pos_int(index)
  assert_leq(index, length(project$parameter_sets))
  assert_scalar_logical(check_delete_output)
  
  # check before overwriting existing output
  if (project$active_set>0 & check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no(sprintf("Output for set %s will be deleted. Continue? (Y/N): ", index))) {
      return(project)
    }
  }
  
  # drop chosen parameter set
  project$parameter_sets[[index]] <- NULL
  
  # drop chosen output
  project$output[[index]] <- NULL
  
  # make new final set active
  project$active_set <- length(project$parameter_sets)
  
  # return
  return(project)
}

#------------------------------------------------
#' @title Generate scaffolds
#'
#' @description Generate scaffolds
#'
#' @details TODO
#' 
#' @param project TODO
#' @param n TODO
#' @param coupling_on TODO
#' @param splitmerge_on TODO
#' @param cluster TODO
#' @param output_console TODO
#' 
#' @export
#' @examples
#' # TODO

generate_scaffolds <- function(project, n = 10, coupling_on = TRUE, splitmerge_on = TRUE, cluster = NULL, output_console = is.null(cluster)) {
  
  # check inputs
  assert_mavproject(project)
  assert_scalar_pos_int(n, zero_allowed = FALSE)
  assert_scalar_logical(coupling_on)
  assert_scalar_logical(splitmerge_on)
  assert_scalar_logical(output_console)
  
  # define argument list
  all_args <- list()
  for (i in 1:n) {
    
    # create progress bar
    pb_scaf <- txtProgressBar(min = 0, max = n, initial = NA, style = 3)
    
    # define arguments
    all_args[[i]] <- list(args_data = project$data_processed, args_params = project$parameter_sets, coupling_on = coupling_on, splitmerge_on = splitmerge_on, output_console = output_console, test_convergence = test_convergence, update_progress = update_progress, pb_scaf = pb_scaf)
  }
  
  # run efficient Rcpp function
  if (!is.null(cluster)) {  # run in parallel
    if (!inherits(cluster, "cluster")) {
      stop("expected a cluster object")
    }
    clusterEvalQ(cluster, library(rmaverick))
    output_raw <- clusterApplyLB(cl = cluster, all_args, generate_scaffolds_cpp)
  } else {  # run in serial
    output_raw <- lapply(all_args, generate_scaffolds_cpp)
  }
  
  
  # return
  return(project)
}

#------------------------------------------------
#' @title Run main MCMC
#'
#' @description Run main MCMC
#'
#' @details TODO
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
#' # TODO

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

