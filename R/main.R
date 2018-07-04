
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib rmaverick
#' @import assertthat
#' @import parallel
#' @import coda
#' @importFrom Rcpp evalCpp
#' @import graphics
#' @import stats
#' @import utils
#' @importFrom grDevices colorRampPalette
NULL

#------------------------------------------------
#' @title Bind data to project
#'
#' @description Load data into a \code{mavproject} prior to analysis. Data must be formatted as a dataframe with samples in rows and loci in columns. If individuals are polyploid then multiple rows can be used per sample. Ploidy is allowed to vary between samples, and can be specified in multiple ways (see details below).
#'
#' @details TODO
#' 
#' @param project an rmaverick project, as produced by the function
#'   \code{mavproject()}
#' @param df a dataframe containing genetic information and optional meta-data
#' @param ID_col which column of the input data contains the sample IDs. If NULL
#'   then IDs must be defined seperately through the \code{ID} argument
#' @param pop_col which column of the input data contains the ostensible
#'   population of the samples. If NULL then populations must be defined
#'   seperately through the \code{pop} argument
#' @param ploidy_col which column of the input data contains the ploidy of the
#'   samples. If NULL then ploidy must be defined seperately through the
#'   \code{ploidy} argument
#' @param data_cols which columns of the input data contain genetic information.
#'   Defaults to all remaining columns of the data once special columns have
#'   been accounted for
#' @param ID sample IDs, if not using the \code{ID_col} option
#' @param pop ostensible populations, if not using the \code{pop_col} option
#' @param ploidy sample ploidy, if not using the \code{ploidy_col} option. Can
#'   be a scalar, in which case the same value is used for all samples, or a
#'   vector specifying the ploidy seperately for each sample.
#' @param missing_data which value represents missing data
#' @param name optional name of the data set, to aid in record keeping
#' @param check_delete_output whether to prompt the user before overwriting 
#'   existing data
#'
#' @export
#' @examples
#' # TODO

bind_data <- function(project, df, ID_col = 1, pop_col = NULL, ploidy_col = NULL, data_cols = NULL, ID = NULL, pop = NULL, ploidy = NULL, missing_data = -9, name = NULL, check_delete_output = TRUE) {
  
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
  dat_proc <- process_data(df, ID_col, pop_col, ploidy_col, data_cols, ID, pop, ploidy, missing_data)
  
  # add data to project
  project$data <- df
  project$data_processed <- dat_proc
  
  return(project)
}

#------------------------------------------------
# process data
# (not exported)
#' @noRd
process_data <- function(df, ID_col, pop_col, ploidy_col, data_cols, ID, pop, ploidy, missing_data) {
  
  # checks on input
  assert_dataframe(df)
  assert_noduplicates(c(ID_col, pop_col, ploidy_col, data_cols))
  
  # get ploidy in final form
  if (is.null(ploidy_col)) {
    if (is.null(ploidy)) {
      ploidy <- 1
      message("using default value of ploidy = 1")
    }
    if (length(ploidy)==1) {
      assert_that((nrow(df)%%ploidy)==0)
      ploidy <- rep(ploidy, nrow(df)/ploidy)
    }
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
  assert_pos_int(ploidy)
  assert_that(nrow(df)==sum(ploidy))
  ind_first_row <- cumsum(ploidy) - ploidy + 1
  n <- length(ploidy)
  
  # get sample IDs in final form
  if (is.null(ID_col)) {
    ID <- define_default(ID, paste0("sample", 1:n))
  } else {
    assert_scalar_pos_int(ID_col, zero_allowed = FALSE)
    assert_leq(ID_col, ncol(df))
    ID <- df[ind_first_row, ID_col]
  }
  assert_length(ID,n)
  
  # get pop in final form
  if (is.null(pop_col)) {
    pop <- define_default(pop, rep(1,n))
  } else {
    assert_scalar_pos_int(pop_col, zero_allowed = FALSE)
    assert_leq(pop_col, ncol(df))
    pop <- df[ind_first_row,pop_col]
  }
  assert_pos_int(pop)
  assert_length(pop,n)
  
  # get genetic data in final form
  if (is.null(data_cols)) {
    data_cols <- setdiff(1:ncol(df), c(ID_col, pop_col, ploidy_col))
  }
  assert_pos_int(data_cols, zero_allowed = FALSE)
  assert_leq(data_cols, ncol(df))
  assert_noduplicates(data_cols)
  
  dat <- as.matrix(df[,data_cols,drop = FALSE])
  L <- ncol(dat)
  assert_that(n>0)
  assert_that(L>0)
  assert_that(all(apply(dat, 1, is.numeric)))
  dat[dat==missing_data] <- NA
  
  # recode to remove redundancy
  Jl <- rep(NA, L)
  for (j in 1:L) {
    u <- unique(dat[,j][!is.na(dat[,j])])
    Jl[j] <- length(u)
    dat[,j] <- match(dat[,j], u)
  }
  
  # return list
  ret <- list(dat = dat,
              ID = ID,
              Jl = Jl,
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
#' @param project an rmaverick project, as produced by the function
#'   \code{mavproject()}
#' @param name TODO
#' @param lambda TODO
#' @param admix_on TODO
#' @param alpha TODO
#' @param estimate_alpha TODO
#' 
#' @export
#' @examples
#' # TODO

new_set <- function(project, name = "(no name)", lambda = 1.0, admix_on = FALSE, alpha = 0.1, estimate_alpha = TRUE) {
  
  # TODO - check inputs
  
  # count current parameter sets and add one
  s <- length(project$parameter_sets) + 1
  
  # make new set active
  project$active_set <- s
  
  # create new parameter set
  project$parameter_sets[[s]] <- list(name = name,
                                      lambda = lambda,
                                      admix_on = admix_on,
                                      alpha = alpha,
                                      estimate_alpha = estimate_alpha)
  
  # create new output corresponding to this set
  project$output$single_set[[s]] <- list(single_K = list(),
                                         all_K = list())
  names(project$output$single_set) <- paste0("set", 1:length(project$output$single_set))
  
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
#' @param project an rmaverick project, as produced by the function
#'   \code{mavproject()}
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
  project$output$single_set[[index]] <- NULL
  
  # TODO - recalculate evidence over sets
  
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
#' @param project an rmaverick project, as produced by the function
#'   \code{mavproject()}
#' @param n TODO
#' @param coupling_on TODO
#' @param splitmerge_on TODO
#' @param cluster TODO
#' @param silent TODO
#' 
#' @export
#' @examples
#' # TODO

generate_scaffolds <- function(project, n = 10, coupling_on = TRUE, splitmerge_on = TRUE, cluster = NULL, silent = !is.null(cluster)) {
  
  # TODO - define K as argument
  
  # check inputs
  #assert_mavproject(project)
  #assert_scalar_pos_int(n, zero_allowed = FALSE)
  #assert_scalar_logical(coupling_on)
  #assert_scalar_logical(splitmerge_on)
  #if (!is.null(cluster)) {
  #  assert_cluster(cluster)
  #}
  #assert_scalar_logical(silent)
  
  # define argument list
  #all_args <- list()
  #for (i in 1:n) {
    
    # create progress bar
    #pb_scaf <- txtProgressBar(min = 0, max = n, initial = NA, style = 3)
    
    # define arguments
    #all_args[[i]] <- list(args_data = project$data_processed, args_params = project$parameter_sets, coupling_on = coupling_on, splitmerge_on = splitmerge_on, output_console = output_console, test_convergence = test_convergence, update_progress = update_progress, pb_scaf = pb_scaf)
  #}
  
  # run efficient Rcpp function
  #if (!is.null(cluster)) {  # run in parallel
  #  clusterEvalQ(cluster, library(rmaverick))
  #  output_raw <- clusterApplyLB(cl = cluster, all_args, generate_scaffolds_cpp)
  #} else {  # run in serial
  #  output_raw <- lapply(all_args, generate_scaffolds_cpp)
  #}
  
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
#' @param project an rmaverick project, as produced by the function
#'   \code{mavproject()}
#' @param K the values of K to run the MCMC on
#' @param burnin the number of burn-in iterations
#' @param samples the number of sampling iterations
#' @param rungs the number of temperature rungs
#' @param GTI_pow the power used in generalised thermodynamic integration
#' @param auto_converge whether convergence should be assessed automatically, in which case burn-in iterations are terminated as soon as convergence is reached. Otherwise the full \code{burnin} iterations are used
#' @param solve_label_switching_on whether to implement the Stevens' solution to the label-switching problem
#' @param coupling_on whether to implement Metropolis-coupling over temperature rungs
#' @param scaffold_on whether to use scaffolds to improve mixing
#' @param splitmerge_on whether to implement a split-merge proposal
#' @param cluster pass in a cluster environment
#' @param pb_markdown whether to run progress bars in markdown mode, in which
#'   case they are updated once at the end to avoid large amounts of output.
#' @param silent whether to suppress all console output
#' 
#' @export
#' @examples
#' # TODO

run_mcmc <- function(project, K = 3, burnin = 1e2, samples = 1e3, rungs = 10, GTI_pow = 3, auto_converge = TRUE, solve_label_switching_on = TRUE, coupling_on = TRUE, scaffold_on = TRUE, splitmerge_on = TRUE, cluster = NULL, pb_markdown = FALSE, silent = !is.null(cluster)) {
  
  # start timer
  t0 <- Sys.time()
  
  # check inputs
  assert_mavproject(project)
  assert_pos_int(K, zero_allowed = FALSE)
  assert_scalar_pos_int(burnin, zero_allowed = FALSE)
  assert_scalar_pos_int(samples, zero_allowed = FALSE)
  assert_scalar_pos_int(rungs, zero_allowed = FALSE)
  assert_scalar_pos(GTI_pow)
  assert_bounded(GTI_pow, 1.5, 10)
  assert_scalar_logical(solve_label_switching_on)
  assert_scalar_logical(coupling_on)
  assert_scalar_logical(scaffold_on)
  assert_scalar_logical(splitmerge_on)
  if (!is.null(cluster)) {
    assert_cluster(cluster)
  }
  assert_scalar_logical(pb_markdown)
  assert_scalar_logical(silent)
  
  # get active set
  s <- project$active_set
  
  # get useful quantities
  ploidy <- project$data_processed$ploidy
  Jl <- project$data_processed$Jl
  n <- length(ploidy)
  L <- length(Jl)
  admix_on <- project$parameter_sets[[s]]$admix_on
  
  # ---------- create argument lists ----------
  
  # data list
  args_data <- list(data = mat_to_rcpp(project$data_processed$dat),
                    Jl = project$data_processed$Jl,
                    ploidy = project$data_processed$ploidy)
  
  # input arguments list
  args_inputs <- list(burnin = burnin,
                      samples = samples,
                      rungs = rungs,
                      GTI_pow = GTI_pow,
                      auto_converge = auto_converge,
                      solve_label_switching_on = solve_label_switching_on,
                      coupling_on = coupling_on,
                      scaffold_on = scaffold_on,
                      splitmerge_on = splitmerge_on,
                      pb_markdown = pb_markdown,
                      silent = silent)
  
  # combine model parameters list with input arguments
  args_model <- c(project$parameter_sets[[s]], args_inputs)
  
  # R functions to pass to Rcpp
  args_functions <- list(test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # define final argument list over all K>1
  parallel_args <- list()
  for (i in 1:length(K)) {
    
    # create progress bars
    pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
    pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
    args_progress <- list(pb_burnin = pb_burnin,
                          pb_samples = pb_samples)
    
    # incporporate arguments unique to this K
    args_model$K <- K[i]
    
    # create argument list
    parallel_args[[i]] <- list(args_data = args_data,
                               args_model = args_model,
                               args_functions = args_functions,
                               args_progress = args_progress)
  }
  
  # ---------- run MCMC ----------
  
  # split into parallel and serial implementations
  if (!is.null(cluster)) { # run in parallel
    clusterEvalQ(cluster, library(rmaverick))
    output_raw <- clusterApplyLB(cl = cluster, parallel_args, run_mcmc_cpp)
  } else { # run in serial
    output_raw <- lapply(parallel_args, run_mcmc_cpp)
  }
  
  #------------------------
  
  # begin processing results
  if (!silent) {
    cat("Processing results\n")
  }
  
  # loop through K
  ret <- list()
  for (i in 1:length(K)) {
    
    # create name lists
    ind_names <- paste0("ind", 1:n)
    locus_names <- paste0("locus", 1:L)
    deme_names <- paste0("deme", 1:K[i])
    rung_names <- paste0("rung", 1:rungs)
    
    # define output manually if K==1
    if (K[i]==1) {
      
      # create qmatrix_ind
      qmatrix_ind <- matrix(1,n,1)
      colnames(qmatrix_ind) <- deme_names
      rownames(qmatrix_ind) <- ind_names
      class(qmatrix_ind) <- "maverick_qmatrix_ind"
      
      # create NULL outputs
      loglike_burnin <- NULL
      loglike_sampling <- NULL
      loglike_quantiles <- NULL
      alpha <- NULL
      ESS <- NULL
      GTI_path <- NULL
      GTI_logevidence <- NULL
      coupling_accept <- NULL
      
    } else { # extract output if K>1
      
      # ---------- raw mcmc results ----------
      
      # get loglikelihood in coda::mcmc format
      loglike_burnin <- mcmc(t(rcpp_to_mat(output_raw[[i]]$loglike_burnin)))
      loglike_sampling <- mcmc(t(rcpp_to_mat(output_raw[[i]]$loglike_sampling)))
      
      # alpha
      alpha <- NULL
      if (admix_on) {
        alpha <- mcmc(output_raw[[i]]$alpha_store)
      }
      
      # ---------- summary results ----------
      
      # get quantiles over sampling loglikelihoods
      loglike_quantiles <- t(apply(loglike_sampling, 2, quantile_95))
      rownames(loglike_quantiles) <- rung_names
      class(loglike_quantiles) <- "maverick_loglike_quantiles"
      
      # process qmatrix_ind
      qmatrix_ind <- rcpp_to_mat(output_raw[[i]]$qmatrix_ind)
      colnames(qmatrix_ind) <- deme_names
      rownames(qmatrix_ind) <- ind_names
      class(qmatrix_ind) <- "maverick_qmatrix_ind"
      
      # ---------- GTI path and model evidence ----------
      
      # get ESS
      ESS <- effectiveSize(loglike_sampling)
      ESS[ESS==0] <- samples # if no variation then assume perfect mixing
      names(ESS) <- rung_names
      
      # weight likelihood according to GTI_pow
      loglike_weighted <- loglike_sampling
      for (j in 1:rungs) {
        beta_j <- j/rungs
        loglike_weighted[,j] <- GTI_pow*beta_j^(GTI_pow-1) * loglike_sampling[,j]
      }
      
      # calculate GTI path mean and SE
      GTI_path_mean <- colMeans(loglike_weighted)
      GTI_path_var <- apply(loglike_weighted, 2, var)
      GTI_path_SE <- sqrt(GTI_path_var/ESS)
      GTI_path <- data.frame(mean = GTI_path_mean, SE = GTI_path_SE)
      rownames(GTI_path) <- rung_names
      class(GTI_path) <- "maverick_GTI_path"
      
      # calculate GTI estimate of log-evidence
      GTI_vec <- 0.5*loglike_weighted[,1]/rungs
      if (rungs>1) {
        for (j in 2:rungs) {
          GTI_vec <- GTI_vec + 0.5*(loglike_weighted[,j]+loglike_weighted[,j-1])/rungs
        }
      }
      GTI_logevidence_mean <- mean(GTI_vec)
      
      # calculate standard error of GTI estimate
      GTI_ESS <- as.numeric(effectiveSize(GTI_vec))
      if (GTI_ESS==0) {
        GTI_ESS <- samples # if no variation then assume perfect mixing
      }
      GTI_logevidence_SE <- sqrt(var(GTI_vec)/GTI_ESS)
      
      # produce final GTI_logevidence object
      GTI_logevidence <- data.frame(estimate = GTI_logevidence_mean,
                                    SE = GTI_logevidence_SE)
      
      # ---------- acceptance rates ----------
      
      # process acceptance rates
      coupling_accept <- output_raw[[i]]$coupling_accept/samples
      
    }
    
    # ---------- save arguments ----------
    
    output_args <- list(burnin = burnin,
                        samples = samples,
                        rungs = rungs,
                        GTI_pow = GTI_pow,
                        auto_converge = auto_converge,
                        coupling_on = coupling_on,
                        scaffold_on = scaffold_on,
                        splitmerge_on = splitmerge_on,
                        solve_label_switching_on = solve_label_switching_on,
                        silent = silent)
    
    # ---------- save results ----------
    
    # add to project
    project$output$single_set[[s]]$single_K[[K[i]]] <- list()
    
    project$output$single_set[[s]]$single_K[[K[i]]]$summary <- list(qmatrix_ind = qmatrix_ind,
                                                                    loglike_quantiles = loglike_quantiles,
                                                                    ESS = ESS,
                                                                    GTI_path = GTI_path,
                                                                    GTI_logevidence = GTI_logevidence)
    
    project$output$single_set[[s]]$single_K[[K[i]]]$raw <- list(loglike_burnin = loglike_burnin,
                                                                loglike_sampling = loglike_sampling,
                                                                alpha = alpha,
                                                                coupling_accept = coupling_accept)
    
    project$output$single_set[[s]]$single_K[[K[i]]]$function_call <- list(args = output_args,
                                                                          call = match.call())
    
  }
  
  # name output over K
  K_all <- length(project$output$single_set[[s]]$single_K)
  names(project$output$single_set[[s]]$single_K) <- paste0("K", 1:K_all)
  
  # ---------- recalculate evidence ----------
  
  # get logevidence over all K
  GTI_logevidence_raw <- mapply(function(x) {
    GTI_logevidence <- x$summary$GTI_logevidence
    if (is.null(GTI_logevidence)) {
      return(rep(NA,2))
    } else {
      return(unlist(GTI_logevidence))
    }
  }, project$output$single_set[[s]]$single_K)
  GTI_logevidence <- as.data.frame(t(GTI_logevidence_raw))
  names(GTI_logevidence) <- c("mean", "SE")
  rownames(GTI_logevidence) <- NULL
  GTI_logevidence <- cbind(K = 1:nrow(GTI_logevidence), GTI_logevidence)
  
  # load GTI_logevidence into project
  project$output$single_set[[s]]$all_K$GTI_logevidence <- GTI_logevidence
  
  # obtain posterior K estimates
  if (any(!is.na(GTI_logevidence$mean))) {
    
    # produce posterior estimates by simulation
    w <- which(!is.na(GTI_logevidence$mean))
    GTI_posterior_raw <- GTI_posterior_K_sim_cpp(list(mean = GTI_logevidence$mean[w],
                                                      SE = GTI_logevidence$SE[w],
                                                      reps = 1e6))$ret
    
    # get posterior quantiles in dataframe
    GTI_posterior_quantiles <- t(mapply(quantile_95, GTI_posterior_raw))
    GTI_posterior <- data.frame(K = 1:nrow(GTI_logevidence),
                                Q2.5 = NA,
                                Q50 = NA,
                                Q97.5 = NA)
    GTI_posterior[w,-1] <- GTI_posterior_quantiles
    
    # load GTI_posterior into project
    project$output$single_set[[s]]$all_K$GTI_posterior <- GTI_posterior
  }
  
  # end timer
  tdiff <- round(Sys.time() - t0, 2)
  message(sprintf("Total run-time: %s seconds", tdiff))
  
  # return invisibly
  invisible(project)
}

