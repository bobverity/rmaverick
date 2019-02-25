
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib rmaverick
#' @import parallel
#' @import coda
#' @import ggplot2
#' @import gridExtra
#' @importFrom Rcpp evalCpp
#' @import graphics
#' @import stats
#' @import utils
#' @importFrom grDevices colorRampPalette
NULL

#------------------------------------------------
#' @title Check that rmaverick package has loaded successfully
#'
#' @description Simple function to check that rmaverick package has loaded
#'   successfully. Prints "rmaverick loaded successfully!" if so.
#'
#' @export

check_rmaverick_loaded <- function() {
  message("rmaverick loaded successfully!")
}

#------------------------------------------------
#' @title Bind data to project
#'   
#' @description Load data into a \code{mavproject} prior to analysis. Data must 
#'   be formatted as a dataframe with samples in rows and loci in columns. If 
#'   individuals are polyploid then multiple rows can be used per sample. Ploidy
#'   is allowed to vary between samples, and can be specified in multiple ways.
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
#' @param wide_format if \code{TRUE} then uses one line per sample, with loci
#'   stacked side-by-side in columns. When using this format the ploidy must be
#'   the same for all samples, and must be specified using the \code{ploidy}
#'   variable rather than as a seperate column
#' @param name optional name of the data set to aid in record keeping
#' @param check_delete_output whether to prompt the user before overwriting 
#'   existing data
#'
#' @export

bind_data <- function(project, df, ID_col = 1, pop_col = NULL, ploidy_col = NULL, data_cols = NULL, ID = NULL, pop = NULL, ploidy = NULL, missing_data = -9, wide_format = FALSE, name = NULL, check_delete_output = TRUE) {
  
  # check inputs (further checks carried out in process_data() function)
  assert_custom_class(project, "mavproject")
  assert_dataframe(df)
  assert_noduplicates(c(ID_col, pop_col, ploidy_col, data_cols))
  
  # check before overwriting existing output
  if (project$active_set>0 && check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no("All existing output and parameter sets for this project will be lost. Continue? (Y/N): ")) {
      return(project)
    }
    
    # replace old project with fresh empty version
    project <- mavproject()
  }
  
  # process and perform checks on data
  dat_processed <- process_data(df, ID_col, pop_col, ploidy_col, data_cols, ID, pop, ploidy, missing_data, wide_format)
  dat_processed$name <- name
  
  # add data to project
  project$data <- df
  project$data_processed <- dat_processed
  
  return(project)
}

#------------------------------------------------
# process data
#' @noRd
process_data <- function(df, ID_col, pop_col, ploidy_col, data_cols, ID, pop, ploidy, missing_data, wide_format) {
  
  # process differently in wide vs. long format
  if (wide_format) {
    ret <- process_data_wide(df = df,
                             ID_col = ID_col,
                             pop_col = pop_col,
                             ploidy_col = ploidy_col,
                             data_cols = data_cols,
                             ID = ID,
                             pop = pop,
                             ploidy = ploidy,
                             missing_data = missing_data)
  } else {
    ret <- process_data_long(df = df,
                             ID_col = ID_col,
                             pop_col = pop_col,
                             ploidy_col = ploidy_col,
                             data_cols = data_cols,
                             ID = ID,
                             pop = pop,
                             ploidy = ploidy,
                             missing_data = missing_data)
  }
  
  return(ret)
}

#------------------------------------------------
# process data in long format
#' @noRd
process_data_long <- function(df, ID_col, pop_col, ploidy_col, data_cols, ID, pop, ploidy, missing_data) {
  
  # get ploidy in final form
  if (is.null(ploidy_col)) {
    if (is.null(ploidy)) {
      message("using default value of ploidy = 1")
      ploidy <- 1
    }
    if (length(ploidy) == 1) {
      assert_eq((nrow(df)%%ploidy), 0)
      ploidy <- rep(ploidy, nrow(df)/ploidy)
    }
  } else {
    assert_single_pos_int(ploidy_col, zero_allowed = FALSE)
    assert_leq(ploidy_col, ncol(df))
    ploidy_raw <- df[,ploidy_col]
    ploidy <- NULL
    i <- 1
    while (i <= nrow(df)) {
      ploidy <- c(ploidy, ploidy_raw[i])
      i <- i + ploidy_raw[i]
    }
  }
  assert_pos_int(ploidy)
  assert_nrow(df, sum(ploidy))
  ind_first_row <- cumsum(ploidy) - ploidy + 1
  n <- length(ploidy)
  
  # get sample IDs in final form
  if (is.null(ID_col)) {
    ID <- define_default(ID, paste0("sample", 1:n))
  } else {
    assert_single_pos_int(ID_col, zero_allowed = FALSE)
    assert_leq(ID_col, ncol(df))
    ID <- df[ind_first_row, ID_col]
  }
  assert_length(ID,n)
  
  # get pop in final form
  if (is.null(pop_col)) {
    pop <- define_default(pop, rep(1,n))
  } else {
    assert_single_pos_int(pop_col, zero_allowed = FALSE)
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
  apply(dat, 1, assert_numeric)
  dat[dat == missing_data] <- NA
  
  # recode to remove redundancy
  Jl <- rep(NA, L)
  for (j in 1:L) {
    u <- unique(dat[,j][!is.na(dat[,j])])
    Jl[j] <- length(u)
    dat[,j] <- match(dat[,j], u)
  }
  dat[is.na(dat)] <- 0
  
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
# process data in wide format
#' @noRd
process_data_wide <- function(df, ID_col, pop_col, ploidy_col, data_cols, ID, pop, ploidy, missing_data) {
  
  # check inputs
  assert_null(ploidy_col, message = "ploidy_col must be null when using wide_format")
  assert_non_null(ploidy, message = "ploidy must be specified when using wide format")
  assert_single_pos_int(ploidy)
  
  # get genetic data columns
  if (is.null(data_cols)) {
    data_cols <- setdiff(1:ncol(df), c(ID_col, pop_col, ploidy_col))
  }
  assert_pos_int(data_cols, zero_allowed = FALSE)
  assert_leq(data_cols, ncol(df))
  assert_noduplicates(data_cols)
  
  # check ploidy compatible with data dimensions
  assert_eq((length(data_cols)%%ploidy), 0)
  
  # get genetic data into long format
  dat <- as.matrix(df[,data_cols,drop = FALSE])
  n <- nrow(dat)
  L <- ncol(dat)/ploidy
  tmp <- apply(dat, 1, function(x){matrix(x, ncol = ploidy, byrow = TRUE)})
  dat_long <- matrix(tmp, ncol = L, byrow = TRUE)
  
  # add meta-data columns back in
  if (is.null(ID_col)) {
    ID <- define_default(ID, paste0("sample", 1:n))
  } else {
    assert_single_pos_int(ID_col, zero_allowed = FALSE)
    assert_leq(ID_col, ncol(df))
    ID <- df[, ID_col]
  }
  assert_length(ID,n)
  df_out <- data.frame(ID = rep(ID, each = ploidy), stringsAsFactors = FALSE)
  if (!is.null(pop_col)) {
    assert_single_pos_int(pop_col, zero_allowed = FALSE)
    assert_leq(pop_col, ncol(df))
    df_out <- cbind(df_out, pop = rep(df[,pop_col], each = ploidy))
  }
  df_out <- cbind(df_out, dat_long)
  
  # now can process data in long format
  ret <- process_data_long(df = df_out,
                           ID_col = ID_col,
                           pop_col = pop_col,
                           ploidy_col = ploidy_col,
                           data_cols = NULL,
                           ID = ID,
                           pop = pop,
                           ploidy = ploidy,
                           missing_data = missing_data)
  
  return(ret)
}

#------------------------------------------------
#' @title Create new parameter set
#'   
#' @description Create a new parameter set within an rmaverick project. The new 
#'   parameter set becomes the active set once created.
#'   
#' @param project an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param name the name of the parameter set
#' @param lambda shape parameter of the dirichlet prior on allele frequencies in
#'   each subpopulation
#' @param admix_on whether to allow admixture between subpopulations, in which
#'   case each gene copy can be assigned to a different subpopulation
#' @param alpha parameter governing the strength of admixture. Higher values
#'   give greater probability of a gene copy being assigned to different
#'   subpopulations, while as alpha tends to zero we converge back at the
#'   no-admixture model in which all gene copies within an individual are
#'   constrained to have originated from the same subpopulation
#' @param estimate_alpha whether the value of alpha should be estimated as part
#'   of the MCMC, in which case the value \code{alpha} is ignored
#' 
#' @export

new_set <- function(project, name = "(no name)", lambda = 1.0, admix_on = FALSE, alpha = 0.1, estimate_alpha = TRUE) {
  
  # check inputs
  assert_custom_class(project, "mavproject")
  assert_string(name)
  assert_single_pos(lambda)
  assert_single_logical(admix_on)
  assert_single_pos(alpha)
  assert_single_logical(estimate_alpha)
  
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
  
  names(project$parameter_sets)[s] <- paste0("set", s)
  
  # create new output corresponding to this set
  GTI_logevidence <- data.frame(K = numeric(),
                                mean = numeric(),
                                SE = numeric())
  class(GTI_logevidence) <- "maverick_GTI_logevidence"
  
  GTI_posterior <- data.frame(K = numeric(),
                              Q2.5 = numeric(),
                              Q50 = numeric(),
                              Q97.5 = numeric())
  class(GTI_posterior) <- "maverick_GTI_posterior"
  
  project$output$single_set[[s]] <- list(single_K = list(),
                                         all_K = list(GTI_logevidence = GTI_logevidence,
                                                      GTI_posterior = GTI_posterior))
  
  names(project$output$single_set) <- paste0("set", 1:length(project$output$single_set))
  
  # expand summary output over all parameter sets
  GTI_logevidence_model <- rbind(project$output$all_sets$GTI_logevidence_model, data.frame(set = s, name = name, mean = NA, SE = NA, stringsAsFactors = FALSE))
  class(GTI_logevidence_model) <- "maverick_GTI_logevidence_model"
  project$output$all_sets$GTI_logevidence_model <- GTI_logevidence_model
  
  GTI_posterior_model <- rbind(project$output$all_sets$GTI_posterior_model, data.frame(set = s, name = name, Q2.5 = NA, Q50 = NA, Q97.5 = NA, stringsAsFactors = FALSE))
  class(GTI_posterior_model) <- "maverick_GTI_posterior_model"
  project$output$all_sets$GTI_posterior_model <- GTI_posterior_model
  
  # return
  return(project)
}

#------------------------------------------------
#' @title Delete parameter set
#'   
#' @description Delete a given parameter set from an rmaverick project.
#'   
#' @param project an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param set which set to delete. Defaults to the current active set
#' @param check_delete_output whether to prompt the user before deleting any
#'   existing output
#'   
#' @export

delete_set <- function(project, set = NULL, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(project, "mavproject")
  assert_single_logical(check_delete_output)
  
  # set index to active_set by default
  set <- define_default(set, project$active_set)
  
  # further checks
  assert_single_pos_int(set, zero_allowed = FALSE)
  assert_leq(set, length(project$parameter_sets))
  
  # check before overwriting existing output
  if (project$active_set>0 & check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no(sprintf("Any existing output for set %s will be deleted. Continue? (Y/N): ", set))) {
      return(project)
    }
  }
  
  # drop chosen parameter set
  project$parameter_sets[[set]] <- NULL
  
  # drop chosen output
  project$output$single_set[[set]] <- NULL
  
  GTI_logevidence_model <- as.data.frame(unclass(project$output$all_sets$GTI_logevidence_model))[-set,]
  class(GTI_logevidence_model) <- "maverick_GTI_logevidence_model"
  project$output$all_sets$GTI_logevidence_model <- GTI_logevidence_model
  
  GTI_posterior_model <- as.data.frame(unclass(project$output$all_sets$GTI_posterior_model))[-set,]
  class(GTI_posterior_model) <- "maverick_GTI_posterior_model"
  project$output$all_sets$GTI_posterior_model <- GTI_posterior_model
  
  # make new final set active
  project$active_set <- length(project$parameter_sets)
  
  # recalculate evidence over sets if needed
  if (project$active_set>0) {
    project <- recalculate_evidence(project)
  }
  
  # return
  return(project)
}

#------------------------------------------------
#' @title Change active parameter set
#'   
#' @description Change the active parameter set within an rmaverick project.
#'   
#' @param project an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param set which set to make the new active set
#'   
#' @export

change_set <- function(project, set) {
  
  # check inputs
  assert_custom_class(project, "mavproject")
  assert_single_pos_int(set)
  assert_leq(set, length(project$parameter_sets))
  
  # change active set
  project$active_set <- set
  
  # return
  return(project)
}

#------------------------------------------------
#' @title Run main MCMC
#'   
#' @description Run the main rmaverick MCMC. Model parameters are taken from the
#'   current active parameter set, and MCMC parameters are passed in as
#'   arguments. All output is stored within the project.
#'   
#' @param project an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K the values of K that the MCMC will explore
#' @param burnin the number of burn-in iterations
#' @param samples the number of sampling iterations
#' @param rungs the number of temperature rungs
#' @param GTI_pow the power used in the generalised thermodynamic integration 
#'   method. Must be greater than 1.1
#' @param auto_converge whether convergence should be assessed automatically 
#'   every \code{converge_test} iterations, leading to termination of the 
#'   burn-in phase. If \code{FALSE} then the full \code{burnin} iterations are 
#'   used
#' @param converge_test test for convergence every \code{convergence_test} 
#'   iterations if \code{auto_converge} is being used
#' @param solve_label_switching_on whether to implement the Stevens' solution to
#'   the label-switching problem. If turned off then Q-matrix output will no 
#'   longer be correct, although evidence estimates will be unaffected.
#' @param coupling_on whether to implement Metropolis-coupling over temperature 
#'   rungs
#' @param cluster option to pass in a cluster environment (see package 
#'   "parallel")
#' @param pb_markdown whether to run progress bars in markdown mode, in which 
#'   case they are updated once at the end to avoid large amounts of output.
#' @param silent whether to suppress all console output
#' 
#' @export

run_mcmc <- function(project, K = 3, burnin = 1e2, samples = 1e3, rungs = 10, GTI_pow = 2, auto_converge = TRUE, converge_test = ceiling(burnin/10), solve_label_switching_on = TRUE, coupling_on = TRUE, cluster = NULL, pb_markdown = FALSE, silent = FALSE) {
  
  # start timer
  t0 <- Sys.time()
  
  # check inputs
  assert_custom_class(project, "mavproject")
  assert_pos_int(K, zero_allowed = FALSE)
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(rungs, zero_allowed = FALSE)
  assert_single_pos(GTI_pow)
  assert_gr(GTI_pow, 1.1)
  assert_single_logical(auto_converge)
  assert_single_pos_int(converge_test, zero_allowed = FALSE)
  assert_single_logical(solve_label_switching_on)
  assert_single_logical(coupling_on)
  if (!is.null(cluster)) {
    assert_cluster(cluster)
  }
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # get active set
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
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
                      converge_test = converge_test,
                      solve_label_switching_on = solve_label_switching_on,
                      coupling_on = coupling_on,
                      pb_markdown = pb_markdown,
                      silent = !is.null(cluster))
  
  # combine model parameters list with input arguments
  args_model <- c(project$parameter_sets[[s]], args_inputs)
  
  # R functions to pass to Rcpp
  args_functions <- list(test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # define final argument list over all K
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
  all_converged <- TRUE
  for (i in 1:length(K)) {
    
    # create name lists
    ind_names <- paste0("ind", 1:n)
    locus_names <- paste0("locus", 1:L)
    deme_names <- paste0("deme", 1:K[i])
    rung_names <- paste0("rung", 1:rungs)
    
    # define output manually if K==1
    if (K[i]==1) {
      
      # get exact log-likelihood
      exact_loglike <- output_raw[[i]]$loglike_sampling[[1]][1]
      
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
      GTI_logevidence <- data.frame(estimate = exact_loglike,
                                    SE = 0)
      coupling_accept <- NULL
      
      # all rungs converged by definition
      converged <- rep(TRUE, rungs)
      
      # get run time
      run_time <- output_raw[[i]]$run_time
      
    } else { # extract output if K>1
      
      # ---------- raw mcmc results ----------
      
      # get loglikelihood in coda::mcmc format
      loglike_burnin <- mapply(function(x){mcmc(x)}, output_raw[[i]]$loglike_burnin)
      loglike_sampling <- mcmc(t(rcpp_to_mat(output_raw[[i]]$loglike_sampling)))
      
      # alpha
      alpha <- NULL
      if (admix_on) {
        alpha <- mcmc(output_raw[[i]]$alpha_store)
      }
      
      # get whether rungs have converged
      converged <- output_raw[[i]]$rung_converged
      if (all_converged && any(!converged)) {
        all_converged <- FALSE
      }
      
      # get run time
      run_time <- output_raw[[i]]$run_time
      
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
      ESS[ESS == 0] <- samples # if no variation then assume zero autocorrelation
      ESS[ESS > samples] <- samples # ESS cannot exceed actual number of samples taken
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
                        converge_test = converge_test,
                        solve_label_switching_on = solve_label_switching_on,
                        coupling_on = coupling_on,
                        pb_markdown = pb_markdown,
                        silent = silent)
    
    # ---------- save results ----------
    
    # add to project
    project$output$single_set[[s]]$single_K[[K[i]]] <- list()
    
    project$output$single_set[[s]]$single_K[[K[i]]]$summary <- list(qmatrix_ind = qmatrix_ind,
                                                                    loglike_quantiles = loglike_quantiles,
                                                                    ESS = ESS,
                                                                    GTI_path = GTI_path,
                                                                    GTI_logevidence = GTI_logevidence,
                                                                    converged = converged,
                                                                    run_time = run_time)
    
    project$output$single_set[[s]]$single_K[[K[i]]]$raw <- list(loglike_burnin = loglike_burnin,
                                                                loglike_sampling = loglike_sampling,
                                                                alpha = alpha,
                                                                coupling_accept = coupling_accept)
    
    project$output$single_set[[s]]$single_K[[K[i]]]$function_call <- list(args = output_args,
                                                                          call = match.call())
    
  } # end loop over K
  
  # name output over K
  K_all <- length(project$output$single_set[[s]]$single_K)
  names(project$output$single_set[[s]]$single_K) <- paste0("K", 1:K_all)
  
  # ---------- tidy up and end ----------
  
  # reorder qmatrices
  project <- align_qmatrix(project)
  
  # recalculate evidence over K
  project <- recalculate_evidence(project)
  
  # end timer
  if (!silent) {
    tdiff <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (tdiff<60) {
      message(sprintf("Total run-time: %s seconds", round(tdiff, 2)))
    } else {
      message(sprintf("Total run-time: %s minutes", round(tdiff/60, 2)))
    }
  }
  
  # warning if any rungs in any MCMCs did not converge
  if (!all_converged && !silent) {
    message("\n**WARNING** at least one MCMC run did not converge within specified burn-in\n")
  }
  
  # return invisibly
  invisible(project)
}

#------------------------------------------------
# extract GTI_logevidence from all K within a given parameter set
#' @noRd
get_GTI_logevidence_K <- function(proj, s) {
  
  # extract objects of interest
  x <- proj$output$single_set[[s]]$single_K
  if (length(x)==0) {
    return(proj)
  }
  
  # get log-evidence over all K
  GTI_logevidence_raw <- mapply(function(y) {
    GTI_logevidence <- y$summary$GTI_logevidence
    if (is.null(GTI_logevidence)) {
      return(rep(NA,2))
    } else {
      return(unlist(GTI_logevidence))
    }
  }, x)
  GTI_logevidence <- as.data.frame(t(GTI_logevidence_raw))
  names(GTI_logevidence) <- c("mean", "SE")
  rownames(GTI_logevidence) <- NULL
  GTI_logevidence <- cbind(K = 1:nrow(GTI_logevidence), GTI_logevidence)
  class(GTI_logevidence) <- "maverick_GTI_logevidence"
  
  # save result in project
  proj$output$single_set[[s]]$all_K$GTI_logevidence <- GTI_logevidence
  
  # return modified project
  return(proj)
}

#------------------------------------------------
# compute posterior over several log-evidence estimates
#' @noRd
get_GTI_posterior <- function(x) {
  
  # return NULL if all NA
  if (length(x$mean)==0 || all(is.na(x$mean))) {
    return(NULL)
  }
  
  # produce posterior estimates by simulation
  w <- which(!is.na(x$mean))
  GTI_posterior_raw <- GTI_posterior_K_sim_cpp(list(mean = x$mean[w],
                                                    SE = x$SE[w],
                                                    reps = 1e6))$ret
  # get posterior quantiles in dataframe
  GTI_posterior_quantiles <- t(mapply(quantile_95, GTI_posterior_raw))
  GTI_posterior <- data.frame(Q2.5 = rep(NA, nrow(GTI_posterior_quantiles)), Q50 = NA, Q97.5 = NA)
  GTI_posterior[w,] <- GTI_posterior_quantiles
  
  return(GTI_posterior)
}

#------------------------------------------------
# call get_GTI_posterior over values of K
#' @noRd
get_GTI_posterior_K <- function(proj, s) {
  
  # calculate posterior K
  GTI_posterior <- get_GTI_posterior(proj$output$single_set[[s]]$all_K$GTI_logevidence)
  if (is.null(GTI_posterior)) {
    return(proj)
  }
  GTI_posterior <- cbind(K = 1:nrow(GTI_posterior), GTI_posterior)
  class(GTI_posterior) <- "maverick_GTI_posterior"
  proj$output$single_set[[s]]$all_K$GTI_posterior <- GTI_posterior
  
  # return modified project
  return(proj)
}

#------------------------------------------------
# call get_GTI_posterior over models
#' @noRd
get_GTI_posterior_model <- function(proj) {
  
  # calculate posterior model
  GTI_posterior_model_raw <- get_GTI_posterior(proj$output$all_sets$GTI_logevidence_model)
  if (is.null(GTI_posterior_model_raw)) {
    return(proj)
  }
  proj$output$all_sets$GTI_posterior_model$Q2.5 <- GTI_posterior_model_raw$Q2.5
  proj$output$all_sets$GTI_posterior_model$Q50 <- GTI_posterior_model_raw$Q50
  proj$output$all_sets$GTI_posterior_model$Q97.5 <- GTI_posterior_model_raw$Q97.5
  
  # return modified project
  return(proj)
}

#------------------------------------------------
# integrate multiple log-evidence estimates by simulation
#' @noRd
integrate_GTI_logevidence <- function(x) {
  
  # return NULL if all NA
  if (length(x$mean)==0 || all(is.na(x$mean))) {
    return(NULL)
  }
  
  # produce integrated estimates by simulation
  w <- which(!is.na(x$mean))
  if (length(w)==1) {
    ret <- list(mean = x$mean[w], SE = x$SE[w])
  } else {
    ret <- GTI_integrated_K_sim_cpp(list(mean = x$mean[w], SE = x$SE[w], reps = 1e6))
  }
  
  return(ret)
}

#------------------------------------------------
# log-evidence estimates over K
#' @noRd
integrate_GTI_logevidence_K <- function(proj, s) {
  
  # integrate over K
  integrated_raw <- integrate_GTI_logevidence(proj$output$single_set[[s]]$all_K$GTI_logevidence)
  if (is.null(integrated_raw)) {
    return(proj)
  }
  proj$output$all_sets$GTI_logevidence_model$mean[s] <- integrated_raw$mean
  proj$output$all_sets$GTI_logevidence_model$SE[s] <- integrated_raw$SE
  
  # return modified project
  return(proj)
}

#------------------------------------------------
# align qmatrices over all K
#' @noRd
align_qmatrix <- function(proj) {
  
  # get active set
  s <- proj$active_set
  
  # extract objects of interest
  x <- proj$output$single_set[[s]]$single_K
  
  # find values with output
  null_output <- mapply(function(y) {is.null(y$summary$qmatrix_ind)}, x)
  w <- which(!null_output)
  
  # set template to first qmatrix
  template_qmatrix <- x[[w[1]]]$summary$qmatrix_ind
  n <- nrow(template_qmatrix)
  c <- ncol(template_qmatrix)
  
  # loop through output
  best_perm <- NULL
  for (i in w) {
    
    # expand template
    qmatrix_ind <- unclass(x[[i]]$summary$qmatrix_ind)
    template_qmatrix <- cbind(template_qmatrix, matrix(0, n, i-c))
    
    # calculate cost matrix
    cost_mat <- matrix(0,i,i)
    for (k1 in 1:i) {
        for (k2 in 1:i) {
          cost_mat[k1,k2] <- sum(qmatrix_ind[,k1] * (log(qmatrix_ind[,k1]+1e-100) - log(template_qmatrix[,k2]+1e-100)))
        }
    }
    
    # reorder qmatrix
    best_perm <- call_hungarian(cost_mat)$best_matching
    best_perm_order <- order(best_perm)
    qmatrix_ind <- qmatrix_ind[, best_perm_order, drop = FALSE]
    
    # qmatrix becomes template for next level up
    template_qmatrix <- qmatrix_ind
    
    # store result
    class(qmatrix_ind) <- "maverick_qmatrix_ind"
    proj$output$single_set[[s]]$single_K[[i]]$summary$qmatrix_ind <- qmatrix_ind
  }
  
  # return modified project
  return(proj)
}

#------------------------------------------------
#' @title Recalculate evidence and posterior estimates
#'
#' @description When a new value of K is added in to the analysis it affects all downstream evidence estimates that depend on this K - for example the overall model evidence integrated over K. This function therefore looks through all values of K in the active set and recalculates all downstream elements as needed.
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' 
#' @export

recalculate_evidence <- function(proj) {
  
  # check inputs
  assert_custom_class(proj, "mavproject")
  
  # get active set
  s <- proj$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # get log-evidence over all K and load into project
  proj <- get_GTI_logevidence_K(proj, s)
  
  # produce posterior estimates of K by simulation and load into project
  proj <- get_GTI_posterior_K(proj, s)
  
  # get log-evidence over all parameter sets
  proj <- integrate_GTI_logevidence_K(proj, s)
  
  # get posterior over all parameter sets
  proj <- get_GTI_posterior_model(proj)
  
  # return modified project
  return(proj)
}

#------------------------------------------------
#' @title Extract q-matrix for a given analysis
#'
#' @description Simple function for extracting the q-matrix output from a given
#'   parameter set (defaults to the active set) and value of K.
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K which value of K to extract
#' @param s which set to extract from. Defaults to the current active set
#'
#' @export

get_qmatrix <- function(proj, K, s = NULL) {
  
  # check inputs
  assert_custom_class(proj, "mavproject")
  
  # default to active set
  s <- define_default(s, proj$active_set)
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # return q-matrix
  return(proj$output$single_set[[s]]$single_K[[K]]$summary$qmatrix_ind)
}
