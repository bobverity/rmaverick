
#------------------------------------------------
#' @title Create a new rmaverick project
#'
#' @description Create a new rmaverick project.
#'
#' @details TODO
#'
#' @export
#' @examples
#' # TODO

mavproject <- function() {
  
  # create some empty data frames for storing results
  GTI_logevidence_model <- data.frame(set = numeric(),
                                      name = character(),
                                      mean = numeric(),
                                      SE = numeric())
  class(GTI_logevidence_model) <- "maverick_GTI_logevidence_model"
  
  GTI_posterior_model <- data.frame(set = numeric(),
                                    name = character(),
                                    Q2.5 = numeric(),
                                    Q50 = numeric(),
                                    Q97.5 = numeric())
  class(GTI_posterior_model) <- "maverick_GTI_posterior_model"
  
  # initialise project with default values
  ret <- list(data = NULL,
              data_processed = NULL,
              parameter_sets = NULL,
              active_set = 0,
              output = list(single_set = list(),
                            all_sets = list(GTI_logevidence_model = GTI_logevidence_model,
                                            GTI_posterior_model = GTI_posterior_model)
                            )
              )
  
  # create class and return
  class(ret) <- "mavproject"
  return(ret)
}

#------------------------------------------------
#' @title Print function for class mavproject
#'
#' @description Overload print function for class mavproject
#'
#' @details TODO
#'
#' @param x object of class \code{mavproject}
#' @param ... other parameters to pass to print function
#'
#' @export
#' @examples
#' # TODO

print.mavproject <- function(x, ...) {
  
  # print as list
  print(unclass(x))
  
  invisible(x)
}

#------------------------------------------------
#' @title Print summary for class mavproject
#'
#' @description Overload summary function for class mavproject
#'
#' @details TODO
#'
#' @param object object of class \code{mavproject}
#' @param ... other parameters to pass to summary function
#'
#' @export
#' @examples
#' # TODO

summary.mavproject <- function(object, ...) {
  
  # print data summary
  cat("DATA:\n")
  if (is.null(object$data)) {
    cat("   (none loaded)\n")
  } else {
    # extract data properties
    name <- object$data_processed$name
    ploidy <- object$data_processed$ploidy
    ploidy_min <- min(ploidy)
    ploidy_max <- max(ploidy)
    n <- length(ploidy)
    loci <- length(object$data_processed$Jl)
    pop <- object$data_processed$pop
    
    if (!is.null(name)) {
      cat(sprintf("   '%s'\n", name))
    }
    cat(sprintf("   individuals = %s\n", n))
    cat(sprintf("   loci = %s\n", loci))
    if (ploidy_min==ploidy_max) {
      cat(sprintf("   ploidy = %s\n", ploidy_min))
    } else {
      cat(sprintf("   ploidy range = %s to %s\n", ploidy_min, ploidy_max))
    }
    if (is.null(pop)) {
      cat("   pops = (none defined)\n")
    } else {
      cat(sprintf("   pops = %s sampled populations\n", length(unique(pop))))
    }
    
    n1 <- sum(is.na(object$data_processed$dat))
    n2 <- length(object$data_processed$dat)
    cat(sprintf("   missing data = %s of %s gene copies (%s%%)\n", n1, n2, round(n1/n2*100)))
  }
  cat("\n")
  
  # print parameter sets summary
  cat("PARAMETER SETS:\n")
  if (length(object$parameter_sets)==0) {
    cat("   (none defined)\n")
  } else {
    # print names of all sets
    s <- object$active_set
    for (i in 1:length(object$parameter_sets)) {
      
      # star next to active set
      if (i==s) {
        cat(" * ")
      } else {
        cat("   ")
      }
      
      # print name of set
      cat(sprintf("SET%s: %s\n", i, object$parameter_sets[[i]]$name))
    }
    cat("\n")
    
    # print details of active set
    name <- object$parameter_sets[[s]]$name
    admix_on <- object$parameter_sets[[s]]$admix_on
    estimate_alpha <- object$parameter_sets[[s]]$estimate_alpha
    alpha <- object$parameter_sets[[s]]$alpha
    
    cat(sprintf("ACTIVE SET: SET%s\n", s))
    if (admix_on) {
      cat(sprintf("   model = %s\n", "admixture"))
      cat(sprintf("   estimate alpha = %s\n", estimate_alpha))
      if (!estimate_alpha) {
        cat(sprintf("   alpha = %s\n", alpha))
      }
    } else {
      cat(sprintf("   model = %s\n", "no-admixture"))
    }
    
  }
  cat("\n")
  
  # print output summary
  cat("OUTPUT:\n")
  if (length(object$parameterSets)==0) {
    cat("   (none saved)\n")
  } else {
    for (i in 1:length(object$parameterSets)) {
    
    #percentComplete <- round(mean(x$output[[i]]$K_complete)*100, digits=2)
    #cat(sprintf("   SET%s: %s%% complete", i, percentComplete))
    
    }
  }

}

#------------------------------------------------
#' @title Determine if object is of class mavproject
#'
#' @description Determine if object is of class mavproject
#'
#' @details TODO
#'
#' @param x object of class \code{mavproject}
#'
#' @export
#' @examples
#' # TODO

is.mavproject <- function(x) {
  inherits(x, "mavproject")
}
