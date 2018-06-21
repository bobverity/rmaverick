
#------------------------------------------------
# define empty mavproject object
# (not exported)
#' @noRd
mavproject <- function() {
  
  # initialise project with default values
  ret <- list(data = NULL,
              data_processed = NULL,
              parameter_sets = NULL,
              active_set = 0,
              output = NULL
              )
  
  # create class and return
  class(ret) <- "mavproject"
  return(ret)
}

#------------------------------------------------
# overload print() function.
# (not exported)
#' @noRd
print.mavproject <- function(x, ...) {
  
  # print selected elements
  print(unclass(x)[c("data", "parameter_sets", "active_set", "output")])
  
  invisible(x)
}

#------------------------------------------------
# overload summary() function.
# (not exported)
#' @noRd
summary.mavproject <- function(x, ...) {
  
  # print data summary
  cat("DATA:\n")
  if (is.null(x$data)) {
    cat("   (none loaded)\n")
  } else {
    cat(sprintf("   '%s'\n", x$hidden$dataProcessed$name))
    cat(sprintf("   individuals = %s\n", x$hidden$dataProcessed$n))
    cat(sprintf("   loci = %s\n", x$hidden$dataProcessed$loci))
    
    ploidy <- x$hidden$dataProcessed$ploidy
    ploidy_min <- min(ploidy)
    ploidy_max <- max(ploidy)
    if (ploidy_min==ploidy_max) {
      cat(sprintf("   ploidy = %s\n", ploidy_min))
    } else {
      cat(sprintf("   ploidy range = %s to %s\n", ploidy_min, ploidy_max))
    }
    
    pop <- x$hidden$dataProcessed$pop
    if (is.null(pop)) {
      cat("   pops = (none defined)\n")
    } else {
      cat(sprintf("   pops = %s sampled populations\n", length(unique(pop))))
    }
    
    n1 <- x$hidden$dataProcessed$missingDataNum
    n2 <- x$hidden$dataProcessed$totalGeneCopies
    cat(sprintf("   missingData = %s of %s (%s%%)\n", n1, n2, round(n1/n2*100)))
  }
  cat("\n")
  
  # print parameter sets summary
  cat("PARAMETER SETS:\n")
  if (length(x$parameterSets)==0) {
    cat("   (none defined)\n")
  } else {
    for (i in 1:length(x$parameterSets)) {
    
    # star next to active set
    if (x$activeSet==i) {
      cat(" * ")
    } else {
      cat("   ")
    }
    
    # print details of set
    cat(sprintf("SET%s: '%s'\n", i, x$parameterSets[[i]]$setName))
    }
  }
  cat("\n")
  
  # print output summary
  cat("OUTPUT:\n")
  if (length(x$parameterSets)==0) {
    cat("   (none saved)\n")
  } else {
    for (i in 1:length(x$parameterSets)) {
    
    #percentComplete <- round(mean(x$output[[i]]$K_complete)*100, digits=2)
    #cat(sprintf("   SET%s: %s%% complete", i, percentComplete))
    
    }
  }

}

#------------------------------------------------
# function for determining if object is of class mavProject
# (not exported)
#' @noRd
is.mavproject <- function(x) {
  inherits(x, "mavproject")
}
