
#------------------------------------------------
# default MavericK colours
# (not exported)
#' @noRd
default_colours <- function(K) {
  
  # generate palette and colours
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  my_palette <- colorRampPalette(raw_cols)
  bar_col <- my_palette(max(K,6))
  
  # if fewer than 6 colours then choose manually
  if (K==5) {
    bar_col = bar_col[c(1,2,3,5,6)]
  } else if (K==4) {
    bar_col = bar_col[c(1,2,3,5)]
  } else if (K==3) {
    bar_col = bar_col[c(1,3,5)]
  } else if (K==2) {
    bar_col = bar_col[c(1,5)]
  } else if (K==1) {
    bar_col = bar_col[1]
  }
  
  return(bar_col)
}

#------------------------------------------------
# whisker plot of quantiles
# (not exported)
#' @noRd
plot_quantiles <- function(q_min, q_mid, q_max, q_x = 1:length(q_min), width = 0.2, connect_points = FALSE, ...) {
  
  # get input arguments
  args <- list(...)
  arg_names <- names(args)
  
  # check inputs
  n <- length(q_min)
  stopifnot(length(q_mid)==n && length(q_max)==n)
  stopifnot(all(q_min<=q_mid, na.rm=TRUE) && all(q_mid<=q_max, na.rm=TRUE))
  
  # get basic data properties
  min_val <- min(q_min, na.rm=TRUE)
  max_val <- max(q_max, na.rm=TRUE)
  
  # stick properties
  stick_left <- 1:n - 0.5
  stick_right <- 1:n + 0.5
  
  # set defaults on undefined arguments
  if (! "xlim" %in% arg_names) {
    args$xlim <- range(q_x)
  }
  if (! "ylim" %in% arg_names) {
    y_mid <- 0.5*(max_val + min_val)
    y_diff <- 0.5*(max_val - min_val)
    args$ylim <- c(y_mid - 1.2*y_diff, y_mid + 1.2*y_diff)
  }
  if (! "lwd" %in% arg_names) {
    args$lwd <- 1
  }
  if (! "col" %in% arg_names) {
    args$col <- 1
  }
  if (! "yaxs" %in% arg_names) {
    args$yaxs <- "i"
  }
  
  # fixed arguments, or arguments that have special meaning
  args$type <- "n"
  axes <- TRUE
  if ("axes" %in% arg_names) {
    axes <- args$axes
  }
  args$axes <- FALSE
  
  # plot with finalised list of parameters
  do.call(plot, c(list(x=0), args))    # create empty plot
  segments(x0=q_x-width, y0=q_mid, x1=q_x+width, y1=q_mid, lwd=args$lwd, col=args$col)  # add horizontal lines
  segments(x0=q_x, y0=q_min, x1=q_x, y1=q_max, lwd=args$lwd, col=args$col)  # add vertical lines
  if (axes) {
    axis(1)
    axis(2)
    box()
  }
  if (connect_points) {
    lines(q_x, q_mid, lwd=args$lwd[1], col=args$col[1])
  }
  
}

#------------------------------------------------
#' @title Default plot for class maverick_loglike_quantiles
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param x TODO
#' @param y TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot.maverick_loglike_quantiles <- function(x, y, ...) {
  
  # get input arguments
  args <- list(...)
  arg_names <- names(args)
  
  # unclass x
  x <- unclass(x)
  n <- nrow(x)
  
  # get basic data properties
  q_min <- as.vector(x[,1])
  q_mid <- as.vector(x[,2])
  q_max <- as.vector(x[,3])
  
  # set defaults on undefined arguments
  if (! "xlab" %in% arg_names) {
    args$xlab <- "rung"
  }
  if (! "ylab" %in% arg_names) {
    args$ylab <- "log-likelihood"
  }
  
  # fixed arguments, or arguments that have special meaning
  if ("q_x" %in% arg_names) {
    q_min <- c(NA, q_min)
    q_mid <- c(NA, q_mid)
    q_max <- c(NA, q_max)
  }
  
  # plot with finalised list of parameters
  do.call(plot_quantiles, c(list(q_min=q_min, q_mid=q_mid, q_max=q_max), args))
}

#------------------------------------------------
#' @title Plot loglike quantiles of current active set
#'
#' @description Plot loglike quantiles of current active set
#'
#' @details TODO
#'
#' @param proj TODO
#' @param K TODO
#' @param axis_type TODO
#' @param connect_points TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot_loglike_quantiles <- function(proj, K = NULL, axis_type = 3, connect_points = FALSE, ...) {
  
  # get input arguments
  args <- list(...)
  arg_names <- names(args)
  
  # get active set and check non-zero
  s <- proj$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(is.null, proj$output[[s]])
  if (all(null_output)) {
    stop("no output for active parameter set")
  }
  K <- define_default(K, which(!null_output)[1])
  
  # check output exists for this K
  if (is.null(proj$output[[s]][K])) {
    stop(sprintf("no output for K = %s of active set", K))
  }
  
  # x axis options
  if (axis_type %in% c(2,3)) {
    rungs <- proj$output[[s]][[K]]$function_call$args$rungs
    if (! "width" %in% arg_names) {
      args$width <- 0.02
    }
    if (! "xlim" %in% arg_names) {
      args$xlim <- c(0,1)
    }
    if (axis_type==2) {
      args$q_x <- (0:rungs)/rungs
      if (! "xlab" %in% arg_names) {
        args$xlab <- parse(text="beta")
      }
    }
    if (axis_type==3) {
      GTI_pow <- proj$output[[s]][[K]]$function_call$args$GTI_pow
      args$q_x <- ((0:rungs)/rungs)^GTI_pow
      if (! "xlab" %in% arg_names) {
        args$xlab <- parse(text="beta^gamma")
      }
    }
  }
  
  # produce quantile plot with finalised list of parameters
  do.call(plot, c(list(proj$output[[s]][[K]]$summary$loglike_quantiles), args))
}

#------------------------------------------------
#' @title Default plot for class maverick_qmatrix_ind
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param x TODO
#' @param y TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot.maverick_qmatrix_ind <- function(x, y = NULL, ...) {
  
  # get input arguments
  args <- list(...)
  arg_names <- names(args)
  
  # unclass x
  x <- t(unclass(x))
  
  # get basic data properties
  n <- ncol(x)
  K <- nrow(x)
  
  # set defaults on undefined arguments
  if (! "ylim" %in% arg_names) {
    args$ylim <- c(0, 1)
  }
  if (! "xlab" %in% arg_names) {
    args$xlab <- "sample"
  }
  if (! "ylab" %in% arg_names) {
    args$ylab <- ""
  }
  if (! "xaxs" %in% arg_names) {
    args$xaxs <- "i"
  }
  if (! "yaxs" %in% arg_names) {
    args$yaxs <- "i"
  }
  if (! "space" %in% arg_names) {
    args$space <- 0
  }
  if (! "names" %in% arg_names) {
    args$names <- rep(NA,n)
  }
  if (! "col" %in% arg_names) {
    args$col <- default_colours(K)
  }
  if (! "border" %in% arg_names) {
    args$border <- NA
  }
  
  # plot with finalised list of parameters
  do.call(barplot, c(list(height=x), args))
  box()
  
  # if y data used, add points above barplot
  # NOTE - RELEGATED UNTIL MOVE TO GGPLOT
  # TODO - replace this code
  if (!is.null(y)) {
    
    # TODO - checks on y data
    
    # variation of above arguments
    args2 <- args
    args2$axes <- FALSE
    args2$ylim <- c(-50,-1)
    args2$col <- args$col[y]
    
    # allow plot outside plotting region
    par_store <- par(xpd = NA, new = TRUE)
    on.exit(par(par_store))
    
    # add secondary barplot
    do.call(barplot, c(list(height=rep(1,length(y))), args2))
  }
  
}

#------------------------------------------------
#' @title Plot individual-level Q-matrix of current active set
#'
#' @description Plot individual-level Q-matrix of current active set
#'
#' @details TODO
#'
#' @param proj TODO
#' @param K TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot_qmatrix_ind <- function(proj, K = NULL, ...) {
  
  # check inputs
  assert_mavproject(proj)
  if (!is.null(K)) {
    assert_scalar_pos_int(K)
  }
  
  # get active set and check non-zero
  s <- proj$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(is.null, proj$output[[s]])
  if (all(null_output)) {
    stop("no output for active parameter set")
  }
  K <- define_default(K, which(!null_output)[1])
  
  # check output exists for this K
  if (is.null(proj$output[[s]][K])) {
    stop(sprintf("no output for K = %s of active set", K))
  }
  
  # produce Q-matrix plot
  plot(proj$output[[s]][[K]]$summary$qmatrix_ind, ...)
}

#------------------------------------------------
#' @title Default plot for class maverick_GTI_path
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param x TODO
#' @param y TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot.maverick_GTI_path <- function(x, y, ...) {
  
  # get input arguments
  args <- list(...)
  arg_names <- names(args)
  
  # unclass x
  n <- nrow(x)
  x <- unclass(x)
  
  # get basic data properties
  q_mid <- x$mean
  q_min <- x$mean - 1.96*x$SE
  q_max <- x$mean + 1.96*x$SE
  
  # set defaults on undefined arguments
  if (! "xlab" %in% arg_names) {
    args$xlab <- "rung"
  }
  if (! "ylab" %in% arg_names) {
    args$ylab <- "GTI path"
  }
  
  # fixed arguments, or arguments that have special meaning
  if ("q_x" %in% arg_names) {
    q_min <- c(0, q_min)
    q_mid <- c(0, q_mid)
    q_max <- c(0, q_max)
  }
  
  # plot with finalised list of parameters
  do.call(plot_quantiles, c(list(q_min=q_min, q_mid=q_mid, q_max=q_max), args))
}

#------------------------------------------------
#' @title Plot GTI path of current active set
#'
#' @description Plot GTI path of current active set
#'
#' @details TODO
#'
#' @param proj TODO
#' @param K TODO
#' @param axis_type TODO
#' @param connect_points TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot_GTI_path <- function(proj, K = NULL, axis_type = 2, connect_points = TRUE, ...) {
  
  # get input arguments
  args <- list(...)
  arg_names <- names(args)
  
  # get active set and check non-zero
  s <- proj$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(is.null, proj$output[[s]])
  if (all(null_output)) {
    stop("no output for active parameter set")
  }
  K <- define_default(K, which(!null_output)[1])
  
  # check output exists for this K
  if (is.null(proj$output[[s]][K])) {
    stop(sprintf("no output for K = %s of active set", K))
  }
  
  # x axis options
  if (axis_type %in% c(2,3)) {
    rungs <- proj$output[[s]][[K]]$function_call$args$rungs
    if (! "width" %in% arg_names) {
      args$width <- 0.02
    }
    if (! "xlim" %in% arg_names) {
      args$xlim <- c(0,1)
    }
    if (axis_type==2) {
      args$q_x <- (0:rungs)/rungs
      if (! "xlab" %in% arg_names) {
        args$xlab <- parse(text="beta")
      }
    }
    if (axis_type==3) {
      GTI_pow <- proj$output[[s]][[K]]$function_call$args$GTI_pow
      args$q_x <- ((0:rungs)/rungs)^GTI_pow
      if (! "xlab" %in% arg_names) {
        args$xlab <- parse(text="beta^gamma")
      }
    }
  }
  
  # default args
  args$connect_points <- connect_points
  
  # produce quantile plot with finalised list of parameters
  do.call(plot, c(list(proj$output[[s]][[K]]$summary$GTI_path), args))
}

#------------------------------------------------
#' @title Plot model log-evidence estimates over K
#'
#' @description Plot model log-evidence estimates over K
#'
#' @details TODO
#'
#' @param proj TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot_logevidence <- function(proj, ...) {
  
  # get input arguments
  args <- list(...)
  arg_names <- names(args)
  
  # get active set and check non-zero
  s <- proj$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set K to all values with output
  null_output <- mapply(is.null, proj$output[[s]])
  if (all(null_output)) {
    stop("no output for active parameter set")
  }
  K <- which(!null_output)
  
  # get logevidence over all K
  GTI_logevidence <- mapply(function(x) {
    GTI_logevidence <- x$summary$GTI_logevidence
    if (is.null(GTI_logevidence)) {
      return(rep(NA,2))
    } else {
      return(GTI_logevidence)
    }
  }, proj$output[[s]])
  
  # set defaults
  if (! "xlab" %in% arg_names) {
    args$xlab <- "K"
  }
  if (! "ylab" %in% arg_names) {
    args$ylab <- "log-evidence"
  }
  
  # get confidence intervals
  q_mid <- unlist(GTI_logevidence[1,])
  q_min <- q_mid - 1.96*unlist(GTI_logevidence[2,])
  q_max <- q_mid + 1.96*unlist(GTI_logevidence[2,])
  
  # plot with finalised list of parameters
  do.call(plot_quantiles, c(list(q_min=q_min, q_mid=q_mid, q_max=q_max, q_x=K), args))
}

#------------------------------------------------
#' @title Plot model evidence estimates over K
#'
#' @description Plot model evidence estimates over K
#'
#' @details TODO
#'
#' @param proj TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot_evidence <- function(proj, ...) {
  
  # get input arguments
  args <- list(...)
  arg_names <- names(args)
  
  # get active set and check non-zero
  s <- proj$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # get evidence over all K
  GTI_evidence <- mapply(function(x) {
    GTI_evidence <- x$summary$GTI_evidence
    if (is.null(GTI_evidence)) {
      return(rep(NA,3))
    } else {
      return(GTI_evidence)
    }
  }, proj$output[[s]])
  
  # set defaults
  if (! "xlab" %in% arg_names) {
    args$xlab <- "K"
  }
  if (! "ylab" %in% arg_names) {
    args$ylab <- "evidence"
  }
  if (! "ylim" %in% arg_names) {
    args$ylim <- c(0,1)
  }
  if (! "names.arg" %in% arg_names) {
    args$names.arg <- 1:ncol(GTI_evidence)
  }
  if (! "lwd" %in% arg_names) {
    args$lwd <- 1
  }
  if (! "col" %in% arg_names) {
    args$col <- "#4575B4"
  }
  args$space <- 0
  
  # get confidence intervals
  q_min <- unlist(GTI_evidence[1,])
  q_mid <- unlist(GTI_evidence[2,])
  q_max <- unlist(GTI_evidence[3,])
  
  # produce plot
  do.call(barplot, c(list(height = q_mid), args))
  segments(x0=1:length(q_mid)-0.5, y0=q_min, x1=1:length(q_mid)-0.5, y1=q_max, lwd=args$lwd)  # add vertical lines
}
