
#------------------------------------------------
# default rmaverick colours
#' @importFrom grDevices colorRampPalette
#' @noRd
default_colours <- function(K) {
  
  # generate palette and colours
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  my_palette <- grDevices::colorRampPalette(raw_cols)
  
  # simple case if small K
  if (K <= 2) {
    return(my_palette(K))
  }
  
  # some logic to choose a palette size and sequence of colours that is
  # consistent across different values of K
  ncol <- 3
  while(ncol<K) {
    ncol <- ncol+(ncol-1)
  }
  dist_mat <- matrix(1:ncol, ncol, ncol)
  dist_mat <- abs(t(dist_mat)-dist_mat)
  x <- rep(FALSE, ncol)
  
  col_index <- 1
  for (i in 2:K) {
    x[col_index] <- TRUE
    s <- apply(dist_mat[which(x),,drop=FALSE], 2, min)
    next_index <- which.max(s)
    col_index <- c(col_index, next_index)
  }
  col_index
  ret <- my_palette(ncol)[col_index]
  
  return(ret)
}

#------------------------------------------------
# ggplot theme with minimal objects
#' @import ggplot2
#' @noRd
theme_empty <- function() {
  ggplot2::theme(axis.line = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 panel.border = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank())
}

#------------------------------------------------
# Default plot for class maverick_loglike_quantiles
#' @import ggplot2
#' @noRd
plot.maverick_loglike_quantiles <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  x_vec <- y
  
  # produce plot
  plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw()
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~x_vec, y = ~Q2.5, xend = ~x_vec, yend = ~Q97.5))
  plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes_(x = ~x_vec, y = ~Q50))
  plot1 <- plot1 + ggplot2::xlab("rung") + ggplot2::ylab("log-likelihood")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot loglike quantiles of current active set
#'   
#' @description Plot loglike quantiles of current active set
#'   
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K which value of K to produce the plot for
#' @param axis_type how to format the x-axis. 1 = integer rungs, 2 = values of
#'   beta, 3 = values of beta raised to the GTI power
#' @param connect_points whether to connect points in the middle of quantiles
#' @param connect_whiskers whether to connect points at the ends of the whiskers
#'
#' @import ggplot2
#' @export

plot_loglike_quantiles <- function(proj, K = NULL, axis_type = 1, connect_points = FALSE, connect_whiskers = FALSE) {
  
  # check inputs
  assert_class(proj, "mavproject")
  if (!is.null(K)) {
    assert_single_pos_int(K)
  }
  assert_in(axis_type, 1:3)
  assert_single_logical(connect_points)
  assert_single_logical(connect_whiskers)
  
  # get active set and check non-zero
  s <- proj$active_set
  assert_neq(s, 0, message = "no active parameter set")
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$loglike_quantiles)}, proj$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no loglike_quantiles output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  assert_leq(K, length(proj$output$single_set[[s]]$single_K), message = "K out of bounds")
  loglike_quantiles <- proj$output$single_set[[s]]$single_K[[K]]$summary$loglike_quantiles
  if (is.null(loglike_quantiles)) {
    stop(sprintf("no loglike_quantiles output for K = %s of active set", K))
  }
  
  # produce plot with different axis options
  rungs <- nrow(loglike_quantiles)
  if (axis_type == 1) {
    x_vec <- 1:rungs
    plot1 <- plot(loglike_quantiles, as.factor(x_vec))
    
  } else if (axis_type == 2) {
    x_vec <- (1:rungs) / rungs
    plot1 <- plot(loglike_quantiles, x_vec)
    plot1 <- plot1 + ggplot2::xlab(parse(text = "beta"))
    plot1 <- plot1 + ggplot2::coord_cartesian(xlim = c(0,1))
    
  } else {
    GTI_pow <- proj$output$single_set[[s]]$single_K[[K]]$function_call$args$GTI_pow
    x_vec <- ((1:rungs) / rungs)^GTI_pow
    plot1 <- plot(loglike_quantiles, x_vec)
    plot1 <- plot1 + ggplot2::xlab(parse(text = "beta^gamma"))
    plot1 <- plot1 + ggplot2::coord_cartesian(xlim = c(0,1))
  }
  
  # optionally add central line
  if (connect_points) {
    df <- as.data.frame(unclass(loglike_quantiles))
    plot1 <- plot1 + ggplot2::geom_line(ggplot2::aes(x = x_vec, y = df$Q50))
  }
  
  # optionally connect whiskers
  if (connect_whiskers) {
    df <- as.data.frame(unclass(loglike_quantiles))
    plot1 <- plot1 + ggplot2::geom_line(ggplot2::aes(x = x_vec, y = df$Q2.5), linetype = "dotted") + ggplot2::geom_line(ggplot2::aes(x = x_vec, y = df$Q97.5), linetype = "dotted")
  }
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class maverick_qmatrix_ind
#' @import ggplot2
#' @noRd
plot.maverick_qmatrix_ind <- function(x, y, ...) {
  
  # get data into ggplot format
  m <- unclass(x)
  n <- nrow(m)
  K <- ncol(m)
  df <- data.frame(ind = rep(1:n,each=K), k = as.factor(rep(1:K,times=n)), val = as.vector(t(m)))
  
  # produce basic plot
  plot1 <- ggplot2::ggplot(df) + theme_empty()
  plot1 <- plot1 + ggplot2::geom_bar(ggplot2::aes_(x = ~ind, y = ~val, fill = ~k), width = 1, stat = "identity")
  plot1 <- plot1 + ggplot2::scale_x_continuous(expand = c(0,0)) + ggplot2::scale_y_continuous(expand = c(0,0))
  plot1 <- plot1 + ggplot2::xlab("sample") + ggplot2::ylab("probability")
  
  # add legends
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = default_colours(K), name = "group")
  plot1 <- plot1 + ggplot2::scale_colour_manual(values = "white")
  plot1 <- plot1 + ggplot2::guides(colour = "none")
  
  # add border
  plot1 <- plot1 + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 2, fill = NA))
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot Q-matrix of current active set
#'
#' @description Plot Q-matrix of current active set
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K which value of K to produce the plot for
#' @param divide_ind_on whether to add dividing lines between bars
#'
#' @import ggplot2
#' @export

plot_qmatrix <- function(proj, K = NULL, divide_ind_on = FALSE) {
  
  # check inputs
  assert_class(proj, "mavproject")
  if (!is.null(K)) {
    assert_pos_int(K)
  }
  assert_single_logical(divide_ind_on)
  
  # get active set and check non-zero
  s <- proj$active_set
  assert_neq(s, 0, message = "no active parameter set")
  
  # set default K to all values with output
  null_output <- mapply(function(x) {is.null(x$summary$qmatrix_ind)}, proj$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no output for active parameter set")
  }
  K <- define_default(K, which(!null_output))
  
  # check output exists for chosen K
  assert_leq(K, length(proj$output$single_set[[s]]$single_K), message = "K out of bounds")
  qmatrix_ind_list <- list()
  for (i in 1:length(K)) {
    qmatrix_ind_list[[i]] <- proj$output$single_set[[s]]$single_K[[K[i]]]$summary$qmatrix_ind
    if (is.null(qmatrix_ind_list[[i]])) {
      stop(sprintf("no qmatrix_ind output for K = %s of active set", K[i]))
    }
  }
  n <- nrow(qmatrix_ind_list[[1]])
  
  # get data into ggplot format
  df <- NULL
  for (i in 1:length(K)) {
    m <- unclass(qmatrix_ind_list[[i]])
    df <- rbind(df, data.frame(K = as.numeric(K[i]), ind = rep(1:n,each=K[i]), k = as.factor(rep(1:K[i],times=n)), val = as.vector(t(m))))
  }
  
  # produce basic plot
  plot1 <- ggplot2::ggplot(df) + theme_empty()
  plot1 <- plot1 + ggplot2::geom_bar(ggplot2::aes_(x = ~ind, y = ~val, fill = ~k), width = 1, stat = "identity")
  plot1 <- plot1 + ggplot2::scale_x_continuous(expand = c(0,0)) + ggplot2::scale_y_continuous(expand = c(0,0))
  
  # arrange in rows
  if (length(K)==1) {
    plot1 <- plot1 + ggplot2::facet_wrap(~K, ncol = 1)
    plot1 <- plot1 + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_blank())
    plot1 <- plot1 + ggplot2::xlab("sample") + ggplot2::ylab("probability")
  } else {
    plot1 <- plot1 + ggplot2::facet_wrap(~K, ncol = 1, strip.position = "left")
    plot1 <- plot1 + ggplot2::theme(strip.background = ggplot2::element_blank())
    plot1 <- plot1 + ggplot2::xlab("sample") + ggplot2::ylab("K")
  }
  
  # add legends
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = default_colours(max(K)), name = "group")
  plot1 <- plot1 + ggplot2::scale_colour_manual(values = "white")
  plot1 <- plot1 + ggplot2::guides(colour = "none")
  
  # add border
  plot1 <- plot1 + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 2, fill = NA))
  
  # optionally add dividing lines
  if (divide_ind_on) {
    plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~x, y = ~y, xend = ~x, yend = ~y+1, col = "white"), size = 0.3, data = data.frame(x = 1:n-0.5, y = rep(0,n)))
  }
  
  return(plot1)
}

#------------------------------------------------
# Default plot for class maverick_GTI_path
#' @import ggplot2
#' @noRd
plot.maverick_GTI_path <- function(x, y, ...) {
  
  # check inputs
  assert_in(y, 1:2)
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  df <- rbind(data.frame(mean = 0, SE = 0), df)
  
  # get quantiles
  df$q_min <- df$mean - 1.96*df$SE
  df$q_mid <- df$mean
  df$q_max <- df$mean + 1.96*df$SE
  
  # produce plot
  plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw()
  if (y == 1) {
    q_x <- 1:(n + 1)
    width <- 0.1
    plot1 <- plot1 + ggplot2::geom_line(ggplot2::aes_(x = ~as.factor(0:n), y = ~q_mid, group = 1))
    plot1 <- plot1 + ggplot2::xlab("rung")
  } else {
    q_x <- seq(0, 1, l = n + 1)
    width <- 0.01
    plot1 <- plot1 + ggplot2::geom_line(ggplot2::aes_(x = ~q_x, y = ~q_mid))
    plot1 <- plot1 + ggplot2::xlab(parse(text = "beta"))
  }
  
  # continue building plot
  plot1 <- plot1 + ggplot2::geom_area(ggplot2::aes_(x = ~q_x, y = ~q_mid, fill = "col1", colour = "col1", alpha = 0.5))
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~q_x, y = ~q_min, xend = ~q_x, yend = ~q_max))
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~q_x-width, y = ~q_min, xend = ~q_x+width, yend = ~q_min))
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~q_x-width, y = ~q_max, xend = ~q_x+width, yend = ~q_max))
  
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = "#4575B4")
  plot1 <- plot1 + ggplot2::scale_colour_manual(values = "black")
  plot1 <- plot1 + ggplot2::guides(fill = "none", colour = "none", alpha = "none")
  plot1 <- plot1 + ggplot2::ylab("weighted log-likelihood")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot GTI path of current active set
#'
#' @description Plot GTI path of current active set
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K which value of K to produce the plot for
#' @param axis_type how to format the x-axis. 1 = integer rungs, 2 = values of
#'   beta
#'
#' @export

plot_GTI_path <- function(proj, K = NULL, axis_type = 1) {
  
  # check inputs
  assert_class(proj, "mavproject")
  if (!is.null(K)) {
    assert_single_pos_int(K)
  }
  assert_in(axis_type, 1:2)
  
  # get active set and check non-zero
  s <- proj$active_set
  assert_neq(s, 0, message = "no active parameter set")
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$GTI_path)}, proj$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  assert_leq(K, length(proj$output$single_set[[s]]$single_K), message = "K out of bounds")
  GTI_path <- proj$output$single_set[[s]]$single_K[[K]]$summary$GTI_path
  if (is.null(GTI_path)) {
    stop(sprintf("no GTI_path output for K = %s of active set", K))
  }
  
  # produce plot with different axis options
  plot1 <- plot(GTI_path, axis_type)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class maverick_GTI_logevidence
#' @import ggplot2
#' @noRd
plot.maverick_GTI_logevidence <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  
  # get quantiles
  df$q_min <- df$mean - 1.96*df$SE
  df$q_mid <- df$mean
  df$q_max <- df$mean + 1.96*df$SE
  q_x <- 1:n
  
  # produce plot
  plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw()
  width <- 0.1
  plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes_(x = ~as.factor(K), y = ~q_mid), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~q_x, y = ~q_min, xend = ~q_x, yend = ~q_max), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~q_x - width, y = ~q_min, xend = ~q_x + width, yend = ~q_min), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~q_x - width, y = ~q_max, xend = ~q_x + width, yend = ~q_max), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::xlab("K") + ggplot2::ylab("log-evidence")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot log-evidence estimates over K
#'
#' @description Plot log-evidence estimates over K
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#'
#' @export

plot_logevidence_K <- function(proj) {
  
  # check inputs
  assert_class(proj, "mavproject")
  
  # get active set and check non-zero
  s <- proj$active_set
  assert_neq(s, 0, message = "no active parameter set")
  
  # check output exists for chosen K
  GTI_logevidence <- proj$output$single_set[[s]]$all_K$GTI_logevidence
  if (is.null(GTI_logevidence)) {
    stop("no GTI_logevidence output for active set")
  }
  
  # produce plot
  plot1 <- plot(GTI_logevidence)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class maverick_GTI_posterior
#' @import ggplot2
#' @noRd
plot.maverick_GTI_posterior <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  width <- 0.1
  
  # produce plot
  plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw()
  plot1 <- plot1 + ggplot2::geom_bar(ggplot2::aes_(x = ~K, y = ~Q50, fill = "blue"), stat = "identity", na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~K, y = ~Q2.5, xend = ~K, yend = ~Q97.5), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~K-width, y = ~Q2.5, xend = ~K+width, yend = ~Q2.5), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~K-width, y = ~Q97.5, xend = ~K+width, yend = ~Q97.5), na.rm = TRUE)
  
  # add legends
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = "#4575B4")
  plot1 <- plot1 + ggplot2::guides(fill = "none")
  
  # modify scales etc.
  plot1 <- plot1 + ggplot2::coord_cartesian(ylim = c(-0.05,1.05))
  plot1 <- plot1 + ggplot2::scale_y_continuous(expand = c(0,0))
  plot1 <- plot1 + ggplot2::xlab("K") + ggplot2::ylab("probability")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot posterior K
#'
#' @description Plot posterior K
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#'
#' @export

plot_posterior_K <- function(proj) {
  
  # check inputs
  assert_class(proj, "mavproject")
  
  # get active set and check non-zero
  s <- proj$active_set
  assert_neq(s, 0, message = "no active parameter set")
  
  # check output exists for chosen K
  GTI_posterior <- proj$output$single_set[[s]]$all_K$GTI_posterior
  if (is.null(GTI_posterior)) {
    stop("no GTI_posterior output for active set")
  }
  
  # produce plot
  plot1 <- plot(GTI_posterior)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class maverick_GTI_logevidence_model
#' @import ggplot2
#' @noRd
plot.maverick_GTI_logevidence_model <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  
  # get quantiles
  df$q_min <- df$mean - 1.96*df$SE
  df$q_mid <- df$mean
  df$q_max <- df$mean + 1.96*df$SE
  q_x <- 1:n
  
  # produce plot
  plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw()
  width <- 0.1
  plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes_(x = ~as.factor(set), y = ~q_mid), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~q_x, y = ~q_min, xend = ~q_x, yend = ~q_max), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~q_x-width, y = ~q_min, xend = ~q_x+width, yend = ~q_min), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~q_x-width, y = ~q_max, xend = ~q_x+width, yend = ~q_max), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::scale_x_discrete(labels = df$name, breaks = 1:n)
  plot1 <- plot1 + ggplot2::xlab("model") + ggplot2::ylab("log-evidence")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot log-evidence estimates over parameter sets
#'
#' @description Plot log-evidence estimates over parameter sets
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#'
#' @export

plot_logevidence_model <- function(proj) {
  
  # check inputs
  assert_class(proj, "mavproject")
  
  # check output exists
  GTI_logevidence_model <- proj$output$all_sets$GTI_logevidence_model
  if (is.null(GTI_logevidence_model)) {
    stop("no GTI_logevidence_model output")
  }
  
  # produce plot
  plot1 <- plot(GTI_logevidence_model)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class maverick_GTI_posterior_model
#' @import ggplot2
#' @noRd
plot.maverick_GTI_posterior_model <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  width <- 0.1
  
  # produce plot
  plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw()
  plot1 <- plot1 + ggplot2::geom_bar(ggplot2::aes_(x = ~as.factor(set), y = ~Q50, fill = "blue"), stat = "identity", na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~set, y = ~Q2.5, xend = ~set, yend = ~Q97.5), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~set-width, y = ~Q2.5, xend = ~set+width, yend = ~Q2.5), na.rm = TRUE)
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~set-width, y = ~Q97.5, xend = ~set+width, yend = ~Q97.5), na.rm = TRUE)
  
  # add legends
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = "#4575B4")
  plot1 <- plot1 + ggplot2::guides(fill = "none")
  
  # modify scales etc.
  plot1 <- plot1 + ggplot2::scale_x_discrete(labels = df$name, breaks = 1:n)
  plot1 <- plot1 + ggplot2::coord_cartesian(ylim = c(-0.05,1.05))
  plot1 <- plot1 + ggplot2::scale_y_continuous(expand = c(0,0))
  plot1 <- plot1 + ggplot2::xlab("model") + ggplot2::ylab("probability")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot posterior model
#'
#' @description Plot posterior model
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#'
#' @export

plot_posterior_model <- function(proj) {
  
  # check inputs
  assert_class(proj, "mavproject")
  
  # check output exists
  GTI_posterior_model <- proj$output$all_sets$GTI_posterior_model
  if (is.null(GTI_posterior_model)) {
    stop("no GTI_posterior_model output")
  }
  
  # produce plot
  plot1 <- plot(GTI_posterior_model)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce MCMC trace plot
#'   
#' @description Produce MCMC trace plot of alpha or log-likelihood
#'   
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K which value of K to produce the plot for
#' @param rung which value of K to produce the plot for. Defaults to the cold rung
#' @param param whether to produce trace plot of admixture parameter
#'   (\code{"alpha"}) or the log-likelihood (\code{"loglike"})
#' @param col colour of the trace
#'   
#' @import ggplot2
#' @export

plot_trace <- function(proj, K = NULL, rung = NULL, param = "alpha", col = "black") {
  
  # check inputs
  assert_class(proj, "mavproject")
  if (!is.null(K)) {
    assert_single_pos_int(K)
  }
  if (!is.null(rung)) {
    assert_single_pos_int(rung)
  }
  assert_in(param, c("alpha", "loglike"))
  
  # get active set and check non-zero
  s <- proj$active_set
  assert_neq(s, 0, message = "no active parameter set")
  
  # split analysis between alpha and loglikelihood
  if (param == "alpha") {
    
    # set default K to first value with output
    null_output <- mapply(function(x) {is.null(x$raw$alpha)}, proj$output$single_set[[s]]$single_K)
    if (all(null_output)) {
      stop("no alpha output for active parameter set")
    }
    if (is.null(K)) {
      K <- which(!null_output)[1]
      message(sprintf("using K = %s by default", K))
    }
    
    # check output exists for chosen K
    assert_leq(K, length(proj$output$single_set[[s]]$single_K), message = "K out of bounds")
    alpha <- as.vector(proj$output$single_set[[s]]$single_K[[K]]$raw$alpha)
    if (is.null(alpha)) {
      stop(sprintf("no alpha output for K = %s of active set", K))
    }
    
    # get into ggplot format
    df <- data.frame(x = 1:length(alpha), y = alpha)
    
    # produce plot
    plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw() + ggplot2::ylab("alpha")
    
  } else {
    
    # set default K to first value with output
    null_output <- mapply(function(x) {is.null(x$raw$loglike_sampling)}, proj$output$single_set[[s]]$single_K)
    if (all(null_output)) {
      stop("no loglike_sampling output for active parameter set")
    }
    if (is.null(K)) {
      K <- which(!null_output)[1]
      message(sprintf("using K = %s by default", K))
    }
    
    # check output exists for chosen K
    assert_leq(K, length(proj$output$single_set[[s]]$single_K), message = "K out of bounds")
    loglike_sampling <- proj$output$single_set[[s]]$single_K[[K]]$raw$loglike_sampling
    if (is.null(loglike_sampling)) {
      stop(sprintf("no loglike_sampling output for K = %s of active set", K))
    }
    
    # use cold rung by default
    rung <- define_default(rung, ncol(loglike_sampling))
    assert_leq(rung, ncol(loglike_sampling))
    loglike <- as.vector(loglike_sampling[,rung])
    
    # get into ggplot format
    df <- data.frame(x = 1:length(loglike), y = loglike)
    
    # produce plot
    plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw() + ggplot2::ylab("log-likelihood")
  }
  
  # complete plot
  plot1 <- plot1 + ggplot2::geom_line(ggplot2::aes_(x = ~x, y = ~y, colour = "col1"))
  plot1 <- plot1 + ggplot2::coord_cartesian(xlim = c(0,nrow(df)))
  plot1 <- plot1 + ggplot2::scale_x_continuous(expand = c(0,0))
  plot1 <- plot1 + ggplot2::scale_colour_manual(values = col)
  plot1 <- plot1 + ggplot2::guides(colour = "none")
  plot1 <- plot1 + ggplot2::xlab("iteration")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce MCMC autocorrelation plot
#'
#' @description Produce MCMC autocorrelation plot of alpha or log-likelihood
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K which value of K to produce the plot for
#' @param rung which value of K to produce the plot for. Defaults to the cold rung
#' @param param whether to produce trace plot of admixture parameter
#'   (\code{"alpha"}) or the log-likelihood (\code{"loglike"})
#' @param col colour of the trace
#'
#' @importFrom coda effectiveSize
#' @import ggplot2
#' @export

plot_acf <- function(proj, K = NULL, rung = NULL, param = "alpha", col = "black") {
  
  # check inputs
  assert_class(proj, "mavproject")
  if (!is.null(K)) {
    assert_single_pos_int(K)
  }
  if (!is.null(rung)) {
    assert_single_pos_int(rung)
  }
  assert_in(param, c("alpha", "loglike"))
  
  # get active set and check non-zero
  s <- proj$active_set
  assert_neq(s, 0, message = "no active parameter set")
  
  # split analysis between alpha and loglikelihood
  if (param=="alpha") {
    
    # set default K to first value with output
    null_output <- mapply(function(x) {is.null(x$raw$alpha)}, proj$output$single_set[[s]]$single_K)
    if (all(null_output)) {
      stop("no alpha output for active parameter set")
    }
    if (is.null(K)) {
      K <- which(!null_output)[1]
      message(sprintf("using K = %s by default", K))
    }
    
    # check output exists for chosen K
    assert_leq(K, length(proj$output$single_set[[s]]$single_K), message = "K out of bounds")
    alpha <- as.vector(proj$output$single_set[[s]]$single_K[[K]]$raw$alpha)
    if (is.null(alpha)) {
      stop(sprintf("no alpha output for K = %s of active set", K))
    }
    
    # store variable to plot
    v <- alpha
    
  } else {
    
    # set default K to first value with output
    null_output <- mapply(function(x) {is.null(x$raw$loglike_sampling)}, proj$output$single_set[[s]]$single_K)
    if (all(null_output)) {
      stop("no loglike_sampling output for active parameter set")
    }
    if (is.null(K)) {
      K <- which(!null_output)[1]
      message(sprintf("using K = %s by default", K))
    }
    
    # check output exists for chosen K
    assert_leq(K, length(proj$output$single_set[[s]]$single_K), message = "K out of bounds")
    loglike_sampling <- proj$output$single_set[[s]]$single_K[[K]]$raw$loglike_sampling
    if (is.null(loglike_sampling)) {
      stop(sprintf("no loglike_sampling output for K = %s of active set", K))
    }
    
    # use cold rung by default
    rung <- define_default(rung, ncol(loglike_sampling))
    assert_leq(rung, ncol(loglike_sampling))
    loglike <- as.vector(loglike_sampling[,rung])
    
    # store variable to plot
    v <- loglike
  }
  
  # get autocorrelation
  lag_max <- round(3 * length(v) / coda::effectiveSize(v))
  lag_max <- max(lag_max, 20)
  lag_max <- min(lag_max, length(v))
  
  # get into ggplot format
  a <- acf(v, lag.max = lag_max, plot = FALSE)
  acf <- as.vector(a$acf)
  df <- data.frame(lag = (1:length(acf))-1, ACF = acf)
  
  # produce plot
  plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw()
  plot1 <- plot1 + ggplot2::geom_segment(ggplot2::aes_(x = ~lag, y = 0, xend = ~lag, yend = ~ACF, colour = "col1"))
  plot1 <- plot1 + ggplot2::scale_colour_manual(values = col)
  plot1 <- plot1 + ggplot2::guides(colour = "none")
  plot1 <- plot1 + ggplot2::xlab("lag") + ggplot2::ylab("ACF")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce MCMC density plot
#'
#' @description Produce MCMC density plot of alpha or log-likelihood
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K which value of K to produce the plot for
#' @param rung which value of K to produce the plot for. Defaults to the cold rung
#' @param param whether to produce trace plot of admixture parameter
#'   (\code{"alpha"}) or the log-likelihood (\code{"loglike"})
#' @param col colour of the trace
#'
#' @import ggplot2
#' @export

plot_density <- function(proj, K = NULL, rung = NULL, param = "alpha", col = "black") {
  
  # check inputs
  assert_class(proj, "mavproject")
  if (!is.null(K)) {
    assert_single_pos_int(K)
  }
  if (!is.null(rung)) {
    assert_single_pos_int(rung)
  }
  assert_in(param, c("alpha", "loglike"))
  
  # get active set and check non-zero
  s <- proj$active_set
  assert_neq(s, 0, message = "no active parameter set")
  
  # split analysis between alpha and loglikelihood
  if (param == "alpha") {
    
    # set default K to first value with output
    null_output <- mapply(function(x) {is.null(x$raw$alpha)}, proj$output$single_set[[s]]$single_K)
    if (all(null_output)) {
      stop("no alpha output for active parameter set")
    }
    if (is.null(K)) {
      K <- which(!null_output)[1]
      message(sprintf("using K = %s by default", K))
    }
    
    # check output exists for chosen K
    assert_leq(K, length(proj$output$single_set[[s]]$single_K), message = "K out of bounds")
    alpha <- as.vector(proj$output$single_set[[s]]$single_K[[K]]$raw$alpha)
    if (is.null(alpha)) {
      stop(sprintf("no alpha output for K = %s of active set", K))
    }
    
    # get into ggplot format
    df <- data.frame(v = alpha)
    
    # produce plot
    plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw() + ggplot2::xlab("alpha")
    
  } else {
    
    # set default K to first value with output
    null_output <- mapply(function(x) {is.null(x$raw$loglike_sampling)}, proj$output$single_set[[s]]$single_K)
    if (all(null_output)) {
      stop("no loglike_sampling output for active parameter set")
    }
    if (is.null(K)) {
      K <- which(!null_output)[1]
      message(sprintf("using K = %s by default", K))
    }
    
    # check output exists for chosen K
    assert_leq(K, length(proj$output$single_set[[s]]$single_K), message = "K out of bounds")
    loglike_sampling <- proj$output$single_set[[s]]$single_K[[K]]$raw$loglike_sampling
    if (is.null(loglike_sampling)) {
      stop(sprintf("no loglike_sampling output for K = %s of active set", K))
    }
    
    # use cold rung by default
    rung <- define_default(rung, ncol(loglike_sampling))
    assert_leq(rung, ncol(loglike_sampling))
    loglike <- as.vector(loglike_sampling[,rung])
    
    # get into ggplot format
    df <- data.frame(v = loglike)
    
    # produce plot
    plot1 <- ggplot2::ggplot(df) + ggplot2::theme_bw() + ggplot2::xlab("log-likelihood")
  }
  
  # produce plot
  #plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + ggplot2::geom_histogram(ggplot2::aes_(x = ~v, y = ~..density.., fill = "col1"), bins = 50)
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = col)
  plot1 <- plot1 + ggplot2::guides(fill = "none")
  plot1 <- plot1 + ggplot2::ylab("density")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce diagnostic plots of parameter alpha
#'
#' @description Produce diagnostic plots of parameter alpha
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K which value of K to produce the plot for
#' @param col colour of the trace
#'
#' @import ggplot2
#' @export

plot_alpha <- function(proj, K = NULL, col = "black") {
  
  # produce individual diagnostic plots and add features
  plot1 <- plot_trace(proj, K = K, col = col)
  plot1 <- plot1 + ggplot2::ggtitle("MCMC trace")
  
  plot2 <- plot_acf(proj, K = K, col = col)
  plot2 <- plot2 + ggplot2::ggtitle("autocorrelation")
  
  plot3 <- plot_density(proj, K = K, col = col)
  plot3 <- plot3 + ggplot2::ggtitle("density")
  
  # produce grid of plots
  ret <- gridExtra::grid.arrange(plot1, plot2, plot3, layout_matrix = rbind(c(1,1), c(2,3)))
}

#------------------------------------------------
#' @title Produce diagnostic plots of log-likelihood
#'
#' @description Produce diagnostic plots of the log-likelihood.
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K which value of K to produce the plot for
#' @param col colour of the trace
#'
#' @import ggplot2
#' @export

plot_loglike <- function(proj, K = NULL, col = "black") {
  
  # produce individual diagnostic plots and add features
  plot1 <- plot_trace(proj, K = K, param = "loglike", col = col)
  plot1 <- plot1 + ggplot2::ggtitle("MCMC trace")
  
  plot2 <- plot_acf(proj, K = K, param = "loglike", col = col)
  plot2 <- plot2 + ggplot2::ggtitle("autocorrelation")
  
  plot3 <- plot_density(proj, K = K, param = "loglike", col = col)
  plot3 <- plot3 + ggplot2::ggtitle("density")
  
  # produce grid of plots
  ret <- gridExtra::grid.arrange(plot1, plot2, plot3, layout_matrix = rbind(c(1,1), c(2,3)))
}

#------------------------------------------------
#' @title Plot Metropolis coupling acceptance rates
#'
#' @description Plot Metropolis coupling acceptance rates between all adjacent
#'   rungs. It is important that acceptance rates are greater than zero between
#'   all pairs of rungs in order to give confidence that the MCMC has thoroughly
#'   explored the full posterior distribution. If there are any points in this
#'   plot where acceptance rates drop close to zero then consider changing the
#'   number of rungs or their distribution (see \code{rungs} and \code{GTI_pow}
#'   arguments in \code{run_mcmc()}, respectively).
#'
#' @param proj an rmaverick project, as produced by the function 
#'   \code{mavproject()}
#' @param K which value of K to produce the plot for
#' @param x_axis_type how to format the x-axis. 1 = integer rungs, 2 = values of
#'   the thermodynamic power.
#'
#' @importFrom grDevices grey
#' @export

plot_mc_acceptance <- function(proj, K = NULL, x_axis_type = 1) {
  
  # check inputs
  assert_class(proj, "mavproject")
  if (!is.null(K)) {
    assert_single_pos_int(K)
  }
  assert_single_pos_int(x_axis_type)
  assert_in(x_axis_type, 1:2)
  
  # get active set and check non-zero
  s <- proj$active_set
  assert_neq(s, 0, message = "no active parameter set")
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$raw$coupling_accept)}, proj$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no coupling acceptance rate output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  assert_leq(K, length(proj$output$single_set[[s]]$single_K), message = "K out of bounds")
  mc_accept <- proj$output$single_set[[s]]$single_K[[K]]$raw$coupling_accept
  if (is.null(mc_accept)) {
    stop(sprintf("no metropolis coupling output for K = %s of active set", K))
  }
  
  # get acceptance rates and thermodynamic powers
  GTI_pow <- proj$output$single_set[[s]]$single_K[[K]]$function_call$args$GTI_pow
  rungs <- length(mc_accept) + 1
  thermo_power <- seq(0, 1.0, l = rungs)^GTI_pow
  thermo_power_mid <- thermo_power[-1] - diff(thermo_power)/2
  
  # exit if rungs = 1
  if (rungs == 1) {
    stop("no metropolis coupling when rungs = 1")
  }
  
  # define x-axis type
  if (x_axis_type == 1) {
    breaks_vec <- 1:rungs
    x_vec <- (2:rungs) - 0.5
    x_lab <- "rung"
  } else {
    breaks_vec <- thermo_power
    x_vec <- thermo_power_mid
    x_lab <- "thermodynamic power"
  }
  
  # get data into ggplot format and define temperature colours
  df <- data.frame(x_vec = x_vec, mc_accept = mc_accept, col = thermo_power_mid)
  
  # produce plot
  plot1 <- ggplot2::ggplot(df) + 
    ggplot2::geom_vline(ggplot2::aes(xintercept = .data$x), col = grey(0.9), data = data.frame(x = breaks_vec)) +
    ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
    ggplot2::geom_point(ggplot2::aes(x = x_vec, y = mc_accept, color = col)) + 
    ggplot2::xlab(x_lab) + ggplot2::ylab("coupling acceptance rate") + 
    ggplot2::scale_colour_gradientn(colours = c("red", "blue"), name = "thermodynamic\npower", limits = c(0,1)) +
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())
  
  return(plot1)
}
