
#------------------------------------------------
# test if integer
# (not exported)

is.int <- function(x) {
  as.integer(x)==x
}

#------------------------------------------------
# replace NULL value with default
# (not exported)

define_default <- function(x, default_value) {
  if (is.null(x)) {
    x <- default_value
  }
  return(x)
}

# -----------------------------------
# user_yes_no
# ask user a yes/no question. Return TRUE/FALSE.
# (not exported)

user_yes_no <- function(x="continue? (Y/N): ") {
  
  user_choice <- NA
  while (!user_choice %in% c("Y", "y" ,"N", "n")) {
    user_choice <- readline(x)
  }
  return(user_choice %in% c("Y", "y"))
}

# -----------------------------------
# rdirichlet
# draw from dirichlet distribution with vector of parameter inputs
# (not exported)

rdirichlet <- function(alpha_vec) {
  z <- rgamma(length(alpha_vec), shape=alpha_vec, scale=1)
  ret <- z/sum(z)
  return(ret)
}

# -----------------------------------
# mat_to_Rcpp
# takes matrix as input, converts to list format for use within Rcpp code
# (not exported)

mat_to_rcpp <- function(x) {
  return(split(x, f=1:nrow(x)))
}

# -----------------------------------
# Rcpp_to_mat
# Takes list format returned from Rcpp and converts to matrix.
# (not exported)

rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  return(ret)
}

#------------------------------------------------
# define default colours for barplots
# (not exported)

default_colours <- function(K) {
  
  # generate palette and colours
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  my_palette <- colorRampPalette(raw_cols)
  ret <- my_palette(max(K,6))
  
  # if fewer than 6 colours then choose manually
  if (K==5) {
    ret = ret[c(1,2,3,5,6)]
  } else if (K==4) {
    ret = ret[c(1,2,3,5)]
  } else if (K==3) {
    ret = ret[c(1,3,5)]
  } else if (K==2) {
    ret = ret[c(1,5)]
  } else if (K==1) {
    ret = ret[1]
  }
  
  return(ret)
}