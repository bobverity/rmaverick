
#------------------------------------------------
# is NULL
assert_null <- function(x, name = deparse(substitute(x))) {
  if (!is.null(x)) {
    stop(sprintf("%s must be null", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is not NULL
assert_non_null <- function(x, name = deparse(substitute(x))) {
  if (is.null(x)) {
    stop(sprintf("%s cannot be null", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is character
assert_character <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x)) {
    stop(sprintf("%s must be character", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is logical
assert_logical <- function(x, name = deparse(substitute(x))) {
  if (!is.logical(x)) {
    stop(sprintf("%s must be logical", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is single logical
assert_scalar_logical <- function(x, name = deparse(substitute(x))) {
  assert_logical(x, name)
  assert_length(x, 1, name)
  return(TRUE)
}

#------------------------------------------------
# is numeric
assert_numeric <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x)) {
    stop(sprintf("%s must be numeric", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is positive (with or without zero allowed)
assert_pos <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_numeric(x, name)
  if (zero_allowed) {
    if (!all(x>=0)) {
      stop(sprintf("%s must be greater than or equal to zero", name), call. = FALSE)
    }
  } else {
    if (!all(x>0)) {
      stop(sprintf("%s must be greater than zero", name), call. = FALSE)
    }
  }
  return(TRUE)
}

#------------------------------------------------
# is single positive value (with or without zero allowed)
assert_scalar_pos <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_pos(x, zero_allowed, name)
  assert_length(x, 1, name)
  return(TRUE)
}

#------------------------------------------------
# is integer
assert_int <- function(x, name = deparse(substitute(x))) {
  assert_numeric(x, name)
  if (!isTRUE(all.equal(x, as.integer(x)))) {
    stop(sprintf("%s must be integer valued", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is single integer
assert_scalar_int <- function(x, name = deparse(substitute(x))) {
  assert_int(x, name)
  assert_length(x, 1, name)
  return(TRUE)
}

#------------------------------------------------
# is positive integer (with or without zero allowed)
assert_pos_int <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_int(x, name)
  assert_pos(x, zero_allowed, name)
  return(TRUE)
}

#------------------------------------------------
# is single positive integer (with or without zero allowed)
assert_scalar_pos_int <- function(x, zero_allowed = TRUE, name = deparse(substitute(x))) {
  assert_pos_int(x, zero_allowed, name)
  assert_length(x, 1, name)
  return(TRUE)
}

#------------------------------------------------
# x and y are equal (accepts any type)
assert_eq <- function(x, y, name = deparse(substitute(x))) {
  if (!isTRUE(all.equal(x,y))) {
    stop(sprintf("%s must equal %s", name, y), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is greater than
assert_gr <- function(x, y, name = deparse(substitute(x))) {
  assert_numeric(x, name)
  if (!all(x>y)) {
    stop(sprintf("%s must be greater than %s", name, y), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is greater than or equal to
assert_greq <- function(x, y, name = deparse(substitute(x))) {
  assert_numeric(x, name)
  if (!all(x>=y)) {
    stop(sprintf("%s must be greater than or equal to %s", name, y), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is less than
assert_le <- function(x, y, name = deparse(substitute(x))) {
  assert_numeric(x, name)
  if (!all(x<y)) {
    stop(sprintf("%s must be less than %s", name, y), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is less than or equal to
assert_leq <- function(x, y, name = deparse(substitute(x))) {
  assert_numeric(x, name)
  if (!all(x<=y)) {
    stop(sprintf("%s must be less than or equal to %s", name, y), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is between bounds (inclusive or exclusive)
assert_bounded <- function(x, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE, name = deparse(substitute(x))) {
  assert_numeric(x, name)
  if (inclusive_left) {
    if (!all(x>=left)) {
      stop(sprintf("%s must be greater than or equal to %s", name, left), call. = FALSE)
    }
  } else {
    if (!all(x>left)) {
      stop(sprintf("%s must be greater than %s", name, left), call. = FALSE)
    }
  }
  if (inclusive_right) {
    if (!all(x<=right)) {
      stop(sprintf("%s must be less than or equal to %s", name, right), call. = FALSE)
    }
  } else {
    if (!all(x<right)) {
      stop(sprintf("%s must be less than %s", name, right), call. = FALSE)
    }
  }
  return(TRUE)
}

#------------------------------------------------
# value matches one of several strings
assert_in <- function(x, s, name = deparse(substitute(x))) {
  if (!all(x %in% s)) {
    stop(sprintf("%s must be one of {%s}", name, paste(s, collapse=", ")))
  }
  return(TRUE)
}

#------------------------------------------------
# no duplicated values
assert_noduplicates <- function(x, name = deparse(substitute(x))) {
  if (any(duplicated(x))) {
    stop(sprintf("%s cannot contain duplicated values", name))
  }
  return(TRUE)
}

#------------------------------------------------
# length equal to given value
assert_length <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) != n) {
    stop(sprintf("%s must be of length %s", name, n), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# objects all same length
assert_same_length <- function(...) {
  l <- mapply(length, list(...))
  if (!length(unique(l)) == 1) {
    dots <- match.call(expand.dots = FALSE)$...
    dot_names <- paste(sapply(dots, deparse), collapse=", ")
    stop(sprintf("variables %s must be the same length", dot_names), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is matrix
assert_matrix <- function(x, name = deparse(substitute(x))) {
  if (!is.matrix(x)) {
    stop(sprintf("%s must be a matrix", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is square matrix
assert_square_matrix <- function(x, name = deparse(substitute(x))) {
  assert_matrix(x, name)
  if (nrow(x) != ncol(x)) {
    stop(sprintf("%s must be a square matrix", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is symmetric matrix
assert_symmetric_matrix <- function(x, name = deparse(substitute(x))) {
  assert_square_matrix(x, name)
  if (!isSymmetric(x)) {
    stop(sprintf("%s must be a symmetric matrix", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is data frame
assert_dataframe <- function(x, name = deparse(substitute(x))) {
  if (!is.data.frame(x)) {
    stop(sprintf("%s must be a data frame", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is class mavproject
assert_mavproject <- function(x, name = deparse(substitute(x))) {
  if (!is.mavproject(x)) {
    stop(sprintf("%s must be of class 'mavproject'", name), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is class cluster
assert_cluster <- function(x, name = deparse(substitute(x))) {
  if (!is.cluster(x)) {
    stop(sprintf("%s must be of class 'cluster'", name), call. = FALSE)
  }
  return(TRUE)
}


