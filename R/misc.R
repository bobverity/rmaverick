
#------------------------------------------------
# replace NULL value with default
# (not exported)

define_default <- function(x, default_value) {
  if (is.null(x)) {
    x <- default_value
  }
  return(x)
}
