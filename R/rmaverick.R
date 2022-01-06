
#------------------------------------------------
#' @title Analysis of genetic population structure via thermodynamic integration.
#'
#' @description Mixture models provide a powerful way of teasing apart
#'   population structure from genetic data, as exemplified by the program
#'   STRUCTURE. However, choosing the number of sub-populations (denoted K) is a
#'   more difficult problem as this requires evaluation of the *evidence* of
#'   each model. In contrast to other programs, *rmaverick* uses thermodynamic
#'   integration to arrive at precise and accurate estimates of the evidence for
#'   each *K*, allowing the user to choose appropriate values of this parameter,
#'   and also to compare between different evolutionary models (e.g.
#'   with/without admixture).
#'
#' @docType package
#' @name rmaverick
NULL

#------------------------------------------------
#' @useDynLib rmaverick, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("rmaverick", libpath)  # nocov
}