#' Collect the Design Results Passed to oc_results() and oc_compare()
#'
#' @description
#'   Internal helper. Apply \code{oc_extract()} to every design result, with
#'   the argument names as design names, and check that the designs can be
#'   set side by side.
#'
#' @param args List of the design results.
#'
#' @return A list of the results of \code{oc_extract()}.
#'
#' @noRd
oc_collect <- function(args) {
  if (length(args) == 0L) {
    stop("give at least one design result", call. = FALSE)
  }
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  sets <- lapply(seq_along(args), function(i) oc_extract(args[[i]], nms[i]))
  designs <- vapply(sets, function(z) z$design, character(1))
  if (anyDuplicated(designs)) {
    stop("design names must be unique; name the arguments, for example ",
         "oc_results(A = res_a, B = res_b)", call. = FALSE)
  }
  types <- vapply(sets, function(z) z$type, character(1))
  if (length(unique(types)) > 1L) {
    stop("the designs must be all for a single agent or all for a drug ",
         "combination", call. = FALSE)
  }
  dims <- lapply(sets, function(z) dim(as.matrix(z$scenarios[[1L]]$sel)))
  if (length(unique(dims)) > 1L) {
    stop("the designs must have the same number of doses", call. = FALSE)
  }
  sets
}
