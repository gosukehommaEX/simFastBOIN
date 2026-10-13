#' Coerce Drug-Combination Counts to a Trial Array
#'
#' @description
#'   Turn a matrix of counts for one trial, or an array of counts for several
#'   trials, into an integer array with dimensions trial by level of drug A by
#'   level of drug B. Internal.
#'
#' @param x
#'   Numeric matrix or three-dimensional array of counts.
#'
#' @param name
#'   Character scalar. Name of the argument, used in error messages.
#'
#' @return
#'   An integer array with three dimensions.
#'
#' @keywords internal
#' @noRd
as_comb_count_array <- function(x, name) {
  if (!is.numeric(x)) stop("'", name, "' must be numeric", call. = FALSE)
  dims <- dim(x)
  if (length(dims) == 2L) {
    x <- array(x, c(1L, dims))
  } else if (length(dims) != 3L) {
    stop("'", name, "' must be a matrix for one trial, or an array with ",
         "dimensions trial by level of drug A by level of drug B",
         call. = FALSE)
  }
  if (any(is.na(x))) stop("'", name, "' must not contain missing values", call. = FALSE)
  if (any(x < 0) || any(x != round(x))) {
    stop("'", name, "' must contain non-negative whole numbers", call. = FALSE)
  }
  storage.mode(x) <- "integer"
  x
}
