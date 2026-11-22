#' Total Scores as a Matrix
#'
#' @description
#'   Coerce the total scores of one trial (a vector) or of several trials (a
#'   matrix with one row per trial) to a numeric matrix, and check that they
#'   are non-negative and finite. Internal.
#'
#' @param x
#'   Numeric vector or matrix.
#'
#' @param name
#'   Character scalar. Name of the argument, for the error message.
#'
#' @return
#'   A numeric matrix with one row per trial.
#'
#' @keywords internal
#' @noRd
as_score_matrix <- function(x, name) {
  if (!is.matrix(x)) x <- matrix(x, nrow = 1L)
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0)) {
    stop("'", name, "' must hold non-negative finite numbers", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}
