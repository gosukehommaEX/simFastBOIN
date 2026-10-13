#' Validate a Matrix of True DLT Probabilities
#'
#' @description
#'   Check the \code{p_true} argument of the drug-combination functions: a
#'   numeric matrix with one row per level of drug A and one column per level
#'   of drug B. Warns when toxicity does not increase along every row and
#'   column. Internal.
#'
#' @param p_true
#'   The object to check.
#'
#' @return
#'   \code{TRUE}, invisibly.
#'
#' @keywords internal
#' @noRd
check_p_true_matrix <- function(p_true) {
  if (!is.matrix(p_true) || !is.numeric(p_true) || length(p_true) < 2L ||
      any(!is.finite(p_true))) {
    stop("'p_true' must be a numeric matrix of finite DLT probabilities, ",
         "with one row per level of drug A and one column per level of ",
         "drug B, and at least two combinations", call. = FALSE)
  }
  if (any(p_true < 0) || any(p_true > 1)) {
    stop("'p_true' must contain values between 0 and 1", call. = FALSE)
  }
  rows_sorted <- all(apply(p_true, 1L, function(z) !is.unsorted(z)))
  cols_sorted <- all(apply(p_true, 2L, function(z) !is.unsorted(z)))
  if (!rows_sorted || !cols_sorted) {
    warning("'p_true' is not non-decreasing along every row and column; ",
            "the design assumes that toxicity increases with the dose of ",
            "each drug", call. = FALSE)
  }
  invisible(TRUE)
}
