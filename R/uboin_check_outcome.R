#' Validate a Matrix of Outcome Probabilities of U-BOIN
#'
#' @description
#'   Check that every row of \code{p_outcome} holds four probabilities that
#'   add up to one. Internal.
#'
#' @param p_outcome
#'   The object to check.
#'
#' @return
#'   \code{p_outcome} as a numeric matrix with the names of the outcomes.
#'
#' @keywords internal
#' @noRd
uboin_check_outcome <- function(p_outcome) {
  if (!is.matrix(p_outcome) || !is.numeric(p_outcome) ||
      ncol(p_outcome) != 4L || nrow(p_outcome) < 1L ||
      any(!is.finite(p_outcome)) || any(p_outcome < 0)) {
    stop("'p_outcome' must be a numeric matrix with one row per dose and ",
         "four columns of non-negative probabilities", call. = FALSE)
  }
  if (any(abs(rowSums(p_outcome) - 1) > 1e-8)) {
    stop("each row of 'p_outcome' must add up to 1", call. = FALSE)
  }
  p_outcome <- matrix(as.numeric(p_outcome), nrow = nrow(p_outcome))
  colnames(p_outcome) <- uboin_outcome_names()
  p_outcome
}
