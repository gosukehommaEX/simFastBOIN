#' Check a Matrix of Response Probabilities
#'
#' @param p_resp Numeric matrix of response probabilities.
#' @param p_true Numeric matrix of DLT probabilities, already checked.
#'
#' @return \code{TRUE}, invisibly, or an error.
#'
#' @noRd
check_p_resp_matrix <- function(p_resp, p_true) {
  if (!is.matrix(p_resp) || !is.numeric(p_resp) ||
      !identical(dim(p_resp), dim(p_true)) || any(!is.finite(p_resp)) ||
      any(p_resp < 0) || any(p_resp > 1)) {
    stop("'p_resp' must be a matrix of response probabilities between 0 ",
         "and 1 of the same size as 'p_true'", call. = FALSE)
  }
  invisible(TRUE)
}
