#' Validate the Starting Dose Combination
#'
#' @description
#'   Check that \code{start_dose} gives a dose level of drug A and of drug B
#'   that exist in \code{p_true}. Internal.
#'
#' @param start_dose
#'   The object to check.
#'
#' @param p_true
#'   Matrix of true DLT probabilities, already validated.
#'
#' @return
#'   \code{TRUE}, invisibly.
#'
#' @keywords internal
#' @noRd
check_start_combination <- function(start_dose, p_true) {
  if (!is.numeric(start_dose) || length(start_dose) != 2L ||
      any(!is.finite(start_dose)) || any(start_dose != round(start_dose)) ||
      any(start_dose < 1)) {
    stop("'start_dose' must be two whole numbers of at least 1, the dose ",
         "levels of drug A and drug B", call. = FALSE)
  }
  if (start_dose[1L] > nrow(p_true) || start_dose[2L] > ncol(p_true)) {
    stop("'start_dose' must not exceed the number of levels of either drug ",
         "in 'p_true'", call. = FALSE)
  }
  invisible(TRUE)
}
