#' Suspension Rules of the TITE-BOIN Design
#'
#' @description
#'   Internal helper shared by \code{\link{tite_boin_decision_table}} and
#'   \code{\link{tite_boin_simulate}}. Fill in the defaults of the two
#'   suspension rules for the chosen method and validate them.
#'
#' @param method
#'   Character scalar, \code{"imputation"} or \code{"ess"}.
#'
#' @param max_pending_ratio
#'   Numeric scalar or \code{NULL}. \code{NULL} gives 0.5 for
#'   \code{"imputation"} and 1 for \code{"ess"}.
#'
#' @param min_completed
#'   Integer scalar or \code{NULL}. \code{NULL} gives 0 for \code{"imputation"}
#'   and 2 for \code{"ess"}.
#'
#' @return
#'   A list with components \code{max_pending_ratio} and \code{min_completed},
#'   the latter an integer.
#'
#' @noRd
tite_rules <- function(method, max_pending_ratio, min_completed) {

  if (is.null(max_pending_ratio)) {
    max_pending_ratio <- if (method == "imputation") 0.5 else 1
  }
  if (is.null(min_completed)) {
    min_completed <- if (method == "imputation") 0L else 2L
  }

  if (!is.numeric(max_pending_ratio) || length(max_pending_ratio) != 1L ||
      !is.finite(max_pending_ratio) || max_pending_ratio <= 0 ||
      max_pending_ratio > 1) {
    stop("'max_pending_ratio' must be a single number greater than 0 and at most 1",
         call. = FALSE)
  }
  check_count(min_completed, "min_completed", 0L)

  list(max_pending_ratio = max_pending_ratio,
       min_completed = as.integer(min_completed))
}
