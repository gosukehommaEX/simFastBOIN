#' Suspension Rules of the TITE-BOIN Design
#'
#' @description
#'   Internal helper shared by \code{\link{tite_boin_decision_table}} and
#'   \code{\link{tite_boin_simulate}}. Fill in the defaults of the suspension
#'   rules for the chosen method and validate them.
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
#' @param min_follow_up
#'   Numeric scalar between 0 and 1. Shortest follow-up, as a fraction of the
#'   assessment window, that every pending patient at the current dose must have
#'   reached before the dose is escalated. Defaults to 0.
#'
#' @return
#'   A list with components \code{max_pending_ratio}, \code{min_completed} (an
#'   integer) and \code{min_follow_up}.
#'
#' @noRd
tite_rules <- function(method, max_pending_ratio, min_completed,
                       min_follow_up = 0) {

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
  if (!is.numeric(min_follow_up) || length(min_follow_up) != 1L ||
      !is.finite(min_follow_up) || min_follow_up < 0 || min_follow_up > 1) {
    stop("'min_follow_up' must be a single number between 0 and 1",
         call. = FALSE)
  }

  list(max_pending_ratio = max_pending_ratio,
       min_completed = as.integer(min_completed),
       min_follow_up = as.numeric(min_follow_up))
}
