#' Stopping Reasons of BOIN-ETC Trials
#'
#' @description
#'   The stopping reasons of \code{\link{boin_etc_simulate}}, in the order of
#'   the stop codes of the C++ engine. Internal.
#'
#' @return A character vector.
#'
#' @keywords internal
#' @noRd
boin_etc_stop_reasons <- function() {
  c("lowest_excluded", "subtrials_completed", "max_cohorts")
}
