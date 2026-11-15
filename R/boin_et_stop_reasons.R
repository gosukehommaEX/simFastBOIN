#' Stopping Reasons of the BOIN-ET and TITE-BOIN-ET Trials
#'
#' @description
#'   The values of \code{stop_reason} in \code{\link{boin_et_simulate}} and
#'   \code{\link{tite_boin_et_simulate}}, in the order of the stop codes of
#'   the C++ core.
#'
#' @return Character vector of the stopping reasons.
#'
#' @noRd
boin_et_stop_reasons <- function() {
  c("no_admissible_dose", "no_lower_admissible_dose", "n_earlystop",
    "max_cohorts", "top_dose_excluded")
}
