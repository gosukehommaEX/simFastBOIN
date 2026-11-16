#' Stopping Reasons of the BOIN12 Simulations
#'
#' @return Character vector of the reasons, in the order of the stop codes of
#'   the C++ engine.
#'
#' @noRd
boin12_stop_reasons <- function() {
  c("no_admissible_dose", "n_earlystop", "max_cohorts")
}
