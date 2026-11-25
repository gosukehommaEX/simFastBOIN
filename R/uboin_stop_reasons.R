#' Stopping Reasons of U-BOIN Trials
#'
#' @description
#'   The reasons a simulated U-BOIN trial ends, in the order of the stop codes
#'   of the C++ core. Internal.
#'
#' @return
#'   A character vector.
#'
#' @keywords internal
#' @noRd
uboin_stop_reasons <- function() {
  c("max_sample_size", "n_stop", "lowest_dose_eliminated",
    "no_admissible_dose")
}
