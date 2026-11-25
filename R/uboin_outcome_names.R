#' Names of the Four Outcomes of U-BOIN
#'
#' @description
#'   The four outcomes of a patient in the order of Zhou, Lee and Yuan (2019).
#'   Internal.
#'
#' @return
#'   A character vector of length four.
#'
#' @keywords internal
#' @noRd
uboin_outcome_names <- function() {
  c("no_resp_dlt", "no_resp_no_dlt", "resp_dlt", "resp_no_dlt")
}
