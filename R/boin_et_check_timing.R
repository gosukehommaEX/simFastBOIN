#' Check the Time-to-Event Settings of the BOIN-ET Designs
#'
#' @description
#'   Validate the assessment windows, the accrual rate and the shares of late
#'   events of \code{\link{tite_boin_et_simulate}} and
#'   \code{\link{sim_bf_boin_et}}.
#'
#' @param window,resp_window,accrual_rate Positive numeric scalars.
#' @param late_fraction,resp_late_fraction Numeric scalars strictly between 0
#'   and 1.
#'
#' @return \code{TRUE}, invisibly.
#'
#' @noRd
boin_et_check_timing <- function(window, resp_window, accrual_rate,
                                 late_fraction, resp_late_fraction) {
  for (nm in c("window", "resp_window", "accrual_rate")) {
    value <- get(nm)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) {
      stop("'", nm, "' must be a single positive number", call. = FALSE)
    }
  }
  check_scalar_prob(late_fraction, "late_fraction")
  check_scalar_prob(resp_late_fraction, "resp_late_fraction")
  invisible(TRUE)
}
