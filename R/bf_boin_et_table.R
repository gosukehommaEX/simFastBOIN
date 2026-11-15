#' Table of the Operating Characteristics of the BF-BOIN-ET Design
#'
#' @description
#'   The rows shown by \code{\link{print.bf_boin_et_oc}} and, for every
#'   scenario, by \code{\link{print.bf_boin_et_oc_multi}}.
#'
#' @param x A \code{bf_boin_et_oc} object.
#' @param percent Logical scalar. Express the averages at each dose as
#'   percentages of the trial totals.
#'
#' @return A data frame with a column \code{Item}, one column per dose and a
#'   column \code{"Total / No OBD"}.
#'
#' @noRd
bf_boin_et_table <- function(x, percent) {
  dose_names <- paste0("DL", seq_len(x$n_doses))
  pts <- x$n_pts_dose
  bf <- x$n_bf_dose
  tox <- x$n_tox_dose
  resp <- x$n_resp_dose
  labels <- c("Patients treated", "Patients backfilled", "Patients with DLT",
              "Patients with response")
  if (percent) {
    pts <- pts / x$total_n_pts * 100
    bf <- bf / x$total_n_bf * 100
    tox <- tox / x$total_n_tox * 100
    resp <- resp / x$total_n_resp * 100
    labels <- paste(labels, "(%)")
  }
  values <- rbind(x$p_true * 100, x$p_resp * 100, x$sel_percent, pts, bf, tox,
                  resp)
  dimnames(values) <- list(NULL, dose_names)
  out <- data.frame(
    Item = c("True DLT rate (%)", "True response rate (%)", "OBD selected (%)",
             labels),
    stringsAsFactors = FALSE
  )
  out <- cbind(out, as.data.frame(values))
  out[["Total / No OBD"]] <- c(NA_real_, NA_real_, x$percent_no_obd,
                               x$total_n_pts, x$total_n_bf, x$total_n_tox,
                               x$total_n_resp)
  out
}
