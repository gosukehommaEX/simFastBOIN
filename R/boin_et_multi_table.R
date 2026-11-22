#' Assemble the Cross-Scenario Summary Table of BOIN-ET
#'
#' @description
#'   Build the data frame shown by \code{\link{print.boin_et_oc_multi}} from
#'   the per-scenario results of \code{\link{sim_boin_et_multi}}. Internal.
#'
#' @param results
#'   List of \code{boin_et_oc} objects, one per scenario.
#'
#' @param scenario_names
#'   Character vector of scenario names.
#'
#' @param n_doses
#'   Integer scalar. Number of dose levels.
#'
#' @param percent
#'   Logical scalar. Express the average numbers of patients, DLTs and
#'   responses at each dose as percentages of the trial totals. Defaults to
#'   \code{FALSE}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @return
#'   A data frame with six rows per scenario.
#'
#' @keywords internal
#' @noRd
boin_et_multi_table <- function(results, scenario_names, n_doses,
                                percent = FALSE, digits = 1) {

  dose_names <- paste0("DL", seq_len(n_doses))
  lab <- et_oc_labels(results[[1L]]$design)
  labels <- c("Patients treated", lab$tox, lab$resp)
  if (percent) labels <- paste(labels, "(%)")
  item_labels <- c(lab$rate, lab$rate_resp, "OBD selected (%)", labels)

  blocks <- lapply(seq_along(results), function(i) {
    res <- results[[i]]
    pts <- res$n_pts_dose
    tox <- res$n_tox_dose
    resp <- res$n_resp_dose
    if (percent) {
      pts <- pts / res$total_n_pts * 100
      tox <- tox / res$total_n_tox * 100
      resp <- resp / res$total_n_resp * 100
    }
    values <- rbind(res$p_true * 100, res$p_resp * 100, res$sel_percent,
                    pts, tox, resp)
    dimnames(values) <- list(NULL, dose_names)
    block <- data.frame(
      Scenario = c(scenario_names[i], rep("", 5L)),
      Item = item_labels,
      stringsAsFactors = FALSE
    )
    block <- cbind(block, as.data.frame(round(values, digits)))
    block[["Total / No OBD"]] <- c(
      NA_real_,
      NA_real_,
      round(res$percent_no_obd, digits),
      round(res$total_n_pts, digits),
      round(res$total_n_tox, digits),
      round(res$total_n_resp, digits)
    )
    block
  })

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}
