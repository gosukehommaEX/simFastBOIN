#' Assemble the Cross-Scenario Summary Table
#'
#' @description
#'   Build the data frame shown by \code{\link{print.boin_oc_multi}} from the
#'   per-scenario results of \code{\link{sim_boin_multi}}. Internal.
#'
#' @param results
#'   List of \code{boin_oc} objects, one per scenario.
#'
#' @param scenario_names
#'   Character vector of scenario names.
#'
#' @param n_doses
#'   Integer scalar. Number of dose levels.
#'
#' @param percent
#'   Logical scalar. Express the average number of patients and DLTs at each dose
#'   as a percentage of the trial total instead of a count. Defaults to
#'   \code{FALSE}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @return
#'   A data frame with four rows per scenario.
#'
#' @keywords internal
#' @noRd
oc_multi_table <- function(results, scenario_names, n_doses, percent = FALSE,
                           digits = 1) {

  dose_names <- paste0("DL", seq_len(n_doses))
  pts_label <- if (percent) "Patients treated (%)" else "Patients treated"
  tox_label <- if (percent) "Patients with DLT (%)" else "Patients with DLT"
  item_labels <- c("True DLT rate (%)", "MTD selected (%)", pts_label, tox_label)

  blocks <- lapply(seq_along(results), function(i) {
    res <- results[[i]]
    pts <- res$n_pts_dose
    tox <- res$n_tox_dose
    if (percent) {
      pts <- pts / res$total_n_pts * 100
      tox <- tox / res$total_n_tox * 100
    }
    values <- rbind(res$p_true * 100, res$sel_percent, pts, tox)
    dimnames(values) <- list(NULL, dose_names)
    block <- data.frame(
      Scenario = c(scenario_names[i], "", "", ""),
      Item = item_labels,
      stringsAsFactors = FALSE
    )
    block <- cbind(block, as.data.frame(round(values, digits)))
    block[["Total / No MTD"]] <- c(
      NA_real_,
      round(res$percent_no_mtd, digits),
      round(res$total_n_pts, digits),
      round(res$total_n_tox, digits)
    )
    block
  })

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}
