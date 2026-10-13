#' Assemble the Cross-Scenario Summary Table for Drug Combinations
#'
#' @description
#'   Build the data frame shown by \code{\link{print.comb_boin_oc_multi}} from
#'   the per-scenario results of \code{\link{sim_comb_boin_multi}}, one row per
#'   scenario. Internal.
#'
#' @param results
#'   List of \code{comb_boin_oc} objects, one per scenario.
#'
#' @param scenario_names
#'   Character vector of scenario names.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @return
#'   A data frame with one row per scenario.
#'
#' @keywords internal
#' @noRd
comb_oc_multi_table <- function(results, scenario_names, digits = 1) {

  rows <- lapply(seq_along(results), function(i) {
    res <- results[[i]]
    data.frame(
      Scenario = scenario_names[i],
      "True MTD" = paste0("A", res$true_mtd[, "dose_a"], "B",
                          res$true_mtd[, "dose_b"], collapse = ", "),
      "MTD selected (%)" = round(res$pct_true_mtd_selected, digits),
      "Patients at MTD (%)" = round(res$pct_patients_true_mtd, digits),
      "No MTD (%)" = round(res$percent_no_mtd, digits),
      "Above cutoff selected (%)" = round(res$overdose$pct_trials_mtd_above,
                                          digits),
      "Patients above cutoff (%)" = round(res$overdose$pct_patients, digits),
      "Patients" = round(res$total_n_pts, digits),
      "DLTs" = round(res$total_n_tox, digits),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
