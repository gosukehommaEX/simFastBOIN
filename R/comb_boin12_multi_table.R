#' Summary Table of Comb-BOIN12 Results Across Scenarios
#'
#' @param results Named list of \code{comb_boin12_oc} objects.
#' @param scenario_names Character vector of scenario names.
#' @param digits Integer scalar. Decimal places.
#'
#' @return A data frame with one row per scenario.
#'
#' @noRd
comb_boin12_multi_table <- function(results, scenario_names, digits = 1) {

  rows <- lapply(seq_along(results), function(i) {
    res <- results[[i]]
    data.frame(
      Scenario = scenario_names[i],
      "True OBD" = if (nrow(res$true_obd) == 0L) {
        "none"
      } else {
        paste0("A", res$true_obd[, "dose_a"], "B", res$true_obd[, "dose_b"],
               collapse = ", ")
      },
      "OBD selected (%)" = round(res$pct_true_obd_selected, digits),
      "Patients at OBD (%)" = round(res$pct_patients_true_obd, digits),
      "No OBD (%)" = round(res$percent_no_obd, digits),
      "Above limit selected (%)" = round(res$overdose$pct_selected, digits),
      "Patients above limit (%)" = round(res$overdose$pct_patients, digits),
      "Patients" = round(res$total_n_pts, digits),
      "DLTs" = round(res$total_n_tox, digits),
      "Responses" = round(res$total_n_resp, digits),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
