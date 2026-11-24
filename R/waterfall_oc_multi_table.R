#' Assemble the Cross-Scenario Summary Table of the Waterfall Design
#'
#' @description
#'   Build the data frame shown by \code{\link{print.waterfall_oc_multi}} from
#'   the per-scenario results of \code{\link{sim_waterfall_multi}}, one row per
#'   scenario. Internal.
#'
#' @param results
#'   List of \code{waterfall_oc} objects, one per scenario.
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
waterfall_oc_multi_table <- function(results, scenario_names, digits = 1) {

  rows <- lapply(seq_along(results), function(i) {
    res <- results[[i]]
    contour <- if (nrow(res$true_mtd) > 0L) {
      paste0("A", res$true_mtd[, "dose_a"], "B", res$true_mtd[, "dose_b"],
             collapse = ", ")
    } else {
      "none"
    }
    data.frame(
      Scenario = scenario_names[i],
      "MTD contour" = contour,
      "Contour selected (%)" = round(res$pct_contour_selected, digits),
      "Patients on contour (%)" = round(res$pct_patients_contour, digits),
      "Patients above contour (%)" = round(res$pct_patients_above, digits),
      "No MTD (%)" = round(res$percent_no_mtd, digits),
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
