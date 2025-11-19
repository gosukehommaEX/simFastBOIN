#' Print BOIN Simulation Summary
#'
#' @description
#'   Display formatted simulation results with operating characteristics including
#'   MTD selection rates, patient allocation, and toxicity outcomes.
#'   Designed for easy interpretation in trial protocol documents.
#'
#' @param summary
#'   List. Summary statistics from `summarize_simulation()`.
#'
#' @param n_doses
#'   Numeric. Number of doses evaluated.
#'
#' @param scenario_name
#'   Character. Optional name for the scenario being presented.
#'   If provided, it will be printed as a header.
#'
#' @return
#'   No return value. Results are printed to console.
#'
#' @details
#'   The function prints four main tables:
#'   1. MTD Selected (%): Percentage of simulations selecting each dose as MTD
#'   2. Number of Participants Treated (mean): Average enrollment at each dose
#'   3. Number of Participants w/ DLTs (mean): Average toxicity count at each dose
#'   4. % No MTD Selected (N/S): Percentage of trials without MTD selection
#'
#'   These operating characteristics are essential for protocol submissions and
#'   to demonstrate design performance across different dose-toxicity scenarios.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' # After running simulation, display results
#' print_simulation_summary_boin(result$summary, n_doses = 3,
#'                         scenario_name = "Scenario 1: Base Case")
#' }
#'
#' @export
print_simulation_summary_boin <- function(summary, n_doses, scenario_name = "") {

  if (scenario_name != "") {
    cat("==================================================\n")
    cat("Scenario:", scenario_name, "\n")
    cat("==================================================\n\n")
  }

  dose_labels <- paste0("DL", seq_len(n_doses))

  # Vectorized rounding operations
  mtd_rounded <- round(summary$mtd_selection_percent, 1)
  n_pts_rounded <- round(summary$avg_n_pts, 1)
  n_tox_rounded <- round(summary$avg_n_tox, 1)

  # (1) MTD selection percentage
  cat("MTD Selected (%)\n")
  mtd_row <- c("", mtd_rounded)
  names(mtd_row) <- c("Item", dose_labels)
  print(mtd_row)
  cat("\n")

  # (2) Average number of patients treated
  cat("Number of Participants Treated (mean)\n")
  n_pts_row <- c("", n_pts_rounded, round(summary$avg_total_n_pts, 1))
  names(n_pts_row) <- c("Item", dose_labels, "Total")
  print(n_pts_row)
  cat("\n")

  # (3) Average number of DLTs
  cat("Number of Participants w/ DLTs (mean)\n")
  n_tox_row <- c("", n_tox_rounded, round(summary$avg_total_n_tox, 1))
  names(n_tox_row) <- c("Item", dose_labels, "Total")
  print(n_tox_row)
  cat("\n")

  # (4) Percentage with no MTD selected
  cat("% No MTD Selected (N/S):", round(summary$percent_no_mtd, 1), "%\n\n")
}
