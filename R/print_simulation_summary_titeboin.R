#' Print TITE-BOIN Simulation Summary
#'
#' @description
#' Display a formatted summary of TITE-BOIN simulation results, including MTD
#' selection rates, patient allocation, toxicity counts, and trial duration statistics.
#'
#' @param summary Summary list from \code{sim_titeboin()}, containing:
#'   \itemize{
#'     \item \code{mtd_selection_percent}: Percentage of trials selecting each dose as MTD
#'     \item \code{avg_n_pts}: Average number of patients treated at each dose
#'     \item \code{avg_total_n_pts}: Average total number of patients per trial
#'     \item \code{avg_n_tox}: Average number of toxicities at each dose
#'     \item \code{avg_total_n_tox}: Average total toxicities per trial
#'     \item \code{percent_no_mtd}: Percentage of trials with no MTD selected
#'     \item \code{avg_duration}: Average trial duration
#'     \item \code{sd_duration}: Standard deviation of trial duration
#'   }
#' @param n_doses Number of dose levels in the trial
#' @param scenario_name Optional character string describing the simulation scenario.
#'   If provided, it will be displayed as a header (default: "")
#'
#' @return NULL (prints formatted output to console)
#'
#' @details
#' The function displays four main sections:
#' \enumerate{
#'   \item MTD Selection Rate: Percentage of trials selecting each dose as MTD
#'   \item Average Number of Patients Treated: Mean patient allocation by dose
#'   \item Average Number of DLTs: Mean toxicity counts by dose
#'   \item Summary Statistics: MTD not selected rate and trial duration
#' }
#'
#' All percentages and counts are rounded to one decimal place for readability.
#'
#' @references
#' Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K.E. (2018).
#' Time-to-Event Bayesian Optimal Interval Design to Accelerate Phase I Trials.
#' Clinical Cancer Research, 24(20): 4921-4930.
#'
#' @examples
#' \dontrun{
#' # After running sim_titeboin()
#' result <- sim_titeboin(
#'   n_trials = 1000,
#'   decision_table = decision_table,
#'   stopping_boundaries = stopping_boundaries,
#'   target = 0.30,
#'   p_true = c(0.25, 0.41, 0.45, 0.49, 0.53),
#'   n_doses = 5,
#'   ncohort = 10,
#'   cohort_size = 3,
#'   maxt = 3,
#'   accrual_rate = 1
#' )
#'
#' # Display results
#' print_simulation_summary_titeboin(
#'   result$summary,
#'   n_doses = 5,
#'   scenario_name = "TITE-BOIN (5 doses, 30 patients)"
#' )
#' }
#'
#' @seealso
#' \code{\link{sim_titeboin}} for running TITE-BOIN simulations
#'
#' @export
print_simulation_summary_titeboin <- function(summary, n_doses, scenario_name = "") {
  
  if (scenario_name != "") {
    cat("==================================================\n")
    cat("Scenario:", scenario_name, "\n")
    cat("==================================================\n\n")
  }
  
  dose_labels <- paste("DL", 1:n_doses, sep = "")
  
  # MTD Selection Rate
  cat("MTD Selection Rate (%)\n")
  mtd_row <- c("", round(summary$mtd_selection_percent, 1))
  names(mtd_row) <- c("Item", dose_labels)
  print(mtd_row)
  cat("\n")
  
  # Average Number of Patients Treated
  cat("Average Number of Patients Treated\n")
  n_pts_row <- c("", round(summary$avg_n_pts, 1), round(summary$avg_total_n_pts, 1))
  names(n_pts_row) <- c("Item", dose_labels, "Total")
  print(n_pts_row)
  cat("\n")
  
  # Average Number of DLTs
  cat("Average Number of DLTs\n")
  n_tox_row <- c("", round(summary$avg_n_tox, 1), round(summary$avg_total_n_tox, 1))
  names(n_tox_row) <- c("Item", dose_labels, "Total")
  print(n_tox_row)
  cat("\n")
  
  # Summary statistics
  cat("MTD Not Selected Rate (%):", round(summary$percent_no_mtd, 1), "%\n")
  cat("Average Trial Duration:", round(summary$avg_duration, 1),
      "(SD:", round(summary$sd_duration, 1), ")\n\n")
}