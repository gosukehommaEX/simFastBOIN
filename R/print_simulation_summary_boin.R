#' Print BOIN Simulation Summary
#'
#' @description
#'   Display formatted simulation results with operating characteristics including
#'   MTD selection rates, patient allocation, and toxicity outcomes.
#'   Designed for easy interpretation in trial protocol documents.
#'
#' @param summary
#'   List. Summary statistics from `summarize_simulation_boin()`.
#'
#' @param scenario_name
#'   Character. Optional name for the scenario being presented.
#'   If provided, it will be printed as a header. Default is "".
#'
#' @param kable_output
#'   Logical. If TRUE, output tables using knitr::kable() format.
#'   If FALSE (default), output as simple formatted vectors.
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
#'   The number of doses is automatically extracted from the summary object,
#'   so explicit specification is not required.
#'
#'   These operating characteristics are essential for protocol submissions and
#'   to demonstrate design performance across different dose-toxicity scenarios.
#'
#' @importFrom knitr kable
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' # After running simulation, display results with default format
#' print_simulation_summary_boin(result$summary,
#'                         scenario_name = "Scenario 1: Base Case")
#'
#' # Display results in knitr::kable format
#' print_simulation_summary_boin(result$summary,
#'                         scenario_name = "Scenario 1: Base Case",
#'                         kable_output = TRUE)
#' }
#'
#' @export
print_simulation_summary_boin <- function(summary, scenario_name = "", kable_output = FALSE) {

  # Automatically extract number of doses from summary
  n_doses <- length(summary$mtd_selection_percent)
  dose_labels <- paste0("DL", seq_len(n_doses))

  if (scenario_name != "") {
    cat("==================================================\n")
    cat("Scenario:", scenario_name, "\n")
    cat("==================================================\n\n")
  }

  # Vectorized rounding operations
  mtd_rounded <- round(summary$mtd_selection_percent, 1)
  n_pts_rounded <- round(summary$avg_n_pts, 1)
  n_tox_rounded <- round(summary$avg_n_tox, 1)

  if (kable_output) {
    # Using knitr::kable format

    # (1) MTD selection percentage
    cat("MTD Selected (%)\n")
    mtd_df <- data.frame(matrix(mtd_rounded, nrow = 1, byrow = TRUE))
    colnames(mtd_df) <- dose_labels
    print(knitr::kable(mtd_df, format = "simple"))
    cat("\n")

    # (2) Average number of patients treated
    cat("Number of Participants Treated (mean)\n")
    n_pts_df <- data.frame(matrix(c(n_pts_rounded, round(summary$avg_total_n_pts, 1)), nrow = 1, byrow = TRUE))
    colnames(n_pts_df) <- c(dose_labels, "Total")
    print(knitr::kable(n_pts_df, format = "simple"))
    cat("\n")

    # (3) Average number of DLTs
    cat("Number of Participants w/ DLTs (mean)\n")
    n_tox_df <- data.frame(matrix(c(n_tox_rounded, round(summary$avg_total_n_tox, 1)), nrow = 1, byrow = TRUE))
    colnames(n_tox_df) <- c(dose_labels, "Total")
    print(knitr::kable(n_tox_df, format = "simple"))
    cat("\n")

  } else {
    # Original simple vector format

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
  }

  # (4) Percentage with no MTD selected
  cat("% No MTD Selected (N/S):", round(summary$percent_no_mtd, 1), "%\n\n")
}
