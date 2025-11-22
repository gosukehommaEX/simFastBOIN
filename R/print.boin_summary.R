#' Print Method for BOIN Summary Objects
#'
#' @description
#'   Custom print method for "boin_summary" S3 class objects.
#'   Displays simulation results in a formatted table showing operating
#'   characteristics including MTD selection rates, patient allocation,
#'   and safety metrics.
#'
#' @param x
#'   Object of class "boin_summary" returned by `summarize_simulation_boin()`.
#'
#' @param scenario_name
#'   Character. Optional name for the scenario being displayed.
#'   If provided, it will be shown as a header. Default is NULL.
#'
#' @param kable_output
#'   Logical. If TRUE, output is formatted using `knitr::kable()` for
#'   RMarkdown/knitr documents. If FALSE, uses base R formatting.
#'   Default is FALSE.
#'
#' @param ...
#'   Additional arguments (currently unused, for S3 method consistency).
#'
#' @return
#'   Invisibly returns the input object `x`. The function is called for its
#'   side effect of printing formatted output.
#'
#' @details
#'   This is an S3 method for the generic `print()` function. It is automatically
#'   called when you print or display a "boin_summary" object.
#'
#'   The output includes four key tables:
#'   1. MTD Selected (%): Percentage of trials selecting each dose as MTD
#'   2. Number of Participants Treated (mean): Average patient enrollment
#'   3. Number of Participants w/ DLTs (mean): Average DLT counts
#'   4. % No MTD Selected: Percentage of trials without MTD selection
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' # Run simulation
#' target <- 0.30
#' p_true <- c(0.10, 0.25, 0.40)
#' boin_bound <- get_boin_boundary(target)
#' decision_table <- get_boin_decision(target, boin_bound$lambda_e,
#'                                 boin_bound$lambda_d, 18, 0.95)
#' stopping_boundaries <- get_boin_stopping_boundaries(target, 18, 0.90)
#' result <- sim_boin(
#'   n_trials = 1000,
#'   target = target,
#'   p_true = p_true,
#'   n_doses = 3,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   decision_table = decision_table,
#'   stopping_boundaries = stopping_boundaries,
#'   seed = 123
#' )
#'
#' # Standard output
#' print(result$summary)
#'
#' # Or simply type the object name
#' result$summary
#'
#' # With scenario name
#' print(result$summary, scenario_name = "Base Case")
#'
#' # knitr::kable format for RMarkdown
#' print(result$summary, kable_output = TRUE)
#' }
#'
#' @importFrom knitr kable
#' @export
print.boin_summary <- function(x, scenario_name = NULL, kable_output = FALSE, ...) {

  # Extract summary statistics
  mtd_selection_percent <- x$mtd_selection_percent
  avg_n_pts <- x$avg_n_pts
  avg_n_tox <- x$avg_n_tox
  percent_no_mtd <- x$percent_no_mtd
  avg_total_n_pts <- x$avg_total_n_pts
  avg_total_n_tox <- x$avg_total_n_tox

  n_doses <- length(mtd_selection_percent)

  # Generate dose level labels
  dose_labels <- paste0("DL", 1:n_doses)

  # ===== Display scenario name if provided =====
  if (!is.null(scenario_name)) {
    cat("========================================\n")
    cat("Scenario:", scenario_name, "\n")
    cat("========================================\n\n")
  }

  # ===== Table 1: MTD Selected (%) =====
  if (kable_output) {
    cat("MTD Selected (%)\n\n")
    tbl1 <- data.frame(t(round(mtd_selection_percent, 1)))
    colnames(tbl1) <- dose_labels
    rownames(tbl1) <- NULL
    print(knitr::kable(tbl1, format = "pipe", align = "r"))
  } else {
    cat("MTD Selected (%)\n")
    cat("|", paste(sprintf("%4s", dose_labels), collapse = " | "), "|\n", sep = " ")
    cat("|", paste(rep("------", n_doses), collapse = "|"), "|\n", sep = "")
    cat("|", paste(sprintf("%5.1f", mtd_selection_percent), collapse = " | "), "|\n\n", sep = " ")
  }

  # ===== Table 2: Number of Participants Treated (mean) =====
  if (kable_output) {
    cat("\nNumber of Participants Treated (mean)\n\n")
    tbl2 <- data.frame(t(c(round(avg_n_pts, 1), round(avg_total_n_pts, 1))))
    colnames(tbl2) <- c(dose_labels, "Total")
    rownames(tbl2) <- NULL
    print(knitr::kable(tbl2, format = "pipe", align = "r"))
  } else {
    cat("Number of Participants Treated (mean)\n")
    cat("|", paste(sprintf("%4s", dose_labels), collapse = " | "), "| Total |\n", sep = " ")
    cat("|", paste(rep("------", n_doses + 1), collapse = "|"), "|\n", sep = "")
    cat("|", paste(sprintf("%5.1f", avg_n_pts), collapse = " | "),
        "| ", sprintf("%5.1f", avg_total_n_pts), " |\n\n", sep = " ")
  }

  # ===== Table 3: Number of Participants w/ DLTs (mean) =====
  if (kable_output) {
    cat("\nNumber of Participants w/ DLTs (mean)\n\n")
    tbl3 <- data.frame(t(c(round(avg_n_tox, 1), round(avg_total_n_tox, 1))))
    colnames(tbl3) <- c(dose_labels, "Total")
    rownames(tbl3) <- NULL
    print(knitr::kable(tbl3, format = "pipe", align = "r"))
  } else {
    cat("Number of Participants w/ DLTs (mean)\n")
    cat("|", paste(sprintf("%4s", dose_labels), collapse = " | "), "| Total |\n", sep = " ")
    cat("|", paste(rep("------", n_doses + 1), collapse = "|"), "|\n", sep = "")
    cat("|", paste(sprintf("%5.1f", avg_n_tox), collapse = " | "),
        "| ", sprintf("%5.1f", avg_total_n_tox), " |\n\n", sep = " ")
  }

  # ===== % No MTD Selected =====
  cat("% No MTD Selected (N/S):", sprintf("%.1f", percent_no_mtd), "%\n")

  # Return invisibly for method chaining
  invisible(x)
}
