#' Print Method for BOIN Summary Objects
#'
#' @description
#'   Custom print method for "boin_summary" S3 class objects.
#'   Displays simulation results in a formatted table showing operating
#'   characteristics including true toxicity probabilities, MTD selection rates,
#'   patient allocation, and safety metrics.
#'
#' @param x Object of class "boin_summary" returned by \code{summarize_simulation_boin()}.
#' @param scenario_name Character. Optional name for the scenario being displayed.
#'   If provided, it will be shown as a header. Default is NULL.
#' @param kable_output Logical. If TRUE, output is formatted using \code{knitr::kable()} for
#'   RMarkdown/knitr documents. If FALSE, uses base R formatting. Default is FALSE.
#' @param ... Additional arguments (currently unused, for S3 method consistency).
#'
#' @return Invisibly returns the input object \code{x}. The function is called for its
#'   side effect of printing formatted output.
#'
#' @details
#'   This is an S3 method for the generic \code{print()} function. It is automatically
#'   called when you print or display a "boin_summary" object.
#'
#'   The output includes four key tables:
#'   \enumerate{
#'     \item True Toxicity Probabilities: The true DLT rate at each dose
#'     \item MTD Selected (\%): Percentage of trials selecting each dose as MTD (including No MTD)
#'     \item Number of Participants Treated (mean): Average patient enrollment
#'     \item Number of Participants w/ DLTs (mean): Average DLT counts
#'   }
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
#'
#' result <- sim_boin(
#'   n_trials = 1000,
#'   target = target,
#'   p_true = p_true,
#'   n_doses = 3,
#'   n_cohort = 10,
#'   cohort_size = 3,
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
  p_true <- x$p_true
  mtd_selection_percent <- x$mtd_selection_percent
  avg_n_pts <- x$avg_n_pts
  avg_n_tox <- x$avg_n_tox
  avg_total_n_pts <- x$avg_total_n_pts
  avg_total_n_tox <- x$avg_total_n_tox

  n_doses <- length(p_true)

  # Create column names for doses and "No MTD"
  dose_names <- paste0("DL", 1:n_doses)
  mtd_col_names <- c(dose_names, "No MTD")

  # Print scenario name if provided
  if (!is.null(scenario_name)) {
    cat("\n=== Scenario:", scenario_name, "===\n")
  }

  # ========== Table 1: True Toxicity Probabilities ==========
  cat("True Toxicity Probabilities\n")
  if (kable_output) {
    tox_df <- data.frame(matrix(p_true * 100, nrow = 1))
    colnames(tox_df) <- dose_names
    print(knitr::kable(tox_df, digits = 1, row.names = FALSE))
  } else {
    cat("|")
    for (name in dose_names) {
      cat(sprintf(" %5s |", name))
    }
    cat("\n|")
    for (i in 1:n_doses) {
      cat("------|")
    }
    cat("\n|")
    for (val in p_true * 100) {
      cat(sprintf(" %5.1f |", val))
    }
    cat("\n")
  }
  cat("\n")

  # ========== Table 2: MTD Selection Percentages ==========
  cat("MTD Selected (%)\n")
  if (kable_output) {
    mtd_df <- data.frame(matrix(mtd_selection_percent, nrow = 1))
    colnames(mtd_df) <- mtd_col_names
    print(knitr::kable(mtd_df, digits = 1, row.names = FALSE))
  } else {
    cat("|")
    for (name in mtd_col_names) {
      cat(sprintf(" %6s |", name))
    }
    cat("\n|")
    for (i in 1:(n_doses + 1)) {
      cat("--------|")
    }
    cat("\n|")
    for (val in mtd_selection_percent) {
      cat(sprintf(" %6.1f |", val))
    }
    cat("\n")
  }
  cat("\n")

  # ========== Table 3: Average Number of Participants ==========
  cat("Number of Participants Treated (mean)\n")
  # Add total column
  avg_n_pts_with_total <- c(avg_n_pts, avg_total_n_pts)
  col_names_with_total <- c(dose_names, "Total")

  if (kable_output) {
    n_pts_df <- data.frame(matrix(avg_n_pts_with_total, nrow = 1))
    colnames(n_pts_df) <- col_names_with_total
    print(knitr::kable(n_pts_df, digits = 1, row.names = FALSE))
  } else {
    cat("|")
    for (name in col_names_with_total) {
      cat(sprintf(" %5s |", name))
    }
    cat("\n|")
    for (i in 1:(n_doses + 1)) {
      cat("------|")
    }
    cat("\n|")
    for (val in avg_n_pts_with_total) {
      cat(sprintf(" %5.1f |", val))
    }
    cat("\n")
  }
  cat("\n")

  # ========== Table 4: Average Number of DLTs ==========
  cat("Number of Participants w/ DLTs (mean)\n")
  # Add total column
  avg_n_tox_with_total <- c(avg_n_tox, avg_total_n_tox)

  if (kable_output) {
    n_tox_df <- data.frame(matrix(avg_n_tox_with_total, nrow = 1))
    colnames(n_tox_df) <- col_names_with_total
    print(knitr::kable(n_tox_df, digits = 1, row.names = FALSE))
  } else {
    cat("|")
    for (name in col_names_with_total) {
      cat(sprintf(" %5s |", name))
    }
    cat("\n|")
    for (i in 1:(n_doses + 1)) {
      cat("------|")
    }
    cat("\n|")
    for (val in avg_n_tox_with_total) {
      cat(sprintf(" %5.1f |", val))
    }
    cat("\n")
  }

  # Return invisibly
  invisible(x)
}
