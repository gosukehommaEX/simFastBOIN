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
#'   The output is presented as a unified table with four rows:
#'   \enumerate{
#'     \item True Toxicity (\%): The true DLT rate at each dose
#'     \item MTD Selected (\%): Percentage of trials selecting each dose as MTD
#'     \item Participants Treated (mean): Average patient enrollment at each dose
#'     \item Participants w/ DLTs (mean): Average DLT counts at each dose
#'   }
#'
#'   The last column shows "No MTD" for MTD selection rate, and "Total" for
#'   participant counts and DLT counts.
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

  # Create column names for doses and last column
  dose_names <- paste0("DL", 1:n_doses)

  # Print scenario name if provided
  if (!is.null(scenario_name)) {
    cat("\n=== Scenario:", scenario_name, "===\n\n")
  }

  # ========== Unified Table Format ==========
  if (kable_output) {
    # Create unified data frame
    unified_df <- data.frame(
      Metric = c(
        "True Toxicity (%)",
        "MTD Selected (%)",
        "Participants Treated (mean)",
        "Participants w/ DLTs (mean)"
      ),
      stringsAsFactors = FALSE
    )

    # Add dose columns
    for (i in 1:n_doses) {
      unified_df[[dose_names[i]]] <- c(
        p_true[i] * 100,
        mtd_selection_percent[i],
        avg_n_pts[i],
        avg_n_tox[i]
      )
    }

    # Add last column (No MTD for row 2, Total for rows 3-4, NA for row 1)
    unified_df[["Total/No MTD"]] <- c(
      NA,  # True toxicity has no total
      mtd_selection_percent[n_doses + 1],  # No MTD percentage
      avg_total_n_pts,
      avg_total_n_tox
    )

    # Print using knitr::kable
    cat("BOIN Simulation Summary\n\n")
    print(knitr::kable(unified_df, digits = 1, row.names = FALSE, align = c("l", rep("r", n_doses + 1))))

  } else {
    # Base R format - create unified table
    cat("BOIN Simulation Summary\n")

    # Header row
    cat(sprintf("| %-32s |", "Metric"))
    for (name in dose_names) {
      cat(sprintf(" %6s |", name))
    }
    cat(" Total/No MTD |\n")

    # Separator row
    cat(sprintf("| %s |", paste(rep("-", 32), collapse = "")))
    for (i in 1:n_doses) {
      cat("--------|")
    }
    cat("--------------|\n")

    # Row 1: True Toxicity
    cat(sprintf("| %-32s |", "True Toxicity (%)"))
    for (val in p_true * 100) {
      cat(sprintf(" %6.1f |", val))
    }
    cat(sprintf(" %12s |\n", "-"))

    # Row 2: MTD Selected
    cat(sprintf("| %-32s |", "MTD Selected (%)"))
    for (i in 1:n_doses) {
      cat(sprintf(" %6.1f |", mtd_selection_percent[i]))
    }
    cat(sprintf(" %12.1f |\n", mtd_selection_percent[n_doses + 1]))

    # Row 3: Participants Treated
    cat(sprintf("| %-32s |", "Participants Treated (mean)"))
    for (val in avg_n_pts) {
      cat(sprintf(" %6.1f |", val))
    }
    cat(sprintf(" %12.1f |\n", avg_total_n_pts))

    # Row 4: Participants w/ DLTs
    cat(sprintf("| %-32s |", "Participants w/ DLTs (mean)"))
    for (val in avg_n_tox) {
      cat(sprintf(" %6.1f |", val))
    }
    cat(sprintf(" %12.1f |\n", avg_total_n_tox))
  }

  # Return invisibly
  invisible(x)
}
