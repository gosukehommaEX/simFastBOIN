#' Print Method for BOIN Summary Objects
#'
#' @description
#'   Custom print method for "boin_summary" S3 class objects returned by
#'   \code{sim_boin()}. Displays operating characteristics in a formatted table.
#'
#' @usage
#'   \method{print}{boin_summary}(x, scenario_name = NULL, ...)
#'
#' @param x
#'   Object of class "boin_summary" returned by \code{sim_boin()}.
#'
#' @param scenario_name
#'   Character. Optional name for the scenario. If provided, displayed as header.
#'   Default is NULL.
#'
#' @param ...
#'   Additional arguments (currently unused, for S3 method consistency).
#'
#' @return
#'   Invisibly returns the input object \code{x}. Function is called for printing side effect.
#'
#' @details
#'   Displays operating characteristics in a unified table with rows:
#'   \enumerate{
#'     \item True Toxicity (%): True DLT rate at each dose
#'     \item MTD Selected (%): Percentage selecting each dose as MTD (+ No MTD)
#'     \item Avg Patients: Average enrollment per dose (+ Total)
#'     \item Avg DLTs: Average DLT count per dose (+ Total)
#'   }
#'
#' @examples
#' result <- sim_boin(
#'   n_trials = 1000,
#'   target = 0.30,
#'   p_true = c(0.10, 0.25, 0.40),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   seed = 123
#' )
#'
#' # Automatic print on display
#' result
#'
#' # Or explicit print with scenario name
#' print(result, scenario_name = "Base Case")
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @export
print.boin_summary <- function(x, scenario_name = NULL, ...) {

  p_true <- x$p_true
  mtd_selection_percent <- x$mtd_selection_percent
  avg_n_pts <- x$avg_n_pts
  avg_n_tox <- x$avg_n_tox
  avg_total_n_pts <- x$avg_total_n_pts
  avg_total_n_tox <- x$avg_total_n_tox

  n_doses <- length(p_true)

  # Print scenario name if provided
  if (!is.null(scenario_name)) {
    cat("Scenario:", scenario_name, "\n\n")
  }

  # Create dose labels
  dose_labels <- paste0("DL", 1:n_doses)

  # Define column widths
  col_width <- 10

  # Print table header
  cat(sprintf("|%*s", 15, "Metric"))
  for (i in 1:n_doses) {
    cat(sprintf("|%*s", col_width, dose_labels[i]))
  }
  cat(sprintf("|%*s|\n", col_width, "Total/No MTD"))

  # Print separator
  separator_width <- 15 + (n_doses + 1) * col_width + (n_doses + 2)
  cat(paste(rep("-", separator_width), collapse = ""), "\n")

  # Row 1: True Toxicity
  cat(sprintf("|%*s", 15, "True Tox (%)"))
  for (i in 1:n_doses) {
    cat(sprintf("|%*.1f", col_width, p_true[i] * 100))
  }
  cat(sprintf("|%*s|\n", col_width, ""))

  # Row 2: MTD Selected
  cat(sprintf("|%*s", 15, "MTD Sel (%)"))
  for (i in 1:n_doses) {
    cat(sprintf("|%*.1f", col_width, mtd_selection_percent[i]))
  }
  cat(sprintf("|%*.1f|\n", col_width, mtd_selection_percent[n_doses + 1]))

  # Row 3: Avg Patients
  cat(sprintf("|%*s", 15, "Avg Pts"))
  for (i in 1:n_doses) {
    cat(sprintf("|%*.1f", col_width, avg_n_pts[i]))
  }
  cat(sprintf("|%*.1f|\n", col_width, avg_total_n_pts))

  # Row 4: Avg DLTs
  cat(sprintf("|%*s", 15, "Avg DLTs"))
  for (i in 1:n_doses) {
    cat(sprintf("|%*.1f", col_width, avg_n_tox[i]))
  }
  cat(sprintf("|%*.1f|\n", col_width, avg_total_n_tox))

  # Print bottom separator
  cat(paste(rep("-", separator_width), collapse = ""), "\n")

  invisible(x)
}
