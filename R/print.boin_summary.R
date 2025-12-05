#' Print Method for BOIN Summary Objects
#'
#' @description
#'   Custom print method for "boin_summary" S3 class objects returned by
#'   \code{sim_boin()}. Displays operating characteristics in a formatted table.
#'
#' @usage
#'   \method{print}{boin_summary}(x, scenario_name = NULL, percent = FALSE, kable = FALSE, kable_format = "pipe", ...)
#'
#' @param x
#'   Object of class "boin_summary" returned by \code{sim_boin()}.
#'
#' @param scenario_name
#'   Character. Optional name for the scenario. If provided, displayed as header.
#'   Default is NULL.
#'
#' @param percent
#'   Logical. If TRUE, display average patients and DLTs as percentages of totals.
#'   If FALSE (default), display as absolute numbers.
#'
#' @param kable
#'   Logical. If TRUE, format output as knitr::kable table. If FALSE (default),
#'   display as plain text table.
#'
#' @param kable_format
#'   Character. Format for kable output when kable = TRUE. Options include
#'   "pipe" (Markdown pipes, default), "simple" (minimal formatting),
#'   "latex" (LaTeX format), and "html" (HTML format). Default is "pipe".
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
#'     \item True Tox (%): True DLT rate at each dose
#'     \item MTD Sel (%): Percentage selecting each dose as MTD (+ No MTD)
#'     \item Avg Pts: Average enrollment per dose (absolute or percentage)
#'     \item Avg DLTs: Average DLT count per dose (absolute or percentage)
#'   }
#'
#'   When \code{percent = TRUE}, Avg Pts and Avg DLTs are displayed as
#'   percentages of their respective totals for each trial.
#'
#'   When \code{kable = TRUE}, output is formatted using \code{knitr::kable()}.
#'   This is useful for R Markdown documents and reports. The \code{kable_format}
#'   parameter controls the output format: "pipe" for Markdown pipes (default),
#'   "simple" for minimal formatting, "latex" for LaTeX tables, and "html" for
#'   HTML tables with enhanced styling.
#'
#'   When \code{kable_format = "html"}, additional styling is applied including
#'   striped rows, hover effects, and responsive formatting via kableExtra.
#'
#' @examples
#' \dontrun{
#' # Create BOIN simulation results
#' result <- sim_boin(
#'   n_trials = 1000,
#'   target = 0.30,
#'   p_true = c(0.10, 0.25, 0.40),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   seed = 123
#' )
#'
#' # Print with absolute numbers (default)
#' print(result$summary, scenario_name = "Absolute Numbers")
#'
#' # Print with percentages
#' print(result$summary, scenario_name = "Percentages", percent = TRUE)
#'
#' # Print as Markdown table
#' print(result$summary, kable = TRUE, kable_format = "pipe")
#'
#' # Print as HTML table with enhanced styling
#' print(result$summary, kable = TRUE, kable_format = "html")
#' }
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @importFrom knitr kable
#'
#' @export
print.boin_summary <- function(x, scenario_name = NULL, percent = FALSE, kable = FALSE, kable_format = "pipe", ...) {

  # Validate kable_format argument
  valid_formats <- c("pipe", "simple", "latex", "html")
  kable_format <- match.arg(kable_format, valid_formats)

  # Extract components from the boin_summary object
  p_true <- x$p_true
  mtd_selection_percent <- x$mtd_selection_percent
  avg_n_pts <- x$avg_n_pts
  avg_n_tox <- x$avg_n_tox
  avg_total_n_pts <- x$avg_total_n_pts
  avg_total_n_tox <- x$avg_total_n_tox

  # Convert to percentages if requested
  if (percent) {
    avg_n_pts_display <- (avg_n_pts / avg_total_n_pts) * 100
    avg_n_tox_display <- (avg_n_tox / avg_total_n_tox) * 100
  } else {
    avg_n_pts_display <- avg_n_pts
    avg_n_tox_display <- avg_n_tox
  }

  # Get number of doses
  n_doses <- length(p_true)

  # Create dose labels (DL1, DL2, ..., DLn)
  dose_labels <- paste0("DL", 1:n_doses)

  # ===== KABLE FORMAT OUTPUT =====
  if (kable) {
    # Create row labels
    row_labels <- c(
      "True Tox (%)",
      "MTD Sel (%)",
      if (percent) "Avg Pts (%)" else "Avg Pts",
      if (percent) "Avg DLTs (%)" else "Avg DLTs"
    )

    # Build data frame for kable
    table_data <- data.frame(
      Metric = row_labels,
      matrix(c(
        round(p_true * 100, 1),
        round(mtd_selection_percent[1:n_doses], 1),
        round(avg_n_pts_display, 1),
        round(avg_n_tox_display, 1),
        NA,
        round(mtd_selection_percent[n_doses + 1], 1),
        round(avg_total_n_pts, 1),
        round(avg_total_n_tox, 1)
      ), nrow = 4, byrow = TRUE),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    # Set column names
    colnames(table_data) <- c("Metric", dose_labels, "Total/No MTD")

    # Convert to character for proper display
    for (i in 2:ncol(table_data)) {
      table_data[, i] <- as.character(table_data[, i])
    }

    # Create kable table with all columns left-aligned
    table_output <- knitr::kable(
      table_data,
      format = kable_format,
      escape = FALSE,
      align = rep("l", ncol(table_data))
    )

    # Apply additional styling for HTML format
    if (kable_format == "html" && requireNamespace("kableExtra", quietly = TRUE)) {
      table_output <- kableExtra::kable_styling(
        table_output,
        bootstrap_options = c("striped", "hover", "condensed", "responsive"),
        full_width = FALSE,
        position = "center"
      )

      # Add header above dose columns
      table_output <- kableExtra::add_header_above(
        table_output,
        c(" " = 1, "Operating Characteristics" = n_doses, " " = 1)
      )

      # Add scenario name header if provided
      if (!is.null(scenario_name)) {
        table_output <- kableExtra::add_header_above(
          table_output,
          c(" " = 1, scenario_name = ncol(table_data) - 1),
          bold = TRUE
        )
      }
    }

    print(table_output)
  } else {
    # ===== PLAIN TEXT FORMAT OUTPUT =====

    # Print scenario name if provided
    if (!is.null(scenario_name)) {
      cat("Scenario: ", scenario_name, "\n\n", sep = "")
    }

    # Set column width for formatting
    col_width <- 10

    # Print table header with dose labels and total/no MTD column
    cat(sprintf("|%*s", 15, "Metric"))
    for (i in 1:n_doses) {
      cat(sprintf("|%*s", col_width, dose_labels[i]))
    }
    cat(sprintf("|%*s|\n", col_width, "Total/No MTD"))

    # Print horizontal separator line
    separator_width <- 15 + (n_doses + 1) * col_width + (n_doses + 2)
    cat(paste(rep("-", separator_width), collapse = ""), "\n")

    # Row 1: True toxicity probability (%)
    cat(sprintf("|%*s", 15, "True Tox (%)"))
    for (i in 1:n_doses) {
      cat(sprintf("|%*.1f", col_width, p_true[i] * 100))
    }
    cat(sprintf("|%*s|\n", col_width, ""))

    # Row 2: MTD selection percentage for each dose and no MTD
    cat(sprintf("|%*s", 15, "MTD Sel (%)"))
    for (i in 1:n_doses) {
      cat(sprintf("|%*.1f", col_width, mtd_selection_percent[i]))
    }
    cat(sprintf("|%*.1f|\n", col_width, mtd_selection_percent[n_doses + 1]))

    # Row 3: Average number of patients per dose and total
    if (percent) {
      cat(sprintf("|%*s", 15, "Avg Pts (%)"))
    } else {
      cat(sprintf("|%*s", 15, "Avg Pts"))
    }
    for (i in 1:n_doses) {
      cat(sprintf("|%*.1f", col_width, avg_n_pts_display[i]))
    }
    cat(sprintf("|%*.1f|\n", col_width, avg_total_n_pts))

    # Row 4: Average number of DLTs per dose and total
    if (percent) {
      cat(sprintf("|%*s", 15, "Avg DLTs (%)"))
    } else {
      cat(sprintf("|%*s", 15, "Avg DLTs"))
    }
    for (i in 1:n_doses) {
      cat(sprintf("|%*.1f", col_width, avg_n_tox_display[i]))
    }
    cat(sprintf("|%*.1f|\n", col_width, avg_total_n_tox))

    # Print horizontal separator line at bottom
    cat(paste(rep("-", separator_width), collapse = ""), "\n")
  }

  # Return input object invisibly (for chaining/consistency)
  invisible(x)
}
