#' Print Method for BOIN Multi-Scenario Summary
#'
#' @description
#'   Custom print method for "boin_multi_summary" objects returned by
#'   \code{sim_boin_multi()}. Displays aggregated operating characteristics
#'   in formatted tables with options for plain text or kable output.
#'
#' @usage
#'   \method{print}{boin_multi_summary}(x, scenario_name = NULL, percent = FALSE, kable = FALSE, kable_format = "pipe", ...)
#'
#' @param x
#'   Object of class "boin_multi_summary" returned by \code{sim_boin_multi()}.
#'
#' @param scenario_name
#'   Character. Optional name for display header (for compatibility with print.boin_summary).
#'   Ignored for multi-scenario results. Default is NULL.
#'
#' @param percent
#'   Logical. If TRUE, display "Avg Pts" and "Avg DLTs" as percentages of totals.
#'   If FALSE (default), display as absolute numbers. Default is FALSE.
#'
#' @param kable
#'   Logical. If TRUE, format output as knitr::kable table. If FALSE (default),
#'   display as plain text table. Default is FALSE.
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
#'   Displays a summary table with:
#'   \itemize{
#'     \item Rows organized by scenario and metric (True Tox (%), MTD Sel (%), Avg Pts, Avg DLTs)
#'     \item Columns for each dose level plus Total/No MTD statistics
#'     \item Scenario names grouped for easy visual comparison
#'     \item All text left-aligned for consistency
#'   }
#'
#'   When \code{percent = FALSE} (default), "Avg Pts" and "Avg DLTs" are displayed
#'   as absolute numbers.
#'
#'   When \code{percent = TRUE}, these rows display percentages of their respective
#'   totals for each trial. This facilitates visual comparison of dose allocation
#'   patterns across scenarios independent of total sample size.
#'
#'   When \code{kable = FALSE} (default), displays as plain text table in console.
#'
#'   When \code{kable = TRUE}, formats the table using \code{knitr::kable()} with
#'   the specified format. This is useful for R Markdown documents and reports.
#'   The \code{kable_format} parameter controls output format:
#'   \itemize{
#'     \item \code{"pipe"}: Markdown pipe format (best for GitHub/Markdown)
#'     \item \code{"simple"}: Minimal formatting (plain text)
#'     \item \code{"latex"}: LaTeX format (for PDF reports)
#'     \item \code{"html"}: HTML format (for web display)
#'   }
#'
#'   All columns are left-aligned for improved readability across all formats.
#'
#' @examples
#' \dontrun{
#' # Create multi-scenario simulation results
#' scenarios <- list(
#'   list(name = "Scenario 1",
#'        p_true = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)),
#'   list(name = "Scenario 2",
#'        p_true = c(0.01, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.30))
#' )
#'
#' result <- sim_boin_multi(
#'   scenarios = scenarios,
#'   target = 0.30,
#'   n_trials = 1000,
#'   n_cohort = 48,
#'   cohort_size = 3,
#'   seed = 123
#' )
#'
#' # Print as plain text with absolute numbers (default)
#' print(result)
#'
#' # Print with absolute numbers as Markdown table
#' print(result, kable = TRUE, kable_format = "pipe")
#'
#' # Print with percentages instead of absolute numbers
#' print(result, percent = TRUE)
#'
#' # Print with percentages as Markdown table
#' print(result, kable = TRUE, kable_format = "pipe", percent = TRUE)
#'
#' # Print as LaTeX table for PDF documents
#' print(result, kable = TRUE, kable_format = "latex")
#'
#' # Print as HTML table
#' print(result, kable = TRUE, kable_format = "html")
#' }
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling add_header_above row_spec
#'
#' @export
print.boin_multi_summary <- function(x, scenario_name = NULL, percent = FALSE, kable = FALSE, kable_format = "pipe", ...) {

  # Validate kable_format argument
  valid_formats <- c("pipe", "simple", "latex", "html")
  kable_format <- match.arg(kable_format, valid_formats)

  # Extract data frame and create working copy
  df <- x$combined_summary_df
  df_display <- df

  # Process percent conversion if requested
  if (percent) {
    n_doses <- x$n_doses
    n_scenarios <- length(x$scenario_names)

    # Iterate through each scenario (4 rows per scenario)
    for (i in seq_len(n_scenarios)) {
      # Rows for this scenario
      row_start <- (i - 1) * 4 + 1
      row_pts <- row_start + 2      # "Number of Participants Treated" row
      row_tox <- row_start + 3      # "Number of Participants w/ DLTs" row

      # Convert "Number of Participants Treated" to percentages
      if (df_display[row_pts, 2] == "Number of Participants Treated") {
        # Extract numeric values (columns 3 to n_doses+2, and total in n_doses+3)
        pts_values <- as.numeric(df_display[row_pts, 3:(n_doses + 2)])
        total_pts <- as.numeric(df_display[row_pts, n_doses + 3])

        # Calculate percentages
        pct_values <- (pts_values / total_pts) * 100
        df_display[row_pts, 3:(n_doses + 2)] <- as.character(round(pct_values, 1))
      }

      # Convert "Number of Participants w/ DLTs" to percentages
      if (df_display[row_tox, 2] == "Number of Participants w/ DLTs") {
        # Extract numeric values (columns 3 to n_doses+2, and total in n_doses+3)
        tox_values <- as.numeric(df_display[row_tox, 3:(n_doses + 2)])
        total_tox <- as.numeric(df_display[row_tox, n_doses + 3])

        # Calculate percentages
        pct_values <- (tox_values / total_tox) * 100
        df_display[row_tox, 3:(n_doses + 2)] <- as.character(round(pct_values, 1))
      }
    }
  }

  if (!kable) {
    # ===== PLAIN TEXT FORMAT OUTPUT =====
    # Standardize metric names for display
    df_display[, 2] <- sapply(df_display[, 2], function(item) {
      if (item == "True DLT Rate (%)") return("True Tox (%)")
      if (item == "MTD Selected (%)") return("MTD Sel (%)")
      if (item == "Number of Participants Treated") {
        if (percent) return("Avg Pts (%)") else return("Avg Pts")
      }
      if (item == "Number of Participants w/ DLTs") {
        if (percent) return("Avg DLTs (%)") else return("Avg DLTs")
      }
      return(item)
    })

    # Print as plain text table
    print(df_display, row.names = FALSE)
  } else {
    # ===== KABLE FORMAT OUTPUT =====
    # Create working copy for kable output
    df_kable <- df_display

    # Standardize metric names for display
    df_kable[, 2] <- sapply(df_kable[, 2], function(item) {
      if (item == "True DLT Rate (%)") return("True Tox (%)")
      if (item == "MTD Selected (%)") return("MTD Sel (%)")
      if (item == "Number of Participants Treated") {
        if (percent) return("Avg Pts (%)") else return("Avg Pts")
      }
      if (item == "Number of Participants w/ DLTs") {
        if (percent) return("Avg DLTs (%)") else return("Avg DLTs")
      }
      return(item)
    })

    # Rename columns for display (replace underscores with spaces)
    colnames(df_kable) <- gsub("_", " ", colnames(df_kable))

    # Create kable table with all columns left-aligned
    table_output <- knitr::kable(
      df_kable,
      format = kable_format,
      escape = FALSE,
      align = rep("l", ncol(df_kable))
    )

    # Apply additional styling for HTML format
    if (kable_format == "html" && requireNamespace("kableExtra", quietly = TRUE)) {
      # Apply base styling
      table_output <- kableExtra::kable_styling(
        table_output,
        bootstrap_options = c("striped", "hover", "condensed", "responsive"),
        full_width = FALSE,
        position = "center"
      )

      # Add header above dose columns with thick border below
      n_doses <- x$n_doses
      table_output <- kableExtra::add_header_above(
        table_output,
        c(" " = 2, "Operating Characteristics" = n_doses + 1),
        bold = TRUE,
        extra_css = "border-bottom: 2px solid #000;"
      )

      # Make column header (Scenario, Item, DL1, ..., Total/No MTD) bold with thick border below
      table_output <- kableExtra::row_spec(
        table_output,
        0,
        bold = TRUE,
        extra_css = "border-bottom: 2px solid #000;"
      )

      # Add bold borders between scenarios (at the end of each scenario's 4th row)
      n_scenarios <- length(x$scenario_names)
      n_data_rows <- nrow(df_kable)

      for (i in seq_len(n_scenarios)) {
        # Last row of each scenario (every 4th row)
        last_row_of_scenario <- i * 4
        if (last_row_of_scenario <= n_data_rows) {
          table_output <- kableExtra::row_spec(
            table_output,
            last_row_of_scenario,
            extra_css = "border-bottom: 2px solid #000;"
          )
        }
      }
    }

    print(table_output)
  }

  invisible(x)
}
