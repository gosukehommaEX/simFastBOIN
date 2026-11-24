#' Print Operating Characteristics of the Waterfall Design
#'
#' @description
#'   Print the design, the true DLT rates, the selection percentages and the
#'   average numbers of patients and DLTs of every combination, and the
#'   summaries of the MTD contour, for a \code{waterfall_oc} object.
#'
#' @param x
#'   An object of class \code{waterfall_oc}, as returned by
#'   \code{\link{sim_waterfall}}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @param kable
#'   Logical scalar. Print the table with \code{knitr::kable()} instead of the
#'   plain console layout. Defaults to \code{FALSE}.
#'
#' @param kable_format
#'   Character scalar. Format passed to \code{knitr::kable()}.
#'   Defaults to \code{"pipe"}.
#'
#' @param ...
#'   Not used.
#'
#' @return
#'   \code{x}, invisibly.
#'
#' @examples
#' p_true <- matrix(c(0.03, 0.10, 0.28,
#'                    0.10, 0.30, 0.50), nrow = 2, byrow = TRUE)
#'
#' oc <- sim_waterfall(
#'   target = 0.30,
#'   p_true = p_true,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#'
#' print(oc)
#' print(oc, digits = 2)
#'
#' @export
print.waterfall_oc <- function(x, digits = 1, kable = FALSE,
                               kable_format = "pipe", ...) {

  tab <- comb_oc_table(x, digits)

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    print(knitr::kable(tab, format = kable_format, row.names = FALSE))
    return(invisible(x))
  }

  set <- x$settings
  cat("Waterfall drug-combination operating characteristics\n")
  cat("  target DLT rate : ", format(x$target), "\n", sep = "")
  cat("  dose levels     : ", x$n_rows, " of drug A, ", x$n_cols,
      " of drug B\n", sep = "")
  cat("  trials          : ", x$n_trials, "\n", sep = "")
  cat("  cohorts         : at most ", paste(set$n_cohort, collapse = "/"),
      " per subtrial, of size ", set$cohort_size, "\n", sep = "")
  cat("  max sample size : ", set$max_total_pts, "\n", sep = "")
  cat("  lambda_e / _d   : ", format(round(x$lambda_e, 4)), " / ",
      format(round(x$lambda_d, 4)), "\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  cat("\nNo MTD selected at any level of drug A : ",
      format(round(x$percent_no_mtd, digits)), "%\n", sep = "")
  cat("Average number of patients             : ",
      format(round(x$total_n_pts, digits)), "\n", sep = "")
  cat("Average number of DLTs                 : ",
      format(round(x$total_n_tox, digits)), "\n", sep = "")
  cat("Average number of subtrials            : ",
      format(round(x$avg_n_subtrials, digits)), "\n", sep = "")

  contour <- if (nrow(x$true_mtd) > 0L) {
    paste0("A", x$true_mtd[, "dose_a"], "B", x$true_mtd[, "dose_b"],
           collapse = ", ")
  } else {
    "none"
  }
  cat("\nTrue MTD contour                       : ", contour, "\n", sep = "")
  cat("  Trials selecting exactly the contour : ",
      format(round(x$pct_contour_selected, digits)), "%\n", sep = "")
  cat("  Patients treated on the contour      : ",
      format(round(x$pct_patients_contour, digits)), "%\n", sep = "")
  cat("  Patients treated above the contour   : ",
      format(round(x$pct_patients_above, digits)), "%\n", sep = "")
  cat("  Patients treated below the contour   : ",
      format(round(x$pct_patients_below, digits)), "%\n", sep = "")

  invisible(x)
}
