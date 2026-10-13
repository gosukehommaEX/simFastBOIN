#' Print Operating Characteristics of the BOIN Drug-Combination Design
#'
#' @description
#'   Display the summary produced by \code{\link{sim_comb_boin}}: the true DLT
#'   rates, the selection percentages and the average numbers of patients and
#'   DLTs at every dose combination, followed by the selection of the true MTD
#'   and the exposure to combinations above the overdose cutoff.
#'
#' @param x
#'   An object of class \code{comb_boin_oc}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @param kable
#'   Logical scalar. Return the table through \code{knitr::kable()} rather than
#'   printing it directly, which is convenient inside a report. Requires the
#'   \pkg{knitr} package. Defaults to \code{FALSE}.
#'
#' @param kable_format
#'   Character scalar passed to \code{knitr::kable()}, for example \code{"pipe"},
#'   \code{"html"} or \code{"latex"}. Defaults to \code{"pipe"}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' p_true <- matrix(c(0.02, 0.04, 0.08, 0.14,
#'                    0.08, 0.25, 0.42, 0.48,
#'                    0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE)
#'
#' oc <- sim_comb_boin(
#'   target = 0.25,
#'   p_true = p_true,
#'   n_cohort = 16,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#'
#' print(oc)
#' print(oc, digits = 2)
#'
#' @export
print.comb_boin_oc <- function(x, digits = 1, kable = FALSE,
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
  cat("BOIN drug-combination operating characteristics\n")
  cat("  target DLT rate : ", format(x$target), "\n", sep = "")
  cat("  dose levels     : ", x$n_rows, " of drug A, ", x$n_cols,
      " of drug B\n", sep = "")
  cat("  trials          : ", x$n_trials, "\n", sep = "")
  cat("  cohorts         : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size : ", set$max_total_pts, "\n", sep = "")
  cat("  lambda_e / _d   : ", format(round(x$lambda_e, 4)), " / ",
      format(round(x$lambda_d, 4)), "\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  cat("\nNo MTD selected                      : ",
      format(round(x$percent_no_mtd, digits)), "%\n", sep = "")
  cat("Average number of patients           : ",
      format(round(x$total_n_pts, digits)), "\n", sep = "")
  cat("Average number of DLTs               : ",
      format(round(x$total_n_tox, digits)), "\n", sep = "")

  cat("\nTrue MTD (closest to the target)     : ",
      paste0("A", x$true_mtd[, "dose_a"], "B", x$true_mtd[, "dose_b"],
             collapse = ", "), "\n", sep = "")
  cat("  Selected as the MTD                : ",
      format(round(x$pct_true_mtd_selected, digits)), "%\n", sep = "")
  cat("  Patients treated there             : ",
      format(round(x$pct_patients_true_mtd, digits)), "%\n", sep = "")

  overdose <- x$overdose
  cat("\nCombinations above a true DLT rate of ", format(overdose$cutoff),
      ": ", nrow(overdose$combinations), " of ", x$n_rows * x$n_cols, "\n",
      sep = "")
  cat("  Patients treated there             : ",
      format(round(overdose$pct_patients, digits)), "%\n", sep = "")
  cat("  Trials treating any patient there  : ",
      format(round(overdose$pct_trials_any, digits)), "%\n", sep = "")
  cat("  Trials treating over 60% there     : ",
      format(round(overdose$pct_trials_over_60, digits)), "%\n", sep = "")
  cat("  Trials treating over 80% there     : ",
      format(round(overdose$pct_trials_over_80, digits)), "%\n", sep = "")
  cat("  Trials selecting an MTD there      : ",
      format(round(overdose$pct_trials_mtd_above, digits)), "%\n", sep = "")

  invisible(x)
}
