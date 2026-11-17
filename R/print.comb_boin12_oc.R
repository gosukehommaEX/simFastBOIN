#' Print the Operating Characteristics of the Comb-BOIN12 Design
#'
#' @param x An object of class \code{comb_boin12_oc} from
#'   \code{\link{sim_comb_boin12}}.
#'
#' @param digits Integer scalar. Decimal places. Defaults to 1.
#'
#' @param percent Logical scalar. Show the patients, DLTs and responses at
#'   each combination as percentages of their totals. Defaults to
#'   \code{FALSE}.
#'
#' @param kable Logical scalar. Print the table with \code{knitr::kable()}.
#'   Defaults to \code{FALSE}.
#'
#' @param kable_format Character scalar. Format passed to
#'   \code{knitr::kable()}. Defaults to \code{"pipe"}.
#'
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
#' p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
#' oc <- sim_comb_boin12(target = 0.33, target_resp = 0.2, p_true = p_true,
#'                       p_resp = p_resp, n_cohort = 14, cohort_size = 3,
#'                       n_trials = 100, est_method = "observed", seed = 1)
#' print(oc)
#' print(oc, percent = TRUE)
#'
#' @seealso \code{\link{sim_comb_boin12}}
#'
#' @export
print.comb_boin12_oc <- function(x, digits = 1, percent = FALSE,
                                 kable = FALSE, kable_format = "pipe", ...) {

  check_flag(percent, "percent")
  tab <- comb_boin12_oc_table(x, digits, percent = percent)

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    print(knitr::kable(tab, format = kable_format, row.names = FALSE))
    return(invisible(x))
  }

  set <- x$settings
  design <- if (is.null(x$design)) "Comb-BOIN12" else x$design
  cat(design, " operating characteristics\n", sep = "")
  cat("  DLT limit            : ", format(x$target), "\n", sep = "")
  cat("  response limit       : ", format(x$target_resp), "\n", sep = "")
  cat("  dose levels          : ", x$n_rows, " of drug A, ", x$n_cols,
      " of drug B\n", sep = "")
  cat("  trials               : ", x$n_trials, "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")
  cat("  lambda_e / lambda_d  : ", format(round(x$lambda_e, 4)), " / ",
      format(round(x$lambda_d, 4)), "\n", sep = "")
  cat("  utility scores       : 100, ", format(round(set$psi00, 4)), ", ",
      format(round(set$psi11, 4)), ", 0; benchmark ",
      format(round(x$u_benchmark, 2)), "\n", sep = "")
  cat("  ties / no candidate  : ", set$tie_break, " / ", set$no_candidate,
      "\n", sep = "")
  cat("  OBD selection        : ", x$est_method, " response, ", x$obd_bound,
      " bound\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  cat("\nNo OBD selected                      : ",
      format(round(x$percent_no_obd, digits)), "%\n", sep = "")
  cat("Trials stopped early                 : ",
      format(round(x$percent_early_stop, digits)), "%\n", sep = "")
  cat("Average number of patients           : ",
      format(round(x$total_n_pts, digits)), "\n", sep = "")

  if (nrow(x$true_obd) == 0L) {
    cat("\nTrue OBD                             : none\n")
  } else {
    cat("\nTrue OBD (highest true utility)      : ",
        paste0("A", x$true_obd[, "dose_a"], "B", x$true_obd[, "dose_b"],
               collapse = ", "), "\n", sep = "")
    cat("  Selected as the OBD                : ",
        format(round(x$pct_true_obd_selected, digits)), "%\n", sep = "")
    cat("  Patients treated there             : ",
        format(round(x$pct_patients_true_obd, digits)), "%\n", sep = "")
  }
  overdose <- x$overdose
  cat("Combinations above the DLT limit     : ",
      nrow(overdose$combinations), " of ", x$n_rows * x$n_cols, "\n",
      sep = "")
  cat("  Selected as the OBD                : ",
      format(round(overdose$pct_selected, digits)), "%\n", sep = "")
  cat("  Patients treated there             : ",
      format(round(overdose$pct_patients, digits)), "%\n", sep = "")
  if (!is.null(x$duration)) {
    cat("Average trial duration               : ",
        format(round(x$duration, digits)), "\n", sep = "")
  }

  invisible(x)
}
