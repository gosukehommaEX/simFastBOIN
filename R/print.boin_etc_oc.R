#' Print the Operating Characteristics of the BOIN-ETC Design
#'
#' @param x An object of class \code{boin_etc_oc} from
#'   \code{\link{sim_boin_etc}}.
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
#' p_true <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
#' p_resp <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
#' oc <- sim_boin_etc(target = 0.30, target_resp = 0.60, p_true = p_true,
#'                    p_resp = p_resp, n_cohort = 8, cohort_size = 3,
#'                    n_trials = 100, seed = 1)
#' print(oc)
#' print(oc, percent = TRUE)
#'
#' @seealso \code{\link{sim_boin_etc}}
#'
#' @export
print.boin_etc_oc <- function(x, digits = 1, percent = FALSE, kable = FALSE,
                              kable_format = "pipe", ...) {

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
  cat(x$design, " operating characteristics\n", sep = "")
  cat("  target DLT rate      : ", format(x$target), "\n", sep = "")
  cat("  target response rate : ", format(x$target_resp), "\n", sep = "")
  cat("  dose levels          : ", x$n_rows, " of drug A, ", x$n_cols,
      " of drug B\n", sep = "")
  cat("  trials               : ", x$n_trials, "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")
  cat("  subtrial end (n_stop): ", set$n_earlystop, "\n", sep = "")
  cat("  lambda1/lambda2/eta1 : ", format(x$lambda1), " / ",
      format(x$lambda2), " / ", format(x$eta1), "\n", sep = "")
  if (set$subtrials == "etc1") {
    cat("  start of a row       : ", set$etc1_start, "\n", sep = "")
  }
  cat("  rule 4               : ", set$rule4, "\n", sep = "")
  cat("  utility scores       : 100, ", format(round(set$psi00, 4)), ", ",
      format(round(set$psi11, 4)), ", 0; ties ", x$odc_tie, "\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  cat("\nNo OBD selected                      : ",
      format(round(x$percent_no_obd, digits)), "%\n", sep = "")
  cat("Trials stopped early                 : ",
      format(round(x$percent_early_stop, digits)), "%\n", sep = "")
  cat("Average number of patients           : ",
      format(round(x$total_n_pts, digits)), "\n", sep = "")
  cat("Average number of subtrials          : ",
      format(round(x$mean_subtrials, 2)), "\n", sep = "")

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
  cat("Combinations above the target DLT    : ",
      nrow(overdose$combinations), " of ", x$n_rows * x$n_cols, "\n",
      sep = "")
  cat("  Selected as the OBD                : ",
      format(round(overdose$pct_selected, digits)), "%\n", sep = "")
  cat("  Patients treated there             : ",
      format(round(overdose$pct_patients, digits)), "%\n", sep = "")

  invisible(x)
}
