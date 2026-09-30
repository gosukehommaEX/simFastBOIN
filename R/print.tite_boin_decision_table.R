#' Print a TITE-BOIN Decision Table
#'
#' @description
#'   Display the decision table produced by \code{\link{tite_boin_decision_table}}
#'   in the layout of the published tables: one row per number of patients
#'   treated, number of DLTs and range of pending patients sharing a decision,
#'   with the boundaries on the follow-up statistic in the Escalate, Stay and
#'   De-escalate columns.
#'
#' @param x
#'   An object of class \code{tite_boin_decision_table}.
#'
#' @param cohort_size
#'   Integer scalar or \code{NULL}. When supplied, only the sample sizes that are
#'   multiples of \code{cohort_size} are shown, which is how the published tables
#'   are laid out.
#'
#' @param digits
#'   Integer scalar. Number of decimal places of the boundaries. Defaults to 2,
#'   as in the published tables.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @details
#'   A boundary shown as \code{">= 2.15"} under Escalate means that the dose is
#'   escalated when the follow-up statistic is at least 2.15; the statistic is
#'   STFT for \code{method = "imputation"} and ESS for \code{method = "ess"}.
#'   "Suspend if >= 4.23" under Escalate marks an escalation blocked by
#'   \code{min_completed}.
#'
#'   A data frame that no longer has the columns or attributes of a decision
#'   table, for example after selecting columns, is printed as a plain data
#'   frame.
#'
#' @examples
#' decisions <- tite_boin_decision_table(target = 0.3, max_n = 15)
#'
#' # Table S1 of Yuan et al. (2018)
#' print(decisions, cohort_size = 3)
#'
#' # Boundaries with three decimal places
#' print(decisions[decisions$n == 9, ], digits = 3)
#'
#' @seealso \code{\link{tite_boin_decision_table}},
#'   \code{\link{plot.tite_boin_decision_table}}
#'
#' @export
print.tite_boin_decision_table <- function(x, cohort_size = NULL, digits = 2, ...) {

  required <- c("n", "n_tox", "n_pending", "decision", "esc_bound", "deesc_bound")
  if (is.null(attr(x, "statistic")) || !all(required %in% names(x))) {
    return(NextMethod())
  }

  # Validate before printing anything, so that a rejected argument does not
  # leave a half-written table behind.
  if (!is.null(cohort_size)) check_count(cohort_size, "cohort_size", 1L)
  check_count(digits, "digits", 0L)

  shown <- tite_boin_display_rows(x, cohort_size = cohort_size,
                                  digits = as.integer(digits))
  if (nrow(shown) == 0L) {
    cat("No sample size in the table is a multiple of the cohort size.\n")
    return(invisible(x))
  }

  statistic <- attr(x, "statistic")
  method_label <- if (identical(attr(x, "method"), "imputation")) {
    "single mean imputation (Yuan et al., 2018)"
  } else {
    "effective sample size (Lin and Yuan, 2020)"
  }
  cat("TITE-BOIN decision table with target ", format(attr(x, "target")), "\n",
      sep = "")
  cat("Method: ", method_label, "; boundaries refer to ", statistic, "\n\n",
      sep = "")

  print(shown, row.names = FALSE, right = FALSE)

  cat("\nPatients = treated at the current dose, DLTs = DLTs observed,",
      "Pending = assessment not yet completed\n")
  cat("STFT = total follow-up time of the pending patients divided by the",
      "length of the assessment window\n")
  if (statistic == "ESS") {
    cat("ESS = number of patients with a completed assessment plus STFT\n")
  }
  ratio <- attr(x, "max_pending_ratio")
  if (!is.null(ratio) && ratio < 1) {
    cat("Accrual is suspended when more than ", format(100 * ratio),
        "% of the patients are pending\n", sep = "")
  }
  min_completed <- attr(x, "min_completed")
  if (!is.null(min_completed) && min_completed > 0L) {
    cat("While patients are pending, escalation requires at least ",
        min_completed, " patients with a completed assessment\n", sep = "")
  }
  cat("Y & Elim = de-escalate and eliminate this dose and all higher doses\n")

  invisible(x)
}
