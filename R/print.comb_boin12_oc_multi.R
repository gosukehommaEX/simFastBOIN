#' Print the Operating Characteristics of Comb-BOIN12 Across Scenarios
#'
#' @param x An object of class \code{comb_boin12_oc_multi} from
#'   \code{\link{sim_comb_boin12_multi}}.
#'
#' @param digits Integer scalar. Decimal places. Defaults to 1.
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
#' @details
#'   One row per scenario: the true OBD combinations, the percentage of
#'   trials selecting one of them and of the patients treated there, the
#'   percentage of trials selecting no OBD combination, the percentages
#'   selecting and treating combinations above the DLT limit, and the
#'   average numbers of patients, DLTs and responses.
#'
#' @examples
#' scenarios <- list(
#'   "Scenario 16" = list(
#'     p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2),
#'     p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)),
#'   "Scenario 17" = list(
#'     p_true = matrix(c(0.05, 0.10, 0.15, 0.20, 0.40, 0.50), nrow = 2),
#'     p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.55), nrow = 2))
#' )
#' oc <- sim_comb_boin12_multi(target = 0.33, target_resp = 0.2,
#'                             scenarios = scenarios, n_cohort = 14,
#'                             cohort_size = 3, n_trials = 100,
#'                             est_method = "observed", seed = 1)
#' print(oc)
#'
#' @seealso \code{\link{sim_comb_boin12_multi}}
#'
#' @export
print.comb_boin12_oc_multi <- function(x, digits = 1, kable = FALSE,
                                       kable_format = "pipe", ...) {

  tab <- comb_boin12_multi_table(x$results, x$scenario_names, digits = digits)

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    print(knitr::kable(tab, format = kable_format, row.names = FALSE))
    return(invisible(x))
  }

  first <- x$results[[1L]]
  design <- if (is.null(first$design)) "Comb-BOIN12" else first$design
  cat(design, " operating characteristics across ", length(x$scenario_names),
      " scenarios\n", sep = "")
  cat("  DLT limit            : ", format(first$target), "\n", sep = "")
  cat("  response limit       : ", format(first$target_resp), "\n", sep = "")
  cat("  dose levels          : ", x$n_rows, " of drug A, ", x$n_cols,
      " of drug B\n", sep = "")
  cat("  trials each          : ", first$n_trials, "\n", sep = "")
  cat("  max sample size      : ", first$settings$max_total_pts, "\n",
      sep = "")
  cat("  ties / no candidate  : ", first$settings$tie_break, " / ",
      first$settings$no_candidate, "\n", sep = "")
  cat("  OBD selection        : ", first$est_method, " response, ",
      first$obd_bound, " bound\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  invisible(x)
}
