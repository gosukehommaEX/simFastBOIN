#' Print the Operating Characteristics of the BOIN-ETC Design Across
#' Scenarios
#'
#' @param x An object of class \code{boin_etc_oc_multi} from
#'   \code{\link{sim_boin_etc_multi}}.
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
#' @examples
#' scenarios <- list(
#'   "Scenario 1" = list(
#'     p_true = matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2),
#'     p_resp = matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)),
#'   "Scenario 2" = list(
#'     p_true = matrix(c(0.05, 0.15, 0.15, 0.35, 0.30, 0.55), nrow = 2),
#'     p_resp = matrix(c(0.25, 0.60, 0.60, 0.65, 0.65, 0.70), nrow = 2))
#' )
#' res <- sim_boin_etc_multi(target = 0.30, target_resp = 0.60,
#'                           scenarios = scenarios, n_cohort = 8,
#'                           cohort_size = 3, n_trials = 100)
#' print(res)
#'
#' @seealso \code{\link{sim_boin_etc_multi}}
#'
#' @export
print.boin_etc_oc_multi <- function(x, digits = 1, kable = FALSE,
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
  set <- first$settings
  cat(first$design, " operating characteristics across ",
      length(x$scenario_names), " scenarios\n", sep = "")
  cat("  target DLT rate      : ", format(first$target), "\n", sep = "")
  cat("  target response rate : ", format(first$target_resp), "\n", sep = "")
  cat("  dose levels          : ", x$n_rows, " of drug A, ", x$n_cols,
      " of drug B\n", sep = "")
  cat("  trials each          : ", first$n_trials, "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")
  cat("  subtrial end (n_stop): ", set$n_earlystop, "\n", sep = "")
  if (set$subtrials == "etc1") {
    cat("  start of a row       : ", set$etc1_start, "\n", sep = "")
  }
  cat("  rule 4               : ", set$rule4, "\n", sep = "")
  cat("  utility scores       : 100, ", format(round(set$psi00, 4)), ", ",
      format(round(set$psi11, 4)), ", 0; ties ", first$odc_tie, "\n",
      sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  invisible(x)
}
