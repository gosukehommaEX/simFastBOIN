#' Print Operating Characteristics of the BOIN Drug-Combination Design Across Scenarios
#'
#' @description
#'   Display the table produced by \code{\link{sim_comb_boin_multi}}, one row per
#'   scenario. The full results of a scenario are printed by printing its
#'   element of \code{x$results}.
#'
#' @param x
#'   An object of class \code{comb_boin_oc_multi}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @param kable
#'   Logical scalar. Return the table through \code{knitr::kable()} rather than
#'   printing it directly. Requires the \pkg{knitr} package. Defaults to
#'   \code{FALSE}.
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
#' oc <- sim_comb_boin_multi(
#'   target = 0.30,
#'   scenarios = list(
#'     "MTD in the middle" = matrix(c(0.10, 0.20, 0.30,
#'                                    0.20, 0.30, 0.45), 2, byrow = TRUE),
#'     "All toxic"         = matrix(c(0.40, 0.50, 0.60,
#'                                    0.50, 0.60, 0.70), 2, byrow = TRUE)
#'   ),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#'
#' print(oc)
#'
#' @export
print.comb_boin_oc_multi <- function(x, digits = 1, kable = FALSE,
                                     kable_format = "pipe", ...) {

  tab <- comb_oc_multi_table(x$results, x$scenario_names, digits = digits)

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    print(knitr::kable(tab, format = kable_format, row.names = FALSE))
    return(invisible(x))
  }

  first <- x$results[[1L]]
  cat("BOIN drug-combination operating characteristics across ",
      length(x$scenario_names), " scenarios\n", sep = "")
  cat("  target DLT rate : ", format(first$target), "\n", sep = "")
  cat("  dose levels     : ", x$n_rows, " of drug A, ", x$n_cols,
      " of drug B\n", sep = "")
  cat("  trials each     : ", first$n_trials, "\n", sep = "")
  cat("  max sample size : ", first$settings$max_total_pts, "\n", sep = "")
  cat("  overdose cutoff : ", format(first$overdose$cutoff), "\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  invisible(x)
}
