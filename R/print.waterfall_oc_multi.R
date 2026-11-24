#' Print Operating Characteristics of the Waterfall Design Across Scenarios
#'
#' @description
#'   Print the design and the cross-scenario summary table of a
#'   \code{waterfall_oc_multi} object.
#'
#' @param x
#'   An object of class \code{waterfall_oc_multi}, as returned by
#'   \code{\link{sim_waterfall_multi}}.
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
#' scenarios <- list(
#'   "Contour A1B3, A2B2" = matrix(c(0.03, 0.10, 0.28,
#'                                   0.10, 0.30, 0.50), nrow = 2, byrow = TRUE),
#'   "Contour A1B1" = matrix(c(0.30, 0.40, 0.50,
#'                             0.42, 0.49, 0.55), nrow = 2, byrow = TRUE)
#' )
#'
#' oc <- sim_waterfall_multi(
#'   target = 0.30,
#'   scenarios = scenarios,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#'
#' print(oc)
#'
#' @export
print.waterfall_oc_multi <- function(x, digits = 1, kable = FALSE,
                                     kable_format = "pipe", ...) {

  tab <- waterfall_oc_multi_table(x$results, x$scenario_names,
                                  digits = digits)

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    print(knitr::kable(tab, format = kable_format, row.names = FALSE))
    return(invisible(x))
  }

  first <- x$results[[1L]]
  cat("Waterfall drug-combination operating characteristics across ",
      length(x$scenario_names), " scenarios\n", sep = "")
  cat("  target DLT rate : ", format(first$target), "\n", sep = "")
  cat("  dose levels     : ", x$n_rows, " of drug A, ", x$n_cols,
      " of drug B\n", sep = "")
  cat("  trials each     : ", first$n_trials, "\n", sep = "")
  cat("  max sample size : ", first$settings$max_total_pts, "\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  invisible(x)
}
