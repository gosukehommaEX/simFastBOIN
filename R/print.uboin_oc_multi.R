#' Print Operating Characteristics of the U-BOIN Design Across Scenarios
#'
#' @description
#'   Print the design and the cross-scenario summary table of a
#'   \code{uboin_oc_multi} object.
#'
#' @param x
#'   An object of class \code{uboin_oc_multi}, as returned by
#'   \code{\link{sim_uboin_multi}}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @param percent
#'   Logical scalar. Show the numbers of patients, DLTs and responses as
#'   percentages of their totals. Defaults to \code{FALSE}.
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
#'   "Scenario 1" = list(p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
#'                       p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65)),
#'   "Scenario 3" = list(p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
#'                       p_resp = c(0.08, 0.46, 0.25, 0.20, 0.10))
#' )
#' oc <- sim_uboin_multi(
#'   tox_upper = 0.30,
#'   resp_lower = 0.20,
#'   scenarios = scenarios,
#'   n_cohort = 18,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#'
#' print(oc)
#'
#' @export
print.uboin_oc_multi <- function(x, digits = 1, percent = FALSE,
                                 kable = FALSE, kable_format = "pipe", ...) {

  check_flag(percent, "percent")

  tab <- boin_et_multi_table(x$results, x$scenario_names, x$n_doses,
                             percent = percent, digits = digits)
  tab[] <- lapply(tab, function(column) {
    if (is.numeric(column)) {
      out <- format(column, nsmall = digits, trim = TRUE)
      out[is.na(column)] <- ""
      out
    } else {
      column
    }
  })

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    print(knitr::kable(tab, format = kable_format, row.names = FALSE))
    return(invisible(x))
  }

  first <- x$results[[1L]]
  cat("U-BOIN operating characteristics across ", length(x$scenario_names),
      " scenarios\n", sep = "")
  cat("  DLT upper limit      : ", format(first$tox_upper), "\n", sep = "")
  cat("  response lower limit : ", format(first$resp_lower), "\n", sep = "")
  cat("  trials each          : ", first$n_trials, "\n", sep = "")
  cat("  max sample size      : ", first$settings$max_total_pts, "\n", sep = "")
  utilities <- paste(format(first$utility, trim = TRUE), collapse = ", ")
  cat("  utilities            : ", utilities, "\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  invisible(x)
}
