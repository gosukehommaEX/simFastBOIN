#' Print the Operating Characteristics of BF-BOIN-ET Across Scenarios
#'
#' @description
#'   Display the table produced by \code{\link{sim_bf_boin_et_multi}}, with
#'   seven rows per scenario.
#'
#' @inheritParams print.bf_boin_et_oc
#'
#' @param x
#'   An object of class \code{bf_boin_et_oc_multi}.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' scenarios <- list(
#'   "Example" = list(p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
#'                    p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70)),
#'   "Plateau" = list(p_true = c(0.05, 0.10, 0.15, 0.30, 0.45),
#'                    p_resp = c(0.10, 0.30, 0.45, 0.50, 0.52))
#' )
#' res <- sim_bf_boin_et_multi(
#'   target = 0.30, target_resp = 0.25, scenarios = scenarios,
#'   n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
#'   n_trials = 100, seed = 1
#' )
#' print(res, digits = 2)
#'
#' @export
print.bf_boin_et_oc_multi <- function(x, digits = 1, percent = FALSE,
                                      kable = FALSE, kable_format = "pipe",
                                      ...) {

  if (!is.logical(percent) || length(percent) != 1L || is.na(percent)) {
    stop("'percent' must be TRUE or FALSE", call. = FALSE)
  }

  tab <- bf_boin_et_multi_table(x$results, x$scenario_names,
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
  cat("BF-BOIN-ET operating characteristics across ",
      length(x$scenario_names), " scenarios\n", sep = "")
  cat("  target DLT rate      : ", format(first$target), "\n", sep = "")
  cat("  target response rate : ", format(first$target_resp), "\n", sep = "")
  cat("  trials each          : ", first$n_trials, "\n", sep = "")
  cat("  cohorts              : ", first$settings$n_cohort, " of size ",
      paste(unique(first$settings$cohort_size), collapse = "/"), "\n",
      sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  invisible(x)
}
