#' Print the Operating Characteristics of BOIN-ET Across Scenarios
#'
#' @description
#'   Display the table produced by \code{\link{sim_boin_et_multi}}, with six
#'   rows per scenario.
#'
#' @inheritParams print.boin_et_oc
#'
#' @param x
#'   An object of class \code{boin_et_oc_multi}.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' scenarios <- list(
#'   list(name = "Monotone",
#'        p_true = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.40),
#'        p_resp = c(0.05, 0.10, 0.25, 0.40, 0.65, 0.70)),
#'   list(name = "Peaked",
#'        p_true = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
#'        p_resp = c(0.30, 0.70, 0.60, 0.50, 0.40, 0.30))
#' )
#' res <- sim_boin_et_multi(
#'   target = 0.30, target_resp = 0.60, scenarios = scenarios,
#'   n_cohort = 12, cohort_size = 3, n_trials = 300, seed = 1
#' )
#' print(res)
#' print(res, percent = TRUE)
#'
#' @seealso \code{\link{sim_boin_et_multi}}
#'
#' @export
print.boin_et_oc_multi <- function(x, digits = 1, percent = FALSE,
                                   kable = FALSE, kable_format = "pipe", ...) {

  if (!is.logical(percent) || length(percent) != 1L || is.na(percent)) {
    stop("'percent' must be TRUE or FALSE", call. = FALSE)
  }

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
  cat("BOIN-ET operating characteristics across ", length(x$scenario_names),
      " scenarios\n", sep = "")
  cat("  target DLT rate      : ", format(first$target), "\n", sep = "")
  cat("  target response rate : ", format(first$target_resp), "\n", sep = "")
  cat("  trials each          : ", first$n_trials, "\n", sep = "")
  cat("  max sample size      : ", first$settings$max_total_pts, "\n", sep = "")
  cat("  OBD selection        : ", first$obd_method, " with ",
      first$est_method, " response rates\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  invisible(x)
}
