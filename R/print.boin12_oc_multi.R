#' Print the Operating Characteristics of the BOIN12 Design Across Scenarios
#'
#' @description
#'   Display the summary produced by \code{\link{sim_boin12_multi}} as one
#'   table with a block of rows per scenario.
#'
#' @inheritParams print.boin12_oc
#'
#' @param x
#'   An object of class \code{boin12_oc_multi}.
#'
#' @return
#'   \code{x}, invisibly.
#'
#' @examples
#' scenarios <- list(
#'   "Plateau" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                    p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
#'   "Increasing" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                       p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
#' )
#' oc <- sim_boin12_multi(
#'   target = 0.35, target_resp = 0.25, scenarios = scenarios,
#'   n_cohort = 12, cohort_size = 3, n_trials = 200, seed = 1
#' )
#' print(oc, percent = TRUE)
#'
#' @export
print.boin12_oc_multi <- function(x, digits = 1, percent = FALSE,
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
  design <- if (is.null(first$design)) "BOIN12" else first$design
  cat(design, " operating characteristics across ", length(x$scenario_names),
      " scenarios\n", sep = "")
  cat("  DLT limit            : ", format(first$target), "\n", sep = "")
  cat("  response limit       : ", format(first$target_resp), "\n", sep = "")
  cat("  trials each          : ", first$n_trials, "\n", sep = "")
  cat("  max sample size      : ", first$settings$max_total_pts, "\n", sep = "")
  cat("  ties / no candidate  : ", first$settings$tie_break, " / ",
      first$settings$no_candidate, "\n", sep = "")
  cat("  OBD selection        : utility by ", first$est_method, "\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  invisible(x)
}
