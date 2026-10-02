#' Print Operating Characteristics of a Backfill Design Across Scenarios
#'
#' @description
#'   Display the result of \code{\link{sim_bf_boin_multi}} or
#'   \code{\link{sim_be_boin_multi}} as one table with seven rows per scenario.
#'
#' @param x
#'   An object of class \code{backfill_oc_multi}.
#'
#' @param digits
#'   Integer scalar. Number of decimal places. Defaults to 1.
#'
#' @param percent
#'   Logical scalar. Show patients, backfilled patients and DLTs as percentages
#'   of their totals. Defaults to \code{FALSE}.
#'
#' @param kable
#'   Logical scalar. Print only the table, formatted by \code{knitr::kable()}.
#'   Defaults to \code{FALSE}.
#'
#' @param kable_format
#'   Character scalar passed to \code{knitr::kable()} as \code{format}.
#'   Defaults to \code{"pipe"}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' scenarios <- list(
#'   list(name = "Rising", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
#'        p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58)),
#'   list(name = "Plateau", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
#'        p_resp = c(0.30, 0.32, 0.35, 0.36, 0.36))
#' )
#' oc <- sim_bf_boin_multi(
#'   target = 0.25, scenarios = scenarios, n_cohort = 10, cohort_size = 3,
#'   window = 1, accrual_rate = 3, n_trials = 100, seed = 1
#' )
#' print(oc, percent = TRUE)
#'
#' @seealso \code{\link{sim_bf_boin_multi}}, \code{\link{sim_be_boin_multi}}
#'
#' @export
print.backfill_oc_multi <- function(x, digits = 1, percent = FALSE,
                                    kable = FALSE, kable_format = "pipe", ...) {

  if (!is.logical(percent) || length(percent) != 1L || is.na(percent)) {
    stop("'percent' must be TRUE or FALSE", call. = FALSE)
  }

  tab <- backfill_oc_multi_table(x$results, x$scenario_names, x$n_doses,
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
  set <- first$settings
  design_label <- if (identical(set$design, "be_boin")) "BE-BOIN" else "BF-BOIN"
  cat(design_label, " operating characteristics across ",
      length(x$scenario_names), " scenarios\n", sep = "")
  cat("  target DLT rate : ", format(first$target), "\n", sep = "")
  cat("  trials each     : ", first$n_trials, "\n", sep = "")
  cat("  escalation size : ", set$max_total_pts, "\n", sep = "")
  cat("  backfill cap    : ", set$n_cap, " patients per dose\n", sep = "")
  cat("  window          : ", format(set$window), " (DLT), ",
      format(set$resp_window), " (response)\n", sep = "")
  cat("  accrual         : ", format(set$accrual_rate),
      " patients per unit of time (", set$accrual, ")\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  invisible(x)
}
