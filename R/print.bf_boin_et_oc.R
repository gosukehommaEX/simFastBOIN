#' Print the Operating Characteristics of the BF-BOIN-ET Design
#'
#' @description
#'   Display the summary produced by \code{\link{sim_bf_boin_et}} as a table
#'   with one column per dose.
#'
#' @inheritParams print.boin_et_oc
#'
#' @param x
#'   An object of class \code{bf_boin_et_oc}.
#'
#' @param percent
#'   Logical scalar. Express the average numbers of patients, backfilled
#'   patients, DLTs and responses at each dose as percentages of the trial
#'   totals. Defaults to \code{FALSE}.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' oc <- sim_bf_boin_et(
#'   target = 0.30, target_resp = 0.25,
#'   p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
#'   p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70),
#'   n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
#'   n_trials = 200, seed = 1
#' )
#' print(oc, digits = 2)
#' print(oc, percent = TRUE)
#'
#' @export
print.bf_boin_et_oc <- function(x, digits = 1, percent = FALSE, kable = FALSE,
                                kable_format = "pipe", ...) {

  if (!is.logical(percent) || length(percent) != 1L || is.na(percent)) {
    stop("'percent' must be TRUE or FALSE", call. = FALSE)
  }

  tab <- bf_boin_et_table(x, percent)
  rownames(tab) <- tab$Item
  tab <- as.matrix(tab[, -1L, drop = FALSE])
  tab <- round(tab, digits)

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    print(knitr::kable(tab, format = kable_format, row.names = TRUE))
    return(invisible(x))
  }

  set <- x$settings
  cat("BF-BOIN-ET operating characteristics\n")
  cat("  target DLT rate      : ", format(x$target), "\n", sep = "")
  cat("  target response rate : ", format(x$target_resp), "\n", sep = "")
  cat("  trials               : ", x$n_trials, "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  patients per dose    : ", set$n_earlystop, " at most with backfill\n",
      sep = "")
  cat("  lambda1 / lambda2    : ", format(round(x$lambda1, 4)), " / ",
      format(round(x$lambda2, 4)), "\n", sep = "")
  cat("  eta1                 : ", format(round(x$eta1, 4)), "\n", sep = "")
  cat("\n")

  print(tab, na.print = "")

  cat("\nTrials stopped early: ", format(round(x$percent_early_stop, digits)),
      "%\n", sep = "")
  cat("Average trial duration: ", format(round(x$duration, digits)), "\n",
      sep = "")

  invisible(x)
}
