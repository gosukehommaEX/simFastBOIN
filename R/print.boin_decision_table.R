#' Print a BOIN Decision Table
#'
#' @description
#'   Display the decision table produced by \code{\link{boin_decision_table}},
#'   optionally restricted to the sample sizes reached at the end of a cohort.
#'
#' @param x
#'   An object of class \code{boin_decision_table}.
#'
#' @param cohort_size
#'   Integer scalar or \code{NULL}. When supplied, only the sample sizes that are
#'   multiples of \code{cohort_size} are shown, which is how the table is
#'   consulted during a trial with equally sized cohorts.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' decisions <- boin_decision_table(target = 0.30, max_n = 18)
#'
#' print(decisions, cohort_size = 3)
#'
#' @seealso \code{\link{boin_decision_table}}, \code{\link{plot.boin_decision_table}}
#'
#' @export
print.boin_decision_table <- function(x, cohort_size = NULL, ...) {

  # Validate before printing anything, so that a rejected argument does not
  # leave a half-written table behind.
  n_pts <- as.integer(colnames(x))
  keep <- seq_along(n_pts)
  if (!is.null(cohort_size)) {
    check_count(cohort_size, "cohort_size", 1L)
    keep <- keep[n_pts %% as.integer(cohort_size) == 0L]
  }

  if (length(keep) == 0L) {
    cat("No sample size in the table is a multiple of the cohort size.\n")
    return(invisible(x))
  }

  out <- unclass(x)[, keep, drop = FALSE]
  names(dimnames(out)) <- c("DLTs", "Patients")
  print(out, quote = FALSE, na.print = "")

  cat("\nE = escalate, S = stay, D = de-escalate,",
      "DE = de-escalate and eliminate this dose and above\n")

  invisible(x)
}
