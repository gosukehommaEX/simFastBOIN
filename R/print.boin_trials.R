#' Print Simulated BOIN Trials
#'
#' @description
#'   Show a compact summary of the object returned by \code{\link{boin_simulate}}
#'   instead of printing the full trial matrices.
#'
#' @param x
#'   An object of class \code{boin_trials}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' trials <- boin_simulate(
#'   target = 0.30,
#'   p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_trials = 100,
#'   seed = 1
#' )
#' print(trials)
#'
#' @export
print.boin_trials <- function(x, ...) {

  set <- x$settings
  n_doses <- length(set$p_true)

  cat("Simulated BOIN trials\n")
  cat("  trials          : ", set$n_trials, "\n", sep = "")
  cat("  doses           : ", n_doses, "\n", sep = "")
  cat("  target DLT rate : ", format(set$target), "\n", sep = "")
  cat("  cohorts         : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size : ", set$max_total_pts, "\n", sep = "")
  cat("\n")

  avg <- rbind(
    "Patients treated" = colMeans(x$n_pts),
    "DLTs observed" = colMeans(x$n_tox)
  )
  colnames(avg) <- paste0("DL", seq_len(n_doses))
  cat("Average per trial:\n")
  print(round(avg, 2))

  cat("\nStopping reason (%):\n")
  print(round(100 * table(x$stop_reason) / set$n_trials, 1))

  invisible(x)
}
