#' Print Simulated BOIN Drug-Combination Trials
#'
#' @description
#'   Show a compact summary of the object returned by
#'   \code{\link{comb_boin_simulate}} instead of printing the full trial arrays.
#'
#' @param x
#'   An object of class \code{comb_boin_trials}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' p_true <- matrix(c(0.05, 0.15, 0.30,
#'                    0.15, 0.30, 0.45), nrow = 2, byrow = TRUE)
#'
#' trials <- comb_boin_simulate(
#'   target = 0.30,
#'   p_true = p_true,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_trials = 100,
#'   seed = 1
#' )
#' print(trials)
#'
#' @export
print.comb_boin_trials <- function(x, ...) {

  set <- x$settings
  n_rows <- nrow(set$p_true)
  n_cols <- ncol(set$p_true)
  labels <- list(paste0("A", seq_len(n_rows)), paste0("B", seq_len(n_cols)))

  cat("Simulated BOIN drug-combination trials\n")
  cat("  trials          : ", set$n_trials, "\n", sep = "")
  cat("  dose levels     : ", n_rows, " of drug A, ", n_cols, " of drug B\n",
      sep = "")
  cat("  target DLT rate : ", format(set$target), "\n", sep = "")
  cat("  cohorts         : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size : ", set$max_total_pts, "\n", sep = "")
  cat("\n")

  pts <- matrix(apply(x$n_pts, c(2L, 3L), mean), n_rows, n_cols,
                dimnames = labels)
  tox <- matrix(apply(x$n_tox, c(2L, 3L), mean), n_rows, n_cols,
                dimnames = labels)

  cat("Average patients treated per trial:\n")
  print(round(pts, 2))
  cat("\nAverage DLTs observed per trial:\n")
  print(round(tox, 2))

  cat("\nStopping reason (%):\n")
  print(round(100 * table(x$stop_reason) / set$n_trials, 1))

  invisible(x)
}
