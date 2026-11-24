#' Print Simulated Waterfall Trials
#'
#' @description
#'   Print the design and the averages over the trials of a
#'   \code{waterfall_trials} object.
#'
#' @param x
#'   An object of class \code{waterfall_trials}, as returned by
#'   \code{\link{waterfall_simulate}}.
#'
#' @param ...
#'   Not used.
#'
#' @return
#'   \code{x}, invisibly.
#'
#' @examples
#' p_true <- matrix(c(0.03, 0.10, 0.28,
#'                    0.10, 0.30, 0.50), nrow = 2, byrow = TRUE)
#'
#' trials <- waterfall_simulate(
#'   target = 0.30,
#'   p_true = p_true,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#' print(trials)
#'
#' @export
print.waterfall_trials <- function(x, ...) {

  set <- x$settings
  n_rows <- nrow(set$p_true)
  n_cols <- ncol(set$p_true)
  labels <- list(paste0("A", seq_len(n_rows)), paste0("B", seq_len(n_cols)))

  cat("Simulated waterfall drug-combination trials\n")
  cat("  trials          : ", set$n_trials, "\n", sep = "")
  cat("  dose levels     : ", n_rows, " of drug A, ", n_cols, " of drug B\n",
      sep = "")
  cat("  target DLT rate : ", format(set$target), "\n", sep = "")
  cat("  cohorts         : at most ", paste(set$n_cohort, collapse = "/"),
      " per subtrial, of size ", set$cohort_size, "\n", sep = "")
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

  cat("\nReason for ending (%):\n")
  print(round(100 * table(x$stop_reason) / set$n_trials, 1))

  invisible(x)
}
