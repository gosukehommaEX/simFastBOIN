#' Print Simulated U-BOIN Trials
#'
#' @description
#'   Print the design and the averages over the trials of a
#'   \code{uboin_trials} object.
#'
#' @param x
#'   An object of class \code{uboin_trials}, as returned by
#'   \code{\link{uboin_simulate}}.
#'
#' @param ...
#'   Not used.
#'
#' @return
#'   \code{x}, invisibly.
#'
#' @examples
#' trials <- uboin_simulate(
#'   tox_upper = 0.30,
#'   resp_lower = 0.20,
#'   p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
#'   p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65),
#'   n_cohort = 18,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#' print(trials)
#'
#' @export
print.uboin_trials <- function(x, ...) {

  set <- x$settings
  cat("Simulated U-BOIN trials\n")
  cat("  trials               : ", set$n_trials, "\n", sep = "")
  cat("  doses                : ", length(set$p_true), "\n", sep = "")
  cat("  DLT upper limit      : ", format(set$tox_upper), "\n", sep = "")
  cat("  response lower limit : ", format(set$resp_lower), "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ", set$cohort_size,
      "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")
  cat("\n")

  tab <- rbind(colMeans(x$n_pts), colMeans(x$n_tox), colMeans(x$n_resp))
  rownames(tab) <- c("Patients treated", "DLTs observed", "Responses observed")
  cat("Averages per trial:\n")
  print(round(tab, 2))

  cat("\nTrials reaching stage II: ", format(round(100 * mean(x$stage2), 1)),
      "%\n", sep = "")
  cat("Reason for ending (%):\n")
  reasons <- uboin_stop_reasons()
  print(round(100 * table(factor(x$stop_reason, levels = reasons)) /
                set$n_trials, 1))

  invisible(x)
}
