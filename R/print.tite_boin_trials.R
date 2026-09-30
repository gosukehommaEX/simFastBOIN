#' Print Simulated TITE-BOIN Trials
#'
#' @description
#'   Show a compact summary of the object returned by
#'   \code{\link{tite_boin_simulate}} instead of printing the full trial
#'   matrices.
#'
#' @param x
#'   An object of class \code{tite_boin_trials}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' trials <- tite_boin_simulate(
#'   target = 0.30,
#'   p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   window = 3,
#'   accrual_rate = 2,
#'   n_trials = 100,
#'   seed = 1
#' )
#' print(trials)
#'
#' @export
print.tite_boin_trials <- function(x, ...) {

  set <- x$settings
  n_doses <- length(set$p_true)
  method_label <- if (identical(set$method, "imputation")) {
    "single mean imputation (Yuan et al., 2018)"
  } else {
    "effective sample size (Lin and Yuan, 2020)"
  }

  cat("Simulated TITE-BOIN trials\n")
  cat("  method          : ", method_label, "\n", sep = "")
  cat("  trials          : ", set$n_trials, "\n", sep = "")
  cat("  doses           : ", n_doses, "\n", sep = "")
  cat("  target DLT rate : ", format(set$target), "\n", sep = "")
  cat("  cohorts         : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size : ", set$max_total_pts, "\n", sep = "")
  cat("  window          : ", format(set$window), "\n", sep = "")
  cat("  accrual         : ", format(set$accrual_rate),
      " patients per unit of time (", set$accrual, ")\n", sep = "")
  cat("\n")

  avg <- rbind(
    "Patients treated" = colMeans(x$n_pts),
    "DLTs observed" = colMeans(x$n_tox)
  )
  colnames(avg) <- paste0("DL", seq_len(n_doses))
  cat("Average per trial:\n")
  print(round(avg, 2))

  cat("\nTrial duration        : mean ", format(round(mean(x$duration), 2)),
      ", median ", format(round(stats::median(x$duration), 2)), "\n", sep = "")
  cat("Trials with suspension: ",
      format(round(100 * mean(x$n_suspensions > 0), 1)), "%\n", sep = "")
  cat("Time suspended        : mean ",
      format(round(mean(x$time_suspended), 2)), " per trial\n", sep = "")

  cat("\nStopping reason (%):\n")
  print(round(100 * table(x$stop_reason) / set$n_trials, 1))

  invisible(x)
}
