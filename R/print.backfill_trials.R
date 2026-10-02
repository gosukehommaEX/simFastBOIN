#' Print Simulated Trials of a BOIN Design with Backfilling
#'
#' @description
#'   Summarize the trials returned by \code{\link{bf_boin_simulate}} or
#'   \code{\link{be_boin_simulate}}: the design, the average numbers of
#'   patients, backfilled patients, DLTs and responses at each dose, the trial
#'   duration, the waiting of the dose escalation and the reasons for stopping.
#'
#' @param x
#'   An object of class \code{backfill_trials}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' trials <- bf_boin_simulate(
#'   target = 0.25,
#'   p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
#'   p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58),
#'   n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
#'   n_trials = 100, seed = 1
#' )
#' print(trials)
#'
#' @seealso \code{\link{bf_boin_simulate}}, \code{\link{be_boin_simulate}}
#'
#' @export
print.backfill_trials <- function(x, ...) {

  set <- x$settings
  n_doses <- length(set$p_true)
  design_label <- if (identical(set$design, "be_boin")) {
    "BE-BOIN (Chen et al., 2026)"
  } else {
    "BF-BOIN (Zhao et al., 2024)"
  }

  cat("Simulated ", design_label, " trials\n", sep = "")
  cat("  trials          : ", set$n_trials, "\n", sep = "")
  cat("  doses           : ", n_doses, "\n", sep = "")
  cat("  target DLT rate : ", format(set$target), "\n", sep = "")
  cat("  cohorts         : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  escalation size : ", set$max_total_pts, "\n", sep = "")
  cat("  backfill cap    : ", set$n_cap, " patients per dose\n", sep = "")
  cat("  window          : ", format(set$window), " (DLT), ",
      format(set$resp_window), " (response)\n", sep = "")
  cat("  accrual         : ", format(set$accrual_rate),
      " patients per unit of time (", set$accrual, ")\n", sep = "")
  cat("\n")

  avg <- rbind(
    "Patients treated" = colMeans(x$n_pts),
    "  of which backfilled" = colMeans(x$n_bf),
    "DLTs observed" = colMeans(x$n_tox),
    "Responses" = colMeans(x$n_resp)
  )
  colnames(avg) <- paste0("DL", seq_len(n_doses))
  cat("Average per trial:\n")
  print(round(avg, 2))

  cat("\nPatients per trial    : mean ", format(round(mean(rowSums(x$n_pts)), 2)),
      ", of which backfilled ", format(round(mean(rowSums(x$n_bf)), 2)), "\n",
      sep = "")
  cat("Trials with backfill  : ",
      format(round(100 * mean(rowSums(x$n_bf) > 0), 1)), "%\n", sep = "")
  cat("Trial duration        : mean ", format(round(mean(x$duration), 2)),
      ", median ", format(round(stats::median(x$duration), 2)), "\n", sep = "")
  cat("Escalation waiting    : ",
      format(round(100 * mean(x$n_suspensions > 0), 1)),
      "% of trials, mean time ", format(round(mean(x$time_suspended), 2)),
      " per trial\n", sep = "")
  if (identical(set$no_slot, "leave")) {
    cat("Patients turned away  : mean ",
        format(round(mean(x$n_turned_away), 2)), " per trial\n", sep = "")
  }

  cat("\nStopping reason (%):\n")
  print(round(100 * table(x$stop_reason) / set$n_trials, 1))

  invisible(x)
}
