#' Print Simulated BOIN-ET Trials
#'
#' @description
#'   Show a compact summary of the object returned by
#'   \code{\link{boin_et_simulate}} instead of printing the full trial
#'   matrices.
#'
#' @param x
#'   An object of class \code{boin_et_trials}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' trials <- boin_et_simulate(
#'   target = 0.30, target_resp = 0.60,
#'   p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
#'   p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
#'   n_cohort = 12, cohort_size = 3, n_trials = 100, seed = 1
#' )
#' print(trials)
#'
#' @seealso \code{\link{boin_et_simulate}}
#'
#' @export
print.boin_et_trials <- function(x, ...) {

  set <- x$settings
  n_doses <- length(set$p_true)

  design <- if (is.null(x$design)) "BOIN-ET" else x$design
  lab <- et_oc_labels(x$design)
  cat("Simulated ", design, " trials\n", sep = "")
  cat("  trials               : ", set$n_trials, "\n", sep = "")
  cat("  doses                : ", n_doses, "\n", sep = "")
  cat(lab$target, format(set$target), "\n", sep = "")
  cat(lab$target_resp, format(set$target_resp), "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")
  cat("\n")

  avg <- rbind(colMeans(x$n_pts), colMeans(x$n_tox), colMeans(x$n_resp))
  rownames(avg) <- c("Patients treated", lab$trial_tox, lab$trial_resp)
  colnames(avg) <- paste0("DL", seq_len(n_doses))
  cat("Average per trial:\n")
  print(round(avg, 2))

  cat("\nStopping reason (%):\n")
  print(round(100 * table(x$stop_reason) / set$n_trials, 1))

  if (!is.null(x$duration)) {
    cat("\nAverage trial duration: ", format(round(mean(x$duration), 1)), "\n",
        sep = "")
  }

  invisible(x)
}
