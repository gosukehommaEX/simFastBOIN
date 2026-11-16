#' Print Simulated BOIN12 Trials
#'
#' @description
#'   Display the settings, the average numbers of patients, DLTs and
#'   responses at each dose and the stopping reasons of trials simulated by
#'   \code{\link{boin12_simulate}}.
#'
#' @param x
#'   An object of class \code{boin12_trials}.
#'
#' @param ...
#'   Not used.
#'
#' @return
#'   \code{x}, invisibly.
#'
#' @examples
#' trials <- boin12_simulate(
#'   target = 0.35, target_resp = 0.25,
#'   p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'   p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
#'   n_cohort = 12, cohort_size = 3, n_trials = 100, seed = 1
#' )
#' print(trials)
#'
#' @export
print.boin12_trials <- function(x, ...) {

  set <- x$settings
  n_doses <- length(set$p_true)

  design <- if (is.null(x$design)) "BOIN12" else x$design
  cat("Simulated ", design, " trials\n", sep = "")
  cat("  trials               : ", set$n_trials, "\n", sep = "")
  cat("  doses                : ", n_doses, "\n", sep = "")
  cat("  DLT limit            : ", format(set$target), "\n", sep = "")
  cat("  response limit       : ", format(set$target_resp), "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")
  cat("\n")

  avg <- rbind(
    "Patients treated" = colMeans(x$n_pts),
    "DLTs observed" = colMeans(x$n_tox),
    "Responses observed" = colMeans(x$n_resp)
  )
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
