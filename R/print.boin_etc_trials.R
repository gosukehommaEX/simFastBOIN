#' Print Simulated BOIN-ETC Trials
#'
#' @param x An object of class \code{boin_etc_trials} from
#'   \code{\link{boin_etc_simulate}}.
#'
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' p_true <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
#' p_resp <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
#' trials <- boin_etc_simulate(target = 0.30, target_resp = 0.60,
#'                             p_true = p_true, p_resp = p_resp,
#'                             n_cohort = 8, cohort_size = 3,
#'                             n_trials = 100, seed = 1)
#' print(trials)
#'
#' @seealso \code{\link{boin_etc_simulate}}
#'
#' @export
print.boin_etc_trials <- function(x, ...) {

  set <- x$settings
  n_rows <- nrow(set$p_true)
  n_cols <- ncol(set$p_true)

  cat("Simulated ", x$design, " trials\n", sep = "")
  cat("  trials               : ", set$n_trials, "\n", sep = "")
  cat("  dose levels          : ", n_rows, " of drug A, ", n_cols,
      " of drug B\n", sep = "")
  cat("  target DLT rate      : ", format(set$target), "\n", sep = "")
  cat("  target response rate : ", format(set$target_resp), "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")
  cat("  subtrial end (n_stop): ", set$n_earlystop, "\n", sep = "")

  labels <- list(paste0("A", seq_len(n_rows)), paste0("B", seq_len(n_cols)))
  for (item in c("n_pts", "n_tox", "n_resp")) {
    title <- c(n_pts = "Patients treated", n_tox = "DLTs observed",
               n_resp = "Responses observed")[[item]]
    avg <- matrix(apply(array(as.numeric(x[[item]]), dim(x[[item]])),
                        c(2L, 3L), mean), n_rows, n_cols, dimnames = labels)
    cat("\n", title, ", average per trial:\n", sep = "")
    print(round(avg, 2))
  }

  cat("\nAverage number of subtrials: ", format(round(mean(x$n_subtrials), 2)),
      "\n", sep = "")
  cat("\nStopping reason (%):\n")
  print(round(100 * table(x$stop_reason) / set$n_trials, 1))

  invisible(x)
}
