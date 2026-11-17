#' Print Simulated Comb-BOIN12 Trials
#'
#' @param x An object of class \code{comb_boin12_trials} from
#'   \code{\link{comb_boin12_simulate}}.
#'
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
#' p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
#' trials <- comb_boin12_simulate(target = 0.33, target_resp = 0.2,
#'                                p_true = p_true, p_resp = p_resp,
#'                                n_cohort = 14, cohort_size = 3,
#'                                n_trials = 100, seed = 1)
#' print(trials)
#'
#' @seealso \code{\link{comb_boin12_simulate}}
#'
#' @export
print.comb_boin12_trials <- function(x, ...) {

  set <- x$settings
  n_rows <- nrow(set$p_true)
  n_cols <- ncol(set$p_true)

  design <- if (is.null(x$design)) "Comb-BOIN12" else x$design
  cat("Simulated ", design, " trials\n", sep = "")
  cat("  trials               : ", set$n_trials, "\n", sep = "")
  cat("  dose levels          : ", n_rows, " of drug A, ", n_cols,
      " of drug B\n", sep = "")
  cat("  DLT limit            : ", format(set$target), "\n", sep = "")
  cat("  response limit       : ", format(set$target_resp), "\n", sep = "")
  cat("  cohorts              : ", set$n_cohort, " of size ",
      paste(unique(set$cohort_size), collapse = "/"), "\n", sep = "")
  cat("  max sample size      : ", set$max_total_pts, "\n", sep = "")

  labels <- list(paste0("A", seq_len(n_rows)), paste0("B", seq_len(n_cols)))
  for (item in c("n_pts", "n_tox", "n_resp")) {
    title <- c(n_pts = "Patients treated", n_tox = "DLTs observed",
               n_resp = "Responses observed")[[item]]
    avg <- matrix(apply(array(as.numeric(x[[item]]), dim(x[[item]])),
                        c(2L, 3L), mean), n_rows, n_cols, dimnames = labels)
    cat("\n", title, ", average per trial:\n", sep = "")
    print(round(avg, 2))
  }

  cat("\nStopping reason (%):\n")
  print(round(100 * table(x$stop_reason) / set$n_trials, 1))

  if (!is.null(x$duration)) {
    cat("\nAverage trial duration: ", format(round(mean(x$duration), 1)), "\n",
        sep = "")
  }

  invisible(x)
}
