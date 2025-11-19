#' Generate Trial Stopping Rule Table
#'
#' @description
#'   Create a lookup table to determine whether the entire trial should be stopped
#'   if excessive toxicity is observed at the lowest dose. Uses vectorized computation
#'   for efficiency. This table is applied only at the lowest dose level.
#'
#' @param target
#'   Numeric. The target toxicity probability (e.g., 0.30 for 30%).
#'
#' @param max_sample_size
#'   Numeric. Maximum sample size (number of patients) for table columns.
#'   Typically 18-30 for phase I trials.
#'
#' @param cutoff_stop
#'   Numeric. Cutoff probability for trial stopping. Default is 0.90.
#'   If Pr(p > target | data) > cutoff_stop at lowest dose, the trial is stopped.
#'
#' @return
#'   A character matrix with stopping decisions:
#'   - Rows represent cumulative number of DLTs at lowest dose (0 to max_sample_size)
#'   - Columns represent cumulative number of patients at lowest dose (1 to max_sample_size)
#'   - Cell values are "STOP" or "GO"
#'   - NA values for logically impossible combinations (n_tox > n_pts) or n_pts < 3
#'
#' @details
#'   The trial stopping rule is based on Bayesian monitoring of safety at the lowest
#'   dose level. The rationale is that if even the lowest dose shows excessive toxicity
#'   with high posterior probability, all doses in the trial may be too toxic, warranting
#'   early trial termination.
#'
#'   The posterior probability Pr(p > target | data) is computed using Beta-Binomial
#'   conjugate prior with uniform prior (Beta(1,1)). The stopping rule is not evaluated
#'   until at least 3 patients have been treated at the lowest dose.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' # Generate stopping rule table for 30% target toxicity rate
#' STOP_DL1 <- get_boin_stopping_boundaries(
#'   target = 0.30,
#'   max_sample_size = 18,
#'   cutoff_stop = 0.90
#' )
#' # Look up decision: 2 DLTs out of 6 patients at lowest dose
#' stop_decision <- STOP_DL1[3, 6]  # Index [n_tox + 1, n_pts]
#' print(stop_decision)
#'
#' @export
get_boin_stopping_boundaries <- function(target, max_sample_size, cutoff_stop) {

  # Initialize stopping rule table
  stop_trial_table <- matrix(NA_character_, nrow = max_sample_size + 1, ncol = max_sample_size)
  colnames(stop_trial_table) <- 1:max_sample_size
  rownames(stop_trial_table) <- 0:max_sample_size

  # Pre-compute all valid (n_tox, n_pts) pairs with condition n_pts >= 3
  valid_indices <- which(
    outer(0:max_sample_size, 1:max_sample_size,
          FUN = function(x, y) (x <= y) & (y >= 3)),
    arr.ind = TRUE
  )

  # Extract indices
  i_indices <- valid_indices[, 1] - 1   # Convert to 0-based DLT counts
  j_indices <- valid_indices[, 2]       # 1-based patient counts

  # Vectorized computation of beta probabilities
  alpha_shape <- i_indices + 1
  beta_shape <- j_indices - i_indices + 1

  prob_exceed <- 1 - pbeta(target, alpha_shape, beta_shape)

  # Apply stopping rule vectorized
  stop_decisions <- ifelse(prob_exceed > cutoff_stop, "STOP", "GO")

  # Fill matrix with vectorized results
  stop_trial_table[valid_indices] <- stop_decisions

  return(stop_trial_table)
}
