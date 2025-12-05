#' Generate Dose Decision Table for BOIN Design
#'
#' @description
#'   Create a pre-computed lookup table that maps (number of DLTs, number of patients)
#'   pairs to dose decisions (Escalate, De-escalate, Stay, or Eliminate). This table
#'   is generated once before the trial and consulted repeatedly during dose assignment.
#'
#' @param target
#'   Numeric. The target toxicity probability (e.g., 0.30 for 30%).
#'
#' @param lambda_e
#'   Numeric. Escalation boundary from `get_boin_boundary()`.
#'   Doses with observed toxicity rate <= lambda_e trigger escalation.
#'
#' @param lambda_d
#'   Numeric. De-escalation boundary from `get_boin_boundary()`.
#'   Doses with observed toxicity rate >= lambda_d trigger de-escalation.
#'
#' @param max_sample_size
#'   Numeric. Maximum sample size (number of patients) for table columns.
#'   Typically 18-30 for phase I trials.
#'
#' @param cutoff_eli
#'   Numeric. Cutoff probability for dose elimination. Default is 0.95.
#'   If Pr(p > target | data) > cutoff_eli, dose is marked for elimination.
#'
#' @return
#'   A character matrix with dose decisions:
#'   - Rows represent cumulative number of DLTs (0 to max_sample_size)
#'   - Columns represent cumulative number of patients (1 to max_sample_size)
#'   - Cell values are decisions: "E" (Escalate), "D" (De-escalate),
#'     "DE" (De-escalate & Eliminate), "S" (Stay), or NA (logically impossible)
#'
#' @details
#'   Decision rules applied in order:
#'   1. Rule DE: If n >= 3 and Pr(p > target | data) > cutoff_eli, eliminate dose
#'   2. Rule E: If observed toxicity rate <= lambda_e, escalate
#'   3. Rule D: If observed toxicity rate >= lambda_d, de-escalate
#'   4. Rule S: Otherwise, stay at current dose
#'
#'   The posterior probability Pr(p > target | data) is computed using Beta-Binomial
#'   conjugate prior with uniform prior (Beta(1,1)).
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' # Generate decision table for 30% target toxicity rate
#' boin_bound <- get_boin_boundary(target = 0.30)
#' DECISION <- get_boin_decision(
#'   target = 0.30,
#'   lambda_e = boin_bound$lambda_e,
#'   lambda_d = boin_bound$lambda_d,
#'   max_sample_size = 18,
#'   cutoff_eli = 0.95
#' )
#' print(DECISION)
#'
#' @importFrom stats pbeta
#'
#' @export
get_boin_decision <- function(target, lambda_e, lambda_d, max_sample_size, cutoff_eli) {

  # Initialize decision table with NA values
  # Rows: number of DLTs (0 to max_sample_size)
  # Columns: number of patients (1 to max_sample_size)
  decision_table <- matrix(NA_character_, nrow = max_sample_size + 1, ncol = max_sample_size)
  colnames(decision_table) <- 1:max_sample_size
  rownames(decision_table) <- 0:max_sample_size

  # Identify all valid (n_tox, n_pts) pairs where n_tox <= n_pts
  # Invalid pairs (n_tox > n_pts) remain as NA
  valid_indices <- which(
    outer(0:max_sample_size, 1:max_sample_size, FUN = "<="),
    arr.ind = TRUE
  )

  # Extract DLT counts and patient counts from valid indices
  # Adjust for 1-based matrix indexing: DLT counts start at 0
  i_indices <- valid_indices[, 1] - 1  # Convert to 0-based DLT counts
  j_indices <- valid_indices[, 2]      # 1-based patient counts

  # Calculate shape parameters for Beta distribution (posterior)
  # Using uniform prior Beta(1,1), posterior is Beta(n_tox+1, n_not_tox+1)
  alpha_shape <- i_indices + 1
  beta_shape <- j_indices - i_indices + 1

  # Compute posterior probability that true toxicity rate exceeds target
  # Pr(p > target | data) using cumulative distribution function
  prob_exceed <- 1 - pbeta(target, alpha_shape, beta_shape)

  # Calculate observed toxicity rate for each valid pair
  tox_rate <- i_indices / j_indices

  # Initialize decision vector for all valid pairs
  decisions <- character(length(i_indices))

  # Apply decision rules in hierarchical order

  # Rule DE: Eliminate dose if n >= 3 and posterior probability exceeds threshold
  elim_rule <- (j_indices >= 3) & (prob_exceed > cutoff_eli)
  decisions[elim_rule] <- "DE"

  # Rule E: Escalate dose if toxicity rate is below escalation boundary
  # Only apply if not already eliminated
  esc_rule <- !elim_rule & (tox_rate <= lambda_e)
  decisions[esc_rule] <- "E"

  # Rule D: De-escalate dose if toxicity rate is above de-escalation boundary
  # Only apply if not already eliminated or escalated
  deesc_rule <- !elim_rule & !esc_rule & (tox_rate >= lambda_d)
  decisions[deesc_rule] <- "D"

  # Rule S: Stay at current dose (default for remaining cases)
  # Applies when toxicity rate is within acceptable interval
  stay_rule <- !elim_rule & !esc_rule & !deesc_rule
  decisions[stay_rule] <- "S"

  # Fill decision table with computed decisions
  decision_table[valid_indices] <- decisions

  return(decision_table)
}
