#' Calculate Effective Sample Size (ESS) for TITE-BOIN
#'
#' @description
#' Computes the Effective Sample Size (ESS) for a dose level in TITE-BOIN design
#' based on the approximated likelihood method proposed by Lin and Yuan (2020).
#' ESS accounts for both completed patients and pending patients with partial
#' follow-up information.
#'
#' @param n_completed Integer. Number of patients who have completed the full
#'   DLT assessment window at the current dose.
#' @param follow_up_times Numeric vector. Follow-up times for pending patients
#'   (i.e., patients who have not yet completed the assessment window).
#'   Time unit should match assessment_window.
#' @param assessment_window Numeric. Length of the DLT assessment window
#'   (e.g., 28 days). Must be positive.
#'
#' @return Numeric. The Effective Sample Size (ESS), which can be non-integer.
#'   ESS represents the "information content" from all patients at the dose level,
#'   accounting for partial information from pending patients.
#'
#' @details
#' The ESS is calculated using the formula from Lin and Yuan (2020):
#'
#' \deqn{ESS = n_{completed} + \sum_{i \in pending} w_i}
#'
#' where \eqn{w_i = \frac{follow\_up\_time_i}{assessment\_window}} for each
#' pending patient \eqn{i}.
#'
#' This approach weights pending patients by their proportional follow-up time,
#' effectively converting partial information into an equivalent number of
#' completed patients. For example:
#' - A patient followed for 14 days in a 28-day window contributes weight 0.5
#' - A patient followed for 21 days in a 28-day window contributes weight 0.75
#'
#' The ESS is then used to calculate the adjusted toxicity rate:
#' \deqn{\tilde{p}_j = \frac{y_j}{ESS_j}}
#'
#' where \eqn{y_j} is the observed number of DLTs at dose \eqn{j}.
#'
#' @references
#' Lin, R., & Yuan, Y. (2020). Time-to-event model-assisted designs for
#' dose-finding trials with delayed toxicity. \emph{Biostatistics}, 21(4), 807-824.
#'
#' @examples
#' # Example 1: Only completed patients (no pending)
#' ess1 <- calculate_ess(
#'   n_completed = 6,
#'   follow_up_times = numeric(0),
#'   assessment_window = 28
#' )
#' ess1  # Returns 6
#'
#' # Example 2: Mix of completed and pending patients
#' # 6 completed patients
#' # 3 pending patients with 14, 21, and 7 days of follow-up
#' ess2 <- calculate_ess(
#'   n_completed = 6,
#'   follow_up_times = c(14, 21, 7),
#'   assessment_window = 28
#' )
#' ess2  # Returns 6 + (14/28) + (21/28) + (7/28) = 6 + 0.5 + 0.75 + 0.25 = 7.5
#'
#' # Example 3: Only pending patients
#' ess3 <- calculate_ess(
#'   n_completed = 0,
#'   follow_up_times = c(7, 14, 21, 28),
#'   assessment_window = 28
#' )
#' ess3  # Returns 0 + (7+14+21+28)/28 = 2.5
#'
#' # Example 4: Application in decision making
#' # At dose level 3, we have observed 2 DLTs
#' n_tox <- 2
#' ess <- calculate_ess(
#'   n_completed = 4,
#'   follow_up_times = c(10, 15, 20),
#'   assessment_window = 28
#' )
#' # ess = 4 + (10+15+20)/28 ≈ 5.607
#'
#' # Calculate adjusted toxicity rate
#' p_tilde <- n_tox / ess
#' p_tilde  # ≈ 0.357
#'
#' # Compare with BOIN boundaries (e.g., lambda_e = 0.197, lambda_d = 0.355)
#' # Decision: p_tilde > lambda_d, so de-escalate
#'
#' @export
calculate_ess <- function(n_completed, follow_up_times, assessment_window) {

  # Input validation
  if (n_completed < 0) {
    stop("n_completed must be non-negative")
  }
  if (!is.numeric(follow_up_times)) {
    stop("follow_up_times must be a numeric vector")
  }
  if (assessment_window <= 0) {
    stop("assessment_window must be positive")
  }
  if (any(follow_up_times < 0)) {
    stop("All follow_up_times must be non-negative")
  }
  if (any(follow_up_times > assessment_window)) {
    warning("Some follow_up_times exceed assessment_window; these will be capped at assessment_window")
    follow_up_times <- pmin(follow_up_times, assessment_window)
  }

  # Calculate weights for pending patients
  # Weight = follow_up_time / assessment_window
  if (length(follow_up_times) == 0) {
    # No pending patients
    pending_contribution <- 0
  } else {
    # Sum of weights from all pending patients
    pending_contribution <- sum(follow_up_times) / assessment_window
  }

  # ESS = completed patients + weighted contribution from pending patients
  ess <- n_completed + pending_contribution

  return(ess)
}
