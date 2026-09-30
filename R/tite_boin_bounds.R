#' Decision Rule of the TITE-BOIN Design at Given States of the Current Dose
#'
#' @description
#'   Internal workhorse of \code{\link{tite_boin_decision_table}}. For each state
#'   of the current dose, given by the number of patients treated, the number of
#'   DLTs observed and the number of patients whose DLT assessment is pending,
#'   return the decision and the values of the follow-up statistic at which the
#'   decision changes.
#'
#' @param n
#'   Integer vector. Number of patients treated at the current dose.
#'
#' @param n_tox
#'   Integer vector. Number of DLTs observed, all among the patients who
#'   completed the assessment.
#'
#' @param n_pending
#'   Integer vector. Number of patients whose assessment is pending.
#'
#' @param bound
#'   An object of class \code{boin_boundary} that covers \code{max(n)}.
#'
#' @param method
#'   Character scalar, \code{"imputation"} or \code{"ess"}.
#'
#' @param max_pending_ratio
#'   Numeric scalar. Accrual is suspended when the proportion of pending patients
#'   exceeds this value.
#'
#' @param min_completed
#'   Integer scalar. Escalation requires at least this many patients with a
#'   completed assessment while some patients are pending.
#'
#' @return
#'   A data frame with columns \code{n}, \code{n_tox}, \code{n_pending},
#'   \code{decision}, \code{esc_bound} and \code{deesc_bound}, as described in
#'   \code{\link{tite_boin_decision_table}}.
#'
#' @noRd
tite_boin_bounds <- function(n, n_tox, n_pending, bound, method,
                             max_pending_ratio, min_completed) {

  target <- bound$target
  lambda_e <- bound$lambda_e
  lambda_d <- bound$lambda_d
  n_done <- n - n_pending
  m <- length(n)

  decision <- character(m)
  esc_bound <- rep(NA_real_, m)
  deesc_bound <- rep(NA_real_, m)

  # No pending patient: the BOIN rule with its integer boundaries, so that these
  # states agree with boin_decision_table() and with sim_boin().
  complete <- n_pending == 0L
  decision[complete] <- ifelse(
    n_tox[complete] <= bound$b_esc[n[complete]], "E",
    ifelse(n_tox[complete] >= bound$b_deesc[n[complete]], "D", "S")
  )

  # Pending patients. The follow-up statistic ranges over [lower, upper). The
  # dose is escalated when the statistic is at least esc and de-escalated when
  # it is at most deesc.
  if (method == "imputation") {
    # Single mean imputation, Yuan et al. (2018), Supplementary Appendix A,
    # equations (4) and (5). The prior Beta(alpha, 1 - alpha) with
    # alpha = target / 2 has an effective sample size of one, so the posterior
    # mean is (n_tox + alpha) / (n_done + 1).
    alpha <- 0.5 * target
    p_post <- (n_tox + alpha) / (n_done + 1)
    odds <- (1 - p_post) / p_post
    # Long-memory coherence: escalation only while the observed rate is below
    # the target, de-escalation once it reaches the target. Allowing
    # de-escalation at equality reproduces Table 1 of Yuan et al. (2018).
    at_or_above <- n_tox / n >= target - 1e-12
    esc <- ifelse(at_or_above, Inf, n_pending - odds * (n * lambda_e - n_tox))
    deesc <- ifelse(at_or_above, n_pending - odds * (n * lambda_d - n_tox), -Inf)
    lower <- rep(0, m)
    upper <- as.numeric(n_pending)
  } else {
    # Effective sample size, Lin and Yuan (2020): the estimate is n_tox / ESS
    # with ESS = n_done + STFT. Without any DLT the estimate is zero, which
    # never de-escalates.
    esc <- n_tox / lambda_e
    deesc <- ifelse(n_tox > 0L, n_tox / lambda_d, -Inf)
    lower <- as.numeric(n_done)
    upper <- as.numeric(n)
  }

  suspend <- n_pending / n > max_pending_ratio
  esc_action <- ifelse(n_done >= min_completed, "E", "SUS")
  deesc_always <- deesc >= upper
  esc_always <- esc <= lower
  has_esc <- esc < upper
  has_deesc <- deesc > lower

  mixed <- ifelse(has_esc, paste0(esc_action, "/S"), "S")
  mixed <- ifelse(has_deesc, paste0(mixed, "/D"), mixed)
  pending_decision <- ifelse(
    deesc_always, "D",
    ifelse(suspend, "SUS", ifelse(esc_always, esc_action, mixed))
  )

  pend <- !complete
  decision[pend] <- pending_decision[pend]
  open <- pend & !deesc_always & !suspend & !esc_always
  esc_bound[open & has_esc] <- esc[open & has_esc]
  deesc_bound[open & has_deesc] <- deesc[open & has_deesc]

  # Elimination takes precedence over every other rule. It uses all treated
  # patients, the pending ones counted as without DLT.
  b_elim <- bound$b_elim[n]
  eliminate <- !is.na(b_elim) & n_tox >= b_elim
  decision[eliminate] <- "DE"
  esc_bound[eliminate] <- NA_real_
  deesc_bound[eliminate] <- NA_real_

  data.frame(
    n = n,
    n_tox = n_tox,
    n_pending = n_pending,
    decision = decision,
    esc_bound = esc_bound,
    deesc_bound = deesc_bound,
    stringsAsFactors = FALSE
  )
}
