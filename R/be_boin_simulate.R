#' Simulate BE-BOIN Trials
#'
#' @description
#'   Run the BOIN design with backfilling and late-onset toxicity (BE-BOIN) of
#'   Chen et al. (2026) over many simulated trials. Dose escalation follows
#'   TITE-BOIN, which decides while some patients are still pending; patients
#'   who arrive while the escalation is suspended are backfilled to a lower dose
#'   that is safe and has shown a response. Every estimate imputes the pending
#'   patients, backfilled ones included. Return the patient, DLT, backfill and
#'   response counts at every dose together with the duration of each trial.
#'
#' @inheritParams bf_boin_simulate
#'
#' @param conflict_dose
#'   Character scalar. When the data of several backfilled doses conflict with
#'   those of the current dose, which of them anchors the pooled estimate:
#'   \code{"lowest"} (the default here, following Chen et al., 2026) or
#'   \code{"highest"} (the default of \code{\link{bf_boin_simulate}}, following
#'   Zhao et al., 2024).
#'
#' @param max_pending_ratio
#'   Numeric scalar greater than 0 and at most 1. Rule 1 of the design: the
#'   escalation is suspended when the proportion of pending patients at the
#'   current dose exceeds this value, unless the dose is de-escalated whatever
#'   their outcomes. Defaults to 0.49, which suspends when fewer than 51 percent
#'   have completed the assessment, as in Chen et al. (2026). See
#'   \code{\link{tite_boin_decision_table}}.
#'
#' @param min_completed
#'   Integer scalar. While some patients are pending, escalation requires at
#'   least this many patients with a completed assessment. Defaults to 0, which
#'   the design does not use. See \code{\link{tite_boin_decision_table}}.
#'
#' @param min_follow_up
#'   Numeric scalar between 0 and 1. Rule 2 of the design: escalation requires
#'   every pending patient at the current dose to have been followed for at
#'   least this fraction of the window. Defaults to 0.25, as in Chen et al.
#'   (2026). See \code{\link{tite_boin_decision_table}}.
#'
#' @return
#'   An object of class \code{c("be_boin_trials", "backfill_trials")}, with the
#'   components described in \code{\link{bf_boin_simulate}}.
#'
#' @details
#'   \strong{Dose escalation.} The dose for a cohort is decided when its first
#'   patient arrives, with the rules of \code{\link{tite_boin_decision_table}}
#'   for \code{method = "imputation"} and the suspension rules
#'   \code{max_pending_ratio} and \code{min_follow_up}; every patient at the
#'   current dose counts, backfilled ones included. A suspension lasts until a
#'   decision can be taken again, as in \code{\link{tite_boin_simulate}}.
#'   Elimination counts every treated patient, the pending ones as without DLT,
#'   and is checked at every dose.
#'
#'   \strong{Estimates.} The DLT rate at a dose is the single mean imputation
#'   estimate of TITE-BOIN, equation (1) of Chen et al. (2026), and the rate
#'   pooled over several doses adds up the observed and imputed DLTs of each
#'   dose, each imputed with its own posterior mean, and divides by the number of
#'   treated patients, equation (2). Both enter the opening and closing of
#'   backfill doses and the resolution of conflicts, which otherwise follow
#'   \code{\link{bf_boin_simulate}}. The categories of a backfilled dose for
#'   the conflicts compare its estimate with the escalation and de-escalation
#'   boundaries directly.
#'
#'   \strong{End of the trial, arrivals and random numbers} are as in
#'   \code{\link{bf_boin_simulate}}. When no dose is ever opened for
#'   backfilling, for example with every \code{p_resp} equal to 0, and with
#'   \code{no_slot = "wait"}, the trials are identical to those of
#'   \code{\link{tite_boin_simulate}} with \code{method = "imputation"}, the
#'   same suspension rules and the same seed.
#'
#'   Chen et al. (2026) report their operating characteristics for scenarios
#'   that are given only as figures, so they cannot be reproduced.
#'
#' @references
#'   Chen, K., Zhao, Y., Takeda, K. and Yuan, Y. (2026). BE-BOIN: A Dose
#'   Optimization Design Accommodating Backfill and Late-Onset Toxicity.
#'   Therapeutic Innovation and Regulatory Science.
#'   \doi{10.1007/s43441-026-00994-0}
#'
#'   Chen, K., Chen, T.-Y., Zhang, Y., Lin, R. and Yuan, Y. (2025). Practical
#'   Considerations for Using the TITE-BOIN Design to Handle Late-Onset Toxicity
#'   or Fast Accrual in Phase I Trials. Clinical Cancer Research, 31(13),
#'   2573-2580.
#'
#'   Zhao, Y., Yuan, Y., Korn, E. L. and Freidlin, B. (2024). Backfilling
#'   Patients in Phase I Dose-Escalation Trials Using Bayesian Optimal Interval
#'   Design (BOIN). Clinical Cancer Research, 30(4), 673-679.
#'
#' @examples
#' # A three month window with two patients a month
#' trials <- be_boin_simulate(
#'   target = 0.25,
#'   p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
#'   p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   window = 3,
#'   accrual_rate = 2,
#'   n_earlystop = 9,
#'   n_trials = 200,
#'   seed = 123
#' )
#' trials
#'
#' summary(trials$duration)
#'
#' @seealso \code{\link{bf_boin_simulate}}, \code{\link{tite_boin_simulate}}
#'
#' @export
be_boin_simulate <- function(target, p_true, p_resp, n_cohort, cohort_size,
                             window, accrual_rate, n_cap = 12,
                             backfill_dose = c("highest", "lowest"),
                             conflict_dose = c("lowest", "highest"),
                             no_slot = c("wait", "leave"),
                             accrual = c("exponential", "uniform", "fixed"),
                             dlt_time = c("weibull", "uniform"),
                             late_fraction = 0.5, resp_window = window,
                             resp_late_fraction = 0.5, resp_cor = 0,
                             max_pending_ratio = 0.49, min_completed = 0,
                             min_follow_up = 0.25,
                             n_trials = 10000, start_dose = 1, n_earlystop = 18,
                             p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                             extrasafe = FALSE, offset = 0.05,
                             n_earlystop_rule = c("with_stay", "simple"),
                             seed = 123) {

  backfill_simulate(
    design = "be_boin", target = target, p_true = p_true, p_resp = p_resp,
    n_cohort = n_cohort, cohort_size = cohort_size, window = window,
    accrual_rate = accrual_rate, n_cap = n_cap,
    backfill_dose = match.arg(backfill_dose),
    conflict_dose = match.arg(conflict_dose), no_slot = match.arg(no_slot),
    accrual = match.arg(accrual), dlt_time = match.arg(dlt_time),
    late_fraction = late_fraction, resp_window = resp_window,
    resp_late_fraction = resp_late_fraction, resp_cor = resp_cor,
    max_pending_ratio = max_pending_ratio, min_completed = min_completed,
    min_follow_up = min_follow_up, n_trials = n_trials,
    start_dose = start_dose, n_earlystop = n_earlystop, p_saf = p_saf,
    p_tox = p_tox, cutoff_eli = cutoff_eli, extrasafe = extrasafe,
    offset = offset, stay_on_1_of_3 = FALSE,
    n_earlystop_rule = match.arg(n_earlystop_rule), seed = seed
  )
}
