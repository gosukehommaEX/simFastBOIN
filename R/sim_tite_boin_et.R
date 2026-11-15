#' Operating Characteristics of the TITE-BOIN-ET Design
#'
#' @description
#'   Simulate trials of the TITE-BOIN-ET design of Takeda, Morita and Taguri
#'   (2020) under one scenario of DLT and response probabilities, select the
#'   optimal biological dose (OBD) of each trial, and summarize how often each
#'   dose is selected, how many patients, DLTs and responses occur at each
#'   dose, and how long the trials last.
#'
#' @inheritParams tite_boin_et_simulate
#' @inheritParams boin_et_select_obd
#' @inheritParams sim_boin_et
#'
#' @return
#'   An object of class \code{boin_et_oc}, as returned by
#'   \code{\link{sim_boin_et}}, with \code{design} set to
#'   \code{"TITE-BOIN-ET"} and one more component
#'   \item{duration}{Average trial duration, in the unit of \code{window}.}
#'
#' @details
#'   The trials are simulated by \code{\link{tite_boin_et_simulate}} and the
#'   OBD is selected by \code{\link{boin_et_select_obd}} from the data at the
#'   last decision, as in \code{boinet::tite.boinet()}.
#'
#' @references
#'   Takeda, K., Morita, S. and Taguri, M. (2020). TITE-BOIN-ET: Time-to-event
#'   Bayesian optimal interval design to accelerate dose-finding based on both
#'   efficacy and toxicity outcomes. Pharmaceutical Statistics, 19,
#'   335-349. \doi{10.1002/pst.1995}
#'
#' @examples
#' # Scenario 1 of Takeda et al. (2018) with a 30-day DLT window, a 45-day
#' # response window and one patient every 10 days on average
#' oc <- sim_tite_boin_et(
#'   target = 0.30, target_resp = 0.60,
#'   p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
#'   p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
#'   n_cohort = 12, cohort_size = 3, window = 30, resp_window = 45,
#'   accrual_rate = 1 / 10, n_trials = 500, seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{tite_boin_et_simulate}},
#'   \code{\link{sim_tite_boin_et_multi}}, \code{\link{sim_boin_et}}
#' @export
sim_tite_boin_et <- function(target, target_resp, p_true, p_resp, n_cohort,
                             cohort_size, window, accrual_rate,
                             resp_window = window,
                             accrual = c("uniform", "exponential", "fixed"),
                             event_time = c("weibull", "uniform"),
                             late_fraction = 0.5, resp_late_fraction = 0.5,
                             resp_cor = 0.2, min_completed_share = 0.5,
                             n_trials = 10000, start_dose = 1,
                             n_earlystop = NULL, p_saf = NULL, p_tox = NULL,
                             p_saf_resp = NULL, cutoff_eli = 0.95,
                             cutoff_eli_resp = 0.99, p_eli = NULL,
                             p_eli_resp = NULL, boinet_compatible = FALSE,
                             prior = rep(1 / 6, 6),
                             est_method = c("observed", "multi_iso",
                                            "fp_logistic"),
                             obd_method = c("max_resp", "utility_weighted",
                                            "utility_truncated_linear",
                                            "utility_scoring"),
                             w1 = 0.33, w2 = 1.09, t_low = NULL, t_upp = NULL,
                             e_low = NULL, e_upp = NULL, psi00 = 40,
                             psi11 = 60, keep_trials = FALSE, verbose = FALSE,
                             seed = 123) {

  call_expr <- match.call()
  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  est_method <- match.arg(est_method)
  obd_method <- match.arg(obd_method)
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- tite_boin_et_simulate(
    target = target, target_resp = target_resp, p_true = p_true,
    p_resp = p_resp, n_cohort = n_cohort, cohort_size = cohort_size,
    window = window, accrual_rate = accrual_rate, resp_window = resp_window,
    accrual = accrual, event_time = event_time, late_fraction = late_fraction,
    resp_late_fraction = resp_late_fraction, resp_cor = resp_cor,
    min_completed_share = min_completed_share,
    n_trials = n_trials, start_dose = start_dose, n_earlystop = n_earlystop,
    p_saf = p_saf, p_tox = p_tox, p_saf_resp = p_saf_resp,
    cutoff_eli = cutoff_eli, cutoff_eli_resp = cutoff_eli_resp,
    p_eli = p_eli, p_eli_resp = p_eli_resp,
    boinet_compatible = boinet_compatible, prior = prior, seed = seed
  )

  if (verbose) message("Selecting the OBD ...")

  selection <- boin_et_select_obd(
    n_pts = trials$n_pts, n_tox = trials$n_tox, n_resp = trials$n_resp,
    target = target, target_resp = target_resp, p_saf = p_saf,
    p_tox = p_tox, p_saf_resp = p_saf_resp, cutoff_eli = cutoff_eli,
    cutoff_eli_resp = cutoff_eli_resp, p_eli = p_eli, p_eli_resp = p_eli_resp,
    est_method = est_method, obd_method = obd_method, w1 = w1, w2 = w2,
    t_low = t_low, t_upp = t_upp, e_low = e_low, e_upp = e_upp,
    psi00 = psi00, psi11 = psi11, early_stop = trials$early_stop,
    boinet_compatible = boinet_compatible
  )

  boin_et_summarize(
    trials, selection, est_method, obd_method,
    rule = list(w1 = w1, w2 = w2, t_low = t_low, t_upp = t_upp, e_low = e_low,
                e_upp = e_upp, psi00 = psi00, psi11 = psi11),
    keep_trials = keep_trials, call_expr = call_expr
  )
}
