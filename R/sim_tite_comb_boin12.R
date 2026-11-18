#' Operating Characteristics of the TITE-Comb-BOIN12 Design
#'
#' @description
#'   Simulate trials of the TITE-Comb-BOIN12 design of Lu et al. (2025), in
#'   which the DLT and the response may be observed late, with
#'   \code{\link{tite_comb_boin12_simulate}}, select the optimal biological
#'   dose combination (OBDC) of every trial with
#'   \code{\link{comb_boin12_select_obd}} from the complete data, and
#'   summarize the selection percentages, the numbers of patients, DLTs and
#'   responses at every combination and the trial durations.
#'
#' @inheritParams tite_comb_boin12_simulate
#' @inheritParams sim_comb_boin12
#'
#' @return
#'   An object of class \code{comb_boin12_oc}, as returned by
#'   \code{\link{sim_comb_boin12}}, with \code{design} set to
#'   \code{"TITE-Comb-BOIN12"} and one more component
#'   \item{duration}{Average trial duration, in the unit of \code{window}.}
#'
#' @details
#'   The OBDC of each trial is selected after every patient has finished
#'   both assessments, by the rules of \code{\link{comb_boin12_select_obd}}.
#'
#' @references
#'   Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
#'   utility-based Bayesian optimal interval design for dose optimization in
#'   cancer drug-combination trials. Statistics in Biopharmaceutical
#'   Research, 17(2), 266-276. \doi{10.1080/19466315.2024.2370403}
#'
#' @examples
#' # Scenario 16 of Lu et al. (2025), with their design settings, windows of
#' # 1 and 2 months and two patients a month
#' p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
#' p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
#' oc <- sim_tite_comb_boin12(
#'   target = 0.33, target_resp = 0.2, p_true = p_true, p_resp = p_resp,
#'   n_cohort = 14, cohort_size = 3, window = 1, resp_window = 2,
#'   accrual_rate = 2, event_time = "uniform", n_trials = 200,
#'   psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11, benchmark_weight = 0.6,
#'   extrasafe = TRUE, offset = 0.10, est_method = "observed", seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_tite_comb_boin12_multi}},
#'   \code{\link{tite_comb_boin12_simulate}}, \code{\link{sim_comb_boin12}}
#'
#' @export
sim_tite_comb_boin12 <- function(target, target_resp, p_true, p_resp,
                                 n_cohort, cohort_size, window, accrual_rate,
                                 resp_window = window,
                                 accrual = c("uniform", "exponential",
                                             "fixed"),
                                 event_time = c("weibull", "uniform"),
                                 late_fraction = 0.5,
                                 resp_late_fraction = 0.5,
                                 min_completed_share = 0.5,
                                 pending_estimate = c("conditional",
                                                      "marginal"),
                                 n_trials = 10000, start_dose = c(1, 1),
                                 n_earlystop = NULL,
                                 n_earlystop_rule = c("with_stay", "simple"),
                                 p_saf = NULL, p_tox = NULL,
                                 cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
                                 extrasafe = FALSE, offset = 0.05,
                                 psi00 = 40, psi11 = 60,
                                 benchmark_weight = 0.5,
                                 utility_prior = c(1, 1), n_star = 6,
                                 tie_break = c("higher", "random"),
                                 no_candidate = c("nearest", "stop"),
                                 resp_cor = 0,
                                 resp_cor_type = c("copula", "binary"),
                                 est_method = c("logistic", "observed"),
                                 obd_bound = c("rectangle", "contour"),
                                 dose_a = NULL, dose_b = NULL,
                                 logistic_scale = c(2.5, 2.5),
                                 n_draws = 4000, keep_trials = FALSE,
                                 verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  pending_estimate <- match.arg(pending_estimate)
  n_earlystop_rule <- match.arg(n_earlystop_rule)
  tie_break <- match.arg(tie_break)
  no_candidate <- match.arg(no_candidate)
  resp_cor_type <- match.arg(resp_cor_type)
  est_method <- match.arg(est_method)
  obd_bound <- match.arg(obd_bound)
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- tite_comb_boin12_simulate(
    target = target, target_resp = target_resp, p_true = p_true,
    p_resp = p_resp, n_cohort = n_cohort, cohort_size = cohort_size,
    window = window, accrual_rate = accrual_rate, resp_window = resp_window,
    accrual = accrual, event_time = event_time, late_fraction = late_fraction,
    resp_late_fraction = resp_late_fraction,
    min_completed_share = min_completed_share,
    pending_estimate = pending_estimate, n_trials = n_trials,
    start_dose = start_dose, n_earlystop = n_earlystop,
    n_earlystop_rule = n_earlystop_rule, p_saf = p_saf, p_tox = p_tox,
    cutoff_eli = cutoff_eli, cutoff_eli_resp = cutoff_eli_resp,
    extrasafe = extrasafe, offset = offset, psi00 = psi00, psi11 = psi11,
    benchmark_weight = benchmark_weight, utility_prior = utility_prior,
    n_star = n_star, tie_break = tie_break, no_candidate = no_candidate,
    resp_cor = resp_cor, resp_cor_type = resp_cor_type, seed = seed
  )

  if (verbose) message("Selecting the OBD combination ...")

  selection <- comb_boin12_select_obd(
    n_pts = trials$n_pts, n_tox = trials$n_tox, n_resp = trials$n_resp,
    target = target, target_resp = target_resp, cutoff_eli = cutoff_eli,
    cutoff_eli_resp = cutoff_eli_resp, extrasafe = extrasafe,
    offset = offset, psi00 = psi00, psi11 = psi11, est_method = est_method,
    obd_bound = obd_bound, dose_a = dose_a, dose_b = dose_b,
    logistic_scale = logistic_scale, n_draws = n_draws,
    early_stop = trials$early_stop
  )

  out <- comb_boin12_summarize(trials, selection, est_method, obd_bound,
                               keep_trials, call_expr)
  out$settings$dose_a <- dose_a
  out$settings$dose_b <- dose_b
  out$settings$logistic_scale <- logistic_scale
  out$settings$n_draws <- n_draws
  out
}
