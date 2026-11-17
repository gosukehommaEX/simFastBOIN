#' Operating Characteristics of the TITE-BOIN12 Design
#'
#' @description
#'   Simulate trials of the TITE-BOIN12 design of Zhou et al. (2022), in which
#'   the DLT and the response may be observed late, under one scenario of DLT
#'   and response probabilities, select the optimal biological dose (OBD) of
#'   each trial, and summarize how often each dose is selected, how many
#'   patients, DLTs and responses occur at each dose and how long the trials
#'   last.
#'
#' @inheritParams tite_boin12_simulate
#' @inheritParams sim_boin12
#'
#' @return
#'   An object of class \code{boin12_oc}, as returned by
#'   \code{\link{sim_boin12}}, with \code{design} set to \code{"TITE-BOIN12"}
#'   and one more component
#'   \item{duration}{Average trial duration, in the unit of \code{window}.}
#'
#' @details
#'   The trials are simulated by \code{\link{tite_boin12_simulate}}, and the
#'   OBD of each trial is selected by \code{\link{boin12_select_obd}} from the
#'   complete data, after every patient has finished both assessments; see
#'   there for the rules.
#'
#' @references
#'   Zhou, Y., Lin, R., Lee, J. J., Li, D., Wang, L., Li, R. and Yuan, Y.
#'   (2022). TITE-BOIN12: A Bayesian phase I/II trial design to find the
#'   optimal biological dose with late-onset toxicity and efficacy. Statistics
#'   in Medicine, 41(11), 1918-1931. \doi{10.1002/sim.9337}
#'
#' @examples
#' # Scenario 1 of Zhou et al. (2022): windows of 45 and 60 days, three
#' # patients a month, 36 patients in cohorts of 3
#' oc <- sim_tite_boin12(
#'   target = 0.35, target_resp = 0.25,
#'   p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'   p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
#'   n_cohort = 12, cohort_size = 3, window = 45, resp_window = 60,
#'   accrual_rate = 3 / 30, n_trials = 500, seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_tite_boin12_multi}},
#'   \code{\link{tite_boin12_simulate}}, \code{\link{sim_boin12}}
#'
#' @export
sim_tite_boin12 <- function(target, target_resp, p_true, p_resp, n_cohort,
                            cohort_size, window, accrual_rate,
                            resp_window = window,
                            accrual = c("uniform", "exponential", "fixed"),
                            event_time = c("weibull", "uniform"),
                            late_fraction = 0.5, resp_late_fraction = 0.5,
                            min_completed_share = 0.5,
                            pending_estimate = c("conditional", "marginal"),
                            n_trials = 10000, start_dose = 1,
                            n_earlystop = NULL,
                            n_earlystop_rule = c("with_stay", "simple"),
                            p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                            cutoff_eli_resp = 0.90, psi00 = 40, psi11 = 60,
                            benchmark_weight = 0.5, utility_prior = c(1, 1),
                            n_star = 6, n_explore = 9,
                            tie_break = c("higher", "random"),
                            no_candidate = c("nearest", "stop"), resp_cor = 0,
                            resp_cor_type = c("copula", "binary"),
                            est_method = c("posterior_mean", "observed",
                                           "tail_probability"),
                            keep_trials = FALSE, verbose = FALSE,
                            seed = 123) {

  call_expr <- match.call()
  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  pending_estimate <- match.arg(pending_estimate)
  n_earlystop_rule <- match.arg(n_earlystop_rule)
  tie_break <- match.arg(tie_break)
  no_candidate <- match.arg(no_candidate)
  resp_cor_type <- match.arg(resp_cor_type)
  est_method <- match.arg(est_method)
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- tite_boin12_simulate(
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
    psi00 = psi00, psi11 = psi11, benchmark_weight = benchmark_weight,
    utility_prior = utility_prior, n_star = n_star, n_explore = n_explore,
    tie_break = tie_break, no_candidate = no_candidate, resp_cor = resp_cor,
    resp_cor_type = resp_cor_type, seed = seed
  )

  if (verbose) message("Selecting the OBD ...")

  selection <- boin12_select_obd(
    n_pts = trials$n_pts, n_tox = trials$n_tox, n_resp = trials$n_resp,
    n_tox_resp = trials$n_tox_resp, target = target,
    target_resp = target_resp, cutoff_eli = cutoff_eli,
    cutoff_eli_resp = cutoff_eli_resp, psi00 = psi00, psi11 = psi11,
    benchmark_weight = benchmark_weight, utility_prior = utility_prior,
    est_method = est_method, early_stop = trials$early_stop
  )

  boin12_summarize(trials, selection, est_method, keep_trials, call_expr)
}
