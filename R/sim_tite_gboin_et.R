#' Operating Characteristics of the TITE-gBOIN-ET Design
#'
#' @description
#'   Simulate trials of the TITE-gBOIN-ET design of Takeda, Yamaguchi,
#'   Taguri and Morita (2023) under one scenario of probabilities of toxicity
#'   and efficacy categories, select the optimal biological dose (OBD) of
#'   each trial, and summarize how often each dose is selected, how many
#'   patients are treated at each dose, with their total scores, and how long
#'   the trials last.
#'
#' @inheritParams tite_gboin_et_simulate
#' @inheritParams gboin_et_select_obd
#' @inheritParams sim_gboin_et
#'
#' @return
#'   An object of class \code{boin_et_oc}, as returned by
#'   \code{\link{sim_gboin_et}}, with \code{design} set to
#'   \code{"TITE-gBOIN-ET"} and one more component
#'   \item{duration}{Average trial duration, in the unit of \code{window}.}
#'
#' @details
#'   The trials are simulated by \code{\link{tite_gboin_et_simulate}} and the
#'   OBD is selected by \code{\link{gboin_et_select_obd}} from the data at
#'   the last decision, as in \code{boinet::tite.gboinet()}; see there for
#'   the rules. The result prints like that of \code{\link{sim_gboin_et}} and
#'   enters \code{\link{oc_results}} and \code{\link{oc_compare}} under the
#'   name TITE-gBOIN-ET.
#'
#'   Takeda et al. (2023) select, among the doses up to the MTD, the dose
#'   with the highest observed mean efficacy score, which is the default
#'   (\code{est_method = "observed"}, \code{obd_method = "max_resp"}), and
#'   exclude a dose for its efficacy with a cutoff of 0.95
#'   (\code{cutoff_eli_resp = 0.95}); the default of 0.99 is that of the
#'   other functions of the BOIN-ET family and of 'boinet'.
#'   \code{inst/reproduce/reproduce_Takeda2023.R} compares the results with
#'   the tables of the article and lists the two points in which the program
#'   of the authors differs from the design.
#'
#' @references
#'   Takeda, K., Yamaguchi, Y., Taguri, M. and Morita, S. (2023).
#'   TITE-gBOIN-ET: Time-to-event generalized Bayesian optimal interval
#'   design to accelerate dose-finding accounting for ordinal graded efficacy
#'   and toxicity outcomes. Biometrical Journal, 65(7), 2200265.
#'   \doi{10.1002/bimj.202200265}
#'
#' @examples
#' # Toxicity grades 0 or 1, 2, 3 and 4, and progressive disease, stable
#' # disease, partial response and complete response, at four doses, with a
#' # 30-day toxicity window, a 60-day efficacy window and one patient every
#' # 10 days on average
#' p_grade <- rbind(
#'   c(0.94, 0.79, 0.62, 0.50),
#'   c(0.05, 0.15, 0.20, 0.20),
#'   c(0.01, 0.05, 0.15, 0.25),
#'   c(0.00, 0.01, 0.03, 0.05)
#' )
#' p_resp_grade <- rbind(
#'   c(0.64, 0.45, 0.20, 0.05),
#'   c(0.30, 0.40, 0.40, 0.15),
#'   c(0.05, 0.10, 0.20, 0.35),
#'   c(0.01, 0.05, 0.20, 0.45)
#' )
#' oc <- sim_tite_gboin_et(
#'   target = 0.47 / 1.5, target_resp = 1.75 / 3,
#'   p_grade = p_grade, p_resp_grade = p_resp_grade,
#'   score = c(0, 0.5, 1, 1.5), score_resp = c(0, 0.25, 1, 3),
#'   n_cohort = 12, cohort_size = 3, window = 30, resp_window = 60,
#'   accrual_rate = 1 / 10, n_trials = 500, seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_tite_gboin_et_multi}},
#'   \code{\link{tite_gboin_et_simulate}},
#'   \code{\link{gboin_et_select_obd}}, \code{\link{sim_gboin_et}},
#'   \code{\link{sim_tite_boin_et}}
#'
#' @export
sim_tite_gboin_et <- function(target, target_resp, p_grade, p_resp_grade,
                              score, score_resp, n_cohort, cohort_size,
                              window, accrual_rate, resp_window = window,
                              accrual = c("uniform", "exponential", "fixed"),
                              event_time = c("weibull", "uniform"),
                              late_fraction = 0.5, resp_late_fraction = 0.5,
                              resp_cor = 0.2, min_completed_share = 0.5,
                              n_trials = 10000, start_dose = 1,
                              n_earlystop = NULL, p_saf = NULL, p_tox = NULL,
                              p_saf_resp = NULL, lambda1 = NULL,
                              lambda2 = NULL, eta1 = NULL, cutoff_eli = 0.95,
                              cutoff_eli_resp = 0.99, p_eli = NULL,
                              p_eli_resp = NULL, boinet_compatible = FALSE,
                              prior = rep(1 / 6, 6),
                              est_method = c("observed", "fp_linear",
                                             "fp_logistic"),
                              obd_method = c("max_resp", "utility_weighted",
                                             "utility_truncated_linear",
                                             "utility_scoring"),
                              min_obd_sample = 1, w1 = 0.33, w2 = 1.09,
                              t_low = NULL, t_upp = NULL, e_low = NULL,
                              e_upp = NULL, psi00 = 40, psi11 = 60,
                              keep_trials = FALSE, verbose = FALSE,
                              seed = 123) {

  call_expr <- match.call()
  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  est_method <- match.arg(est_method)
  obd_method <- match.arg(obd_method)
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- tite_gboin_et_simulate(
    target = target, target_resp = target_resp, p_grade = p_grade,
    p_resp_grade = p_resp_grade, score = score, score_resp = score_resp,
    n_cohort = n_cohort, cohort_size = cohort_size, window = window,
    accrual_rate = accrual_rate, resp_window = resp_window,
    accrual = accrual, event_time = event_time,
    late_fraction = late_fraction, resp_late_fraction = resp_late_fraction,
    resp_cor = resp_cor, min_completed_share = min_completed_share,
    n_trials = n_trials, start_dose = start_dose, n_earlystop = n_earlystop,
    p_saf = p_saf, p_tox = p_tox, p_saf_resp = p_saf_resp, lambda1 = lambda1,
    lambda2 = lambda2, eta1 = eta1, cutoff_eli = cutoff_eli,
    cutoff_eli_resp = cutoff_eli_resp, p_eli = p_eli,
    p_eli_resp = p_eli_resp, boinet_compatible = boinet_compatible,
    prior = prior, seed = seed
  )

  if (verbose) message("Selecting the OBD ...")

  selection <- gboin_et_select_obd(
    n_pts = trials$n_pts, n_tox = trials$n_tox, n_resp = trials$n_resp,
    target = target, target_resp = target_resp, p_saf = p_saf,
    p_tox = p_tox, p_saf_resp = p_saf_resp, cutoff_eli = cutoff_eli,
    cutoff_eli_resp = cutoff_eli_resp, p_eli = p_eli, p_eli_resp = p_eli_resp,
    est_method = est_method, obd_method = obd_method,
    min_obd_sample = min_obd_sample, w1 = w1, w2 = w2, t_low = t_low,
    t_upp = t_upp, e_low = e_low, e_upp = e_upp, psi00 = psi00,
    psi11 = psi11, early_stop = trials$early_stop,
    boinet_compatible = boinet_compatible
  )

  out <- boin_et_summarize(
    trials, selection, est_method, obd_method,
    rule = list(min_obd_sample = min_obd_sample, w1 = w1, w2 = w2,
                t_low = t_low, t_upp = t_upp, e_low = e_low, e_upp = e_upp,
                psi00 = psi00, psi11 = psi11),
    keep_trials = keep_trials, call_expr = call_expr
  )
  out$n_grade_dose <- t(apply(trials$n_grade, c(2L, 3L), mean))
  out$n_resp_grade_dose <- t(apply(trials$n_resp_grade, c(2L, 3L), mean))
  out
}
