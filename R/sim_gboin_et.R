#' Operating Characteristics of the gBOIN-ET Design
#'
#' @description
#'   Simulate trials of the gBOIN-ET design of Takeda, Morita and Taguri
#'   (2022) under one scenario of probabilities of toxicity and efficacy
#'   categories, select the optimal biological dose (OBD) of each trial, and
#'   summarize how often each dose is selected and how many patients are
#'   treated at each dose, with their total scores.
#'
#' @inheritParams gboin_et_simulate
#' @inheritParams gboin_et_select_obd
#'
#' @param keep_trials
#'   Logical scalar. Keep the simulated trials and the selected OBDs in the
#'   result. Defaults to \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Report progress with messages. Defaults to \code{FALSE}.
#'
#' @return
#'   An object of class \code{boin_et_oc}, which is a list with components
#'   \item{design}{\code{"gBOIN-ET"}.}
#'   \item{target, target_resp}{The targets of the mean normalized scores.}
#'   \item{p_true, p_resp}{The true mean normalized toxicity and efficacy
#'     scores of the doses.}
#'   \item{n_doses, n_trials}{Numbers of doses and of simulated trials.}
#'   \item{lambda1, lambda2, eta1}{The boundaries.}
#'   \item{sel_percent}{Percentage of trials selecting each dose as the OBD.}
#'   \item{percent_no_obd}{Percentage of trials selecting no OBD.}
#'   \item{percent_early_stop}{Percentage of trials stopped early.}
#'   \item{n_pts_dose, n_tox_dose, n_resp_dose}{Average numbers of patients
#'     and average totals of the normalized toxicity and efficacy scores at
#'     each dose.}
#'   \item{total_n_pts, total_n_tox, total_n_resp}{Their averages per trial.}
#'   \item{n_grade_dose, n_resp_grade_dose}{Average numbers of patients of
#'     each toxicity category and of each efficacy category (rows) at each
#'     dose (columns).}
#'   \item{stop_reason_percent}{Percentage of trials ending for each reason,
#'     see \code{\link{boin_et_simulate}}.}
#'   \item{est_method, obd_method}{The OBD selection used.}
#'   \item{settings}{List of the design parameters.}
#'   \item{call}{The matched call.}
#'   \item{trials, obd}{The simulated trials and the selected OBDs, when
#'     \code{keep_trials} is \code{TRUE}.}
#'
#' @details
#'   The trials are simulated by \code{\link{gboin_et_simulate}} and the OBD
#'   is selected by \code{\link{gboin_et_select_obd}}; see there for the
#'   rules. The result prints like that of \code{\link{sim_boin_et}}, with
#'   scores in the place of DLTs and responses, and enters
#'   \code{\link{oc_results}} and \code{\link{oc_compare}} under the name
#'   gBOIN-ET.
#'
#'   The simulations of Takeda et al. (2022) estimate the efficacy scores by
#'   the linear fractional polynomial (\code{est_method = "fp_linear"}) and
#'   select the OBD among the doses with at least six patients
#'   (\code{min_obd_sample = 6}). The exclusion of a dose for its efficacy
#'   in Section 2.2 of the article differs from the default, the rule of
#'   'boinet' (see \code{\link{gboin_et_simulate}}), and the program of the
#'   authors has neither rule. Their program does not draw the categories
#'   from the probabilities of their Table S1: it draws one Bernoulli variate
#'   for each category and takes the highest category that comes up, so that
#'   the probabilities of the lower categories are smaller than tabulated.
#'   \code{inst/reproduce/reproduce_Takeda2022.R} gives the probabilities
#'   that the program draws and compares the results with Table 1 of the
#'   article.
#'
#' @references
#'   Takeda, K., Morita, S. and Taguri, M. (2022). gBOIN-ET: The generalized
#'   Bayesian optimal interval design for optimal dose-finding accounting for
#'   ordinal graded efficacy and toxicity in early clinical trials.
#'   Biometrical Journal, 64(7), 1178-1191. \doi{10.1002/bimj.202100263}
#'
#' @examples
#' # Toxicity grades 0 or 1, 2, 3 and 4, and stable or progressive disease,
#' # partial response and complete response, at four doses
#' p_grade <- rbind(
#'   c(0.83, 0.62, 0.51, 0.34),
#'   c(0.12, 0.18, 0.19, 0.16),
#'   c(0.04, 0.11, 0.14, 0.15),
#'   c(0.01, 0.09, 0.16, 0.35)
#' )
#' p_resp_grade <- rbind(
#'   c(0.90, 0.50, 0.30, 0.25),
#'   c(0.05, 0.35, 0.20, 0.20),
#'   c(0.05, 0.15, 0.50, 0.55)
#' )
#' oc <- sim_gboin_et(
#'   target = 0.47 / 1.5, target_resp = 0.95 / 1.5,
#'   p_grade = p_grade, p_resp_grade = p_resp_grade,
#'   score = c(0, 0.5, 1, 1.5), score_resp = c(0, 1, 1.5),
#'   n_cohort = 10, cohort_size = 3, n_trials = 1000, seed = 123
#' )
#' oc
#'
#' # The same trials with the OBD selection of the article
#' sim_gboin_et(
#'   target = 0.47 / 1.5, target_resp = 0.95 / 1.5,
#'   p_grade = p_grade, p_resp_grade = p_resp_grade,
#'   score = c(0, 0.5, 1, 1.5), score_resp = c(0, 1, 1.5),
#'   n_cohort = 10, cohort_size = 3, n_trials = 1000,
#'   est_method = "fp_linear", min_obd_sample = 6, seed = 123
#' )$sel_percent
#'
#' @seealso \code{\link{sim_gboin_et_multi}},
#'   \code{\link{gboin_et_simulate}}, \code{\link{gboin_et_select_obd}},
#'   \code{\link{sim_boin_et}}
#'
#' @export
sim_gboin_et <- function(target, target_resp, p_grade, p_resp_grade, score,
                         score_resp, n_cohort, cohort_size, n_trials = 10000,
                         start_dose = 1, n_earlystop = NULL, p_saf = NULL,
                         p_tox = NULL, p_saf_resp = NULL, lambda1 = NULL,
                         lambda2 = NULL, eta1 = NULL, cutoff_eli = 0.95,
                         cutoff_eli_resp = 0.99, p_eli = NULL,
                         p_eli_resp = NULL, resp_cor = 0.2,
                         boinet_compatible = FALSE, prior = rep(1 / 6, 6),
                         est_method = c("observed", "fp_linear",
                                        "fp_logistic"),
                         obd_method = c("max_resp", "utility_weighted",
                                        "utility_truncated_linear",
                                        "utility_scoring"),
                         min_obd_sample = 1, w1 = 0.33, w2 = 1.09,
                         t_low = NULL, t_upp = NULL, e_low = NULL,
                         e_upp = NULL, psi00 = 40, psi11 = 60,
                         keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  est_method <- match.arg(est_method)
  obd_method <- match.arg(obd_method)
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- gboin_et_simulate(
    target = target, target_resp = target_resp, p_grade = p_grade,
    p_resp_grade = p_resp_grade, score = score, score_resp = score_resp,
    n_cohort = n_cohort, cohort_size = cohort_size, n_trials = n_trials,
    start_dose = start_dose, n_earlystop = n_earlystop, p_saf = p_saf,
    p_tox = p_tox, p_saf_resp = p_saf_resp, lambda1 = lambda1,
    lambda2 = lambda2, eta1 = eta1, cutoff_eli = cutoff_eli,
    cutoff_eli_resp = cutoff_eli_resp, p_eli = p_eli,
    p_eli_resp = p_eli_resp, resp_cor = resp_cor,
    boinet_compatible = boinet_compatible, prior = prior, seed = seed
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
