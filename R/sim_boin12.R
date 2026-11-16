#' Operating Characteristics of the BOIN12 Design
#'
#' @description
#'   Simulate trials of the BOIN12 design of Lin et al. (2020) under one
#'   scenario of DLT and response probabilities, select the optimal
#'   biological dose (OBD) of each trial, and summarize how often each dose is
#'   selected and how many patients, DLTs and responses occur at each dose.
#'
#' @inheritParams boin12_simulate
#' @inheritParams boin12_select_obd
#'
#' @param keep_trials
#'   Logical scalar. Keep the simulated trials and the selected OBDs in the
#'   result. Defaults to \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Report progress with messages. Defaults to \code{FALSE}.
#'
#' @return
#'   An object of class \code{boin12_oc}, which is a list with components
#'   \item{design}{\code{"BOIN12"}.}
#'   \item{target, target_resp, p_true, p_resp}{The limits and the scenario.}
#'   \item{n_doses, n_trials}{Numbers of doses and of simulated trials.}
#'   \item{lambda_e, lambda_d, u_benchmark}{The BOIN boundaries and the
#'     utility benchmark.}
#'   \item{sel_percent}{Percentage of trials selecting each dose as the OBD.}
#'   \item{percent_no_obd}{Percentage of trials selecting no OBD.}
#'   \item{percent_early_stop}{Percentage of trials stopped because no dose
#'     was admissible.}
#'   \item{n_pts_dose, n_tox_dose, n_resp_dose}{Average numbers of patients,
#'     DLTs and responses at each dose.}
#'   \item{total_n_pts, total_n_tox, total_n_resp}{Their averages per trial.}
#'   \item{stop_reason_percent}{Percentage of trials ending for each reason,
#'     see \code{\link{boin12_simulate}}.}
#'   \item{est_method}{The estimate of the utility used to select the OBD.}
#'   \item{settings}{List of the design parameters.}
#'   \item{call}{The matched call.}
#'   \item{trials, obd}{The simulated trials and the selected OBDs, when
#'     \code{keep_trials} is \code{TRUE}.}
#'
#' @details
#'   The trials are simulated by \code{\link{boin12_simulate}} and the OBD is
#'   selected by \code{\link{boin12_select_obd}}; see there for the rules.
#'
#' @references
#'   Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
#'   optimal interval phase I/II trial design for utility-based dose finding
#'   in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
#'   1393-1402. \doi{10.1200/PO.20.00257}
#'
#' @examples
#' # Scenario 1 of Lin et al. (2020), 36 patients in cohorts of 3
#' oc <- sim_boin12(
#'   target = 0.35, target_resp = 0.25,
#'   p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'   p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
#'   n_cohort = 12, cohort_size = 3, n_trials = 1000, seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_boin12_multi}}, \code{\link{boin12_simulate}},
#'   \code{\link{boin12_select_obd}}
#'
#' @export
sim_boin12 <- function(target, target_resp, p_true, p_resp, n_cohort,
                       cohort_size, n_trials = 10000, start_dose = 1,
                       n_earlystop = NULL,
                       n_earlystop_rule = c("with_stay", "simple"),
                       p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                       cutoff_eli_resp = 0.90, psi00 = 40, psi11 = 60,
                       benchmark_weight = 0.5, utility_prior = c(1, 1),
                       n_star = 6, n_explore = 9, resp_cor = 0,
                       resp_cor_type = c("copula", "binary"),
                       est_method = c("posterior_mean", "observed",
                                      "tail_probability"),
                       keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  n_earlystop_rule <- match.arg(n_earlystop_rule)
  resp_cor_type <- match.arg(resp_cor_type)
  est_method <- match.arg(est_method)
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- boin12_simulate(
    target = target, target_resp = target_resp, p_true = p_true,
    p_resp = p_resp, n_cohort = n_cohort, cohort_size = cohort_size,
    n_trials = n_trials, start_dose = start_dose, n_earlystop = n_earlystop,
    n_earlystop_rule = n_earlystop_rule, p_saf = p_saf, p_tox = p_tox,
    cutoff_eli = cutoff_eli, cutoff_eli_resp = cutoff_eli_resp,
    psi00 = psi00, psi11 = psi11, benchmark_weight = benchmark_weight,
    utility_prior = utility_prior, n_star = n_star, n_explore = n_explore,
    resp_cor = resp_cor, resp_cor_type = resp_cor_type, seed = seed
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
