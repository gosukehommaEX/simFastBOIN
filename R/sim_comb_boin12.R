#' Operating Characteristics of the Comb-BOIN12 Design
#'
#' @description
#'   Simulate trials of the Comb-BOIN12 design of Lu et al. (2025) for two
#'   drugs given together with \code{\link{comb_boin12_simulate}}, select the
#'   optimal biological dose combination (OBDC) of every trial with
#'   \code{\link{comb_boin12_select_obd}}, and summarize the selection
#'   percentages and the numbers of patients, DLTs and responses at every
#'   combination.
#'
#' @inheritParams comb_boin12_simulate
#' @inheritParams comb_boin12_select_obd
#'
#' @param keep_trials
#'   Logical scalar. Keep the simulated trials (\code{trials}) and the
#'   selected combinations (\code{obd}) in the result. Defaults to
#'   \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Print progress messages. Defaults to \code{FALSE}.
#'
#' @return
#'   An object of class \code{comb_boin12_oc}, which is a list with
#'   components
#'   \item{design}{\code{"Comb-BOIN12"}.}
#'   \item{target, target_resp, p_true, p_resp}{The limits and the true
#'     probabilities.}
#'   \item{n_rows, n_cols, n_trials}{Numbers of levels of drug A and drug B,
#'     and of trials.}
#'   \item{lambda_e, lambda_d, u_benchmark}{The BOIN boundaries and the
#'     utility benchmark (scale 0 to 100).}
#'   \item{sel_percent}{Matrix of the percentage of trials selecting each
#'     combination as the OBDC.}
#'   \item{percent_no_obd, percent_early_stop}{Percentages of trials
#'     selecting no OBDC and stopped early.}
#'   \item{n_pts_dose, n_tox_dose, n_resp_dose}{Matrices of the average
#'     numbers of patients, DLTs and responses at each combination.}
#'   \item{total_n_pts, total_n_tox, total_n_resp}{Their totals.}
#'   \item{true_obd}{Two-column matrix of the true OBDCs: the combinations
#'     with the highest true utility among those with a DLT probability at
#'     most \code{target}.}
#'   \item{pct_true_obd_selected, pct_patients_true_obd}{Percentage of
#'     trials selecting a true OBDC and percentage of the patients treated
#'     there.}
#'   \item{overdose}{List with the \code{cutoff} (\code{target}), the
#'     \code{combinations} above it, and the percentages of trials selecting
#'     one of them (\code{pct_selected}) and of patients treated there
#'     (\code{pct_patients}).}
#'   \item{stop_reason_percent}{Percentage of trials by stopping reason.}
#'   \item{est_method, obd_bound, settings, call}{The settings used.}
#'   \item{trials, obd}{With \code{keep_trials = TRUE}, the
#'     \code{comb_boin12_trials} object and a two-column matrix of the
#'     selected combinations.}
#'
#' @details
#'   The percentages of patients are those of the average numbers of
#'   patients. Lu et al. (2025) report the average over trials of the
#'   number of patients as a percentage of the maximum sample size; with
#'   \code{keep_trials = TRUE} it can be computed from \code{trials}, as in
#'   \code{inst/reproduce/reproduce_Lu2025.R}.
#'
#'   With the logistic estimate, each trial takes a few milliseconds to
#'   select; \code{est_method = "observed"} is much faster.
#'
#' @references
#'   Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
#'   utility-based Bayesian optimal interval design for dose optimization in
#'   cancer drug-combination trials. Statistics in Biopharmaceutical
#'   Research, 17(2), 266-276. \doi{10.1080/19466315.2024.2370403}
#'
#' @examples
#' # Scenario 16 of Lu et al. (2025), with their design settings
#' p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
#' p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
#' oc <- sim_comb_boin12(
#'   target = 0.33, target_resp = 0.2, p_true = p_true, p_resp = p_resp,
#'   n_cohort = 14, cohort_size = 3, n_trials = 200,
#'   psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11, benchmark_weight = 0.6,
#'   extrasafe = TRUE, offset = 0.10, dose_a = c(100, 200),
#'   dose_b = c(50, 100, 200), seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_comb_boin12_multi}},
#'   \code{\link{comb_boin12_simulate}}, \code{\link{comb_boin12_select_obd}},
#'   \code{\link{print.comb_boin12_oc}}
#'
#' @export
sim_comb_boin12 <- function(target, target_resp, p_true, p_resp, n_cohort,
                            cohort_size, n_trials = 10000,
                            start_dose = c(1, 1), n_earlystop = NULL,
                            n_earlystop_rule = c("with_stay", "simple"),
                            p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                            cutoff_eli_resp = 0.90, extrasafe = FALSE,
                            offset = 0.05, psi00 = 40, psi11 = 60,
                            benchmark_weight = 0.5, utility_prior = c(1, 1),
                            n_star = 6, tie_break = c("higher", "random"),
                            no_candidate = c("nearest", "stop"),
                            resp_cor = 0,
                            resp_cor_type = c("copula", "binary"),
                            est_method = c("logistic", "observed"),
                            obd_bound = c("rectangle", "contour"),
                            dose_a = NULL, dose_b = NULL,
                            logistic_scale = c(2.5, 2.5), n_draws = 4000,
                            keep_trials = FALSE, verbose = FALSE,
                            seed = 123) {

  call_expr <- match.call()
  n_earlystop_rule <- match.arg(n_earlystop_rule)
  tie_break <- match.arg(tie_break)
  no_candidate <- match.arg(no_candidate)
  resp_cor_type <- match.arg(resp_cor_type)
  est_method <- match.arg(est_method)
  obd_bound <- match.arg(obd_bound)
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- comb_boin12_simulate(
    target = target, target_resp = target_resp, p_true = p_true,
    p_resp = p_resp, n_cohort = n_cohort, cohort_size = cohort_size,
    n_trials = n_trials, start_dose = start_dose, n_earlystop = n_earlystop,
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
