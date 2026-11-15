#' Operating Characteristics of the BOIN-ET Design
#'
#' @description
#'   Simulate trials of the BOIN-ET design of Takeda, Taguri and Morita (2018)
#'   under one scenario of DLT and response probabilities, select the optimal
#'   biological dose (OBD) of each trial, and summarize how often each dose is
#'   selected and how many patients, DLTs and responses occur at each dose.
#'
#' @inheritParams boin_et_simulate
#' @inheritParams boin_et_select_obd
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
#'   \item{target, target_resp, p_true, p_resp}{The targets and the scenario.}
#'   \item{n_doses, n_trials}{Numbers of doses and of simulated trials.}
#'   \item{lambda1, lambda2, eta1}{The boundaries, see
#'     \code{\link{boin_et_boundary}}.}
#'   \item{sel_percent}{Percentage of trials selecting each dose as the OBD.}
#'   \item{percent_no_obd}{Percentage of trials selecting no OBD.}
#'   \item{percent_early_stop}{Percentage of trials stopped early.}
#'   \item{n_pts_dose, n_tox_dose, n_resp_dose}{Average numbers of patients,
#'     DLTs and responses at each dose.}
#'   \item{total_n_pts, total_n_tox, total_n_resp}{Their averages per trial.}
#'   \item{stop_reasons}{Percentage of trials ending for each reason, see
#'     \code{\link{boin_et_simulate}}.}
#'   \item{est_method, obd_method}{The OBD selection used.}
#'   \item{settings}{List of the design parameters.}
#'   \item{call}{The matched call.}
#'   \item{trials, obd}{The simulated trials and the selected OBDs, when
#'     \code{keep_trials} is \code{TRUE}.}
#'
#' @details
#'   The trials are simulated by \code{\link{boin_et_simulate}} and the OBD is
#'   selected by \code{\link{boin_et_select_obd}}; see there for the rules. The
#'   defaults are those of \code{boinet::boinet()}: observed response rates,
#'   and the admissible dose with the highest response rate up to the
#'   estimated MTD.
#'
#' @references
#'   Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
#'   interval design for dose finding based on both efficacy and toxicity
#'   outcomes. Pharmaceutical Statistics, 17(4), 383-395.
#'   \doi{10.1002/pst.1864}
#'
#'   Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
#'   Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
#'   and Toxicity Outcomes. The R Journal, 17(2), 15-34.
#'   \doi{10.32614/RJ-2025-012}
#'
#' @examples
#' # Scenario 1 of Takeda et al. (2018)
#' oc <- sim_boin_et(
#'   target = 0.30, target_resp = 0.60,
#'   p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
#'   p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
#'   n_cohort = 12, cohort_size = 3, n_trials = 1000, seed = 123
#' )
#' oc
#'
#' # The same trials under the scoring utility
#' sim_boin_et(
#'   target = 0.30, target_resp = 0.60,
#'   p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
#'   p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
#'   n_cohort = 12, cohort_size = 3, n_trials = 1000,
#'   obd_method = "utility_scoring", seed = 123
#' )$sel_percent
#'
#' @seealso \code{\link{sim_boin_et_multi}}, \code{\link{boin_et_simulate}},
#'   \code{\link{boin_et_select_obd}}
#'
#' @export
sim_boin_et <- function(target, target_resp, p_true, p_resp, n_cohort,
                        cohort_size, n_trials = 10000, start_dose = 1,
                        n_earlystop = NULL, p_saf = NULL, p_tox = NULL,
                        p_saf_resp = NULL, cutoff_eli = 0.95,
                        cutoff_eli_resp = 0.99, p_eli = NULL,
                        p_eli_resp = NULL, resp_cor = 0.2,
                        resp_cor_type = c("copula", "binary"),
                        prior = rep(1 / 6, 6),
                        est_method = c("observed", "multi_iso", "fp_logistic"),
                        obd_method = c("max_resp", "utility_weighted",
                                       "utility_truncated_linear",
                                       "utility_scoring"),
                        w1 = 0.33, w2 = 1.09, t_low = NULL, t_upp = NULL,
                        e_low = NULL, e_upp = NULL, psi00 = 40, psi11 = 60,
                        keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  resp_cor_type <- match.arg(resp_cor_type)
  est_method <- match.arg(est_method)
  obd_method <- match.arg(obd_method)
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- boin_et_simulate(
    target = target, target_resp = target_resp, p_true = p_true,
    p_resp = p_resp, n_cohort = n_cohort, cohort_size = cohort_size,
    n_trials = n_trials, start_dose = start_dose, n_earlystop = n_earlystop,
    p_saf = p_saf, p_tox = p_tox, p_saf_resp = p_saf_resp,
    cutoff_eli = cutoff_eli, cutoff_eli_resp = cutoff_eli_resp,
    p_eli = p_eli, p_eli_resp = p_eli_resp, resp_cor = resp_cor,
    resp_cor_type = resp_cor_type, prior = prior, seed = seed
  )

  if (verbose) message("Selecting the OBD ...")

  selection <- boin_et_select_obd(
    n_pts = trials$n_pts, n_tox = trials$n_tox, n_resp = trials$n_resp,
    target = target, target_resp = target_resp, p_saf = p_saf,
    p_tox = p_tox, p_saf_resp = p_saf_resp, cutoff_eli = cutoff_eli,
    cutoff_eli_resp = cutoff_eli_resp, p_eli = p_eli, p_eli_resp = p_eli_resp,
    est_method = est_method,
    obd_method = obd_method, w1 = w1, w2 = w2, t_low = t_low, t_upp = t_upp,
    e_low = e_low, e_upp = e_upp, psi00 = psi00, psi11 = psi11,
    early_stop = trials$early_stop
  )

  n_doses <- length(p_true)
  dose_names <- paste0("DL", seq_len(n_doses))

  selected <- factor(selection$obd, levels = seq_len(n_doses))
  sel_percent <- as.numeric(table(selected)) / n_trials * 100
  names(sel_percent) <- dose_names

  n_pts_dose <- colMeans(trials$n_pts)
  n_tox_dose <- colMeans(trials$n_tox)
  n_resp_dose <- colMeans(trials$n_resp)
  names(n_pts_dose) <- dose_names
  names(n_tox_dose) <- dose_names
  names(n_resp_dose) <- dose_names

  reasons <- c("no_admissible_dose", "no_lower_admissible_dose",
               "n_earlystop", "max_cohorts")
  stop_reasons <- as.numeric(table(factor(trials$stop_reason,
                                          levels = reasons))) / n_trials * 100
  names(stop_reasons) <- reasons

  settings <- trials$settings
  settings$w1 <- w1
  settings$w2 <- w2
  settings$t_low <- t_low
  settings$t_upp <- t_upp
  settings$e_low <- e_low
  settings$e_upp <- e_upp
  settings$psi00 <- psi00
  settings$psi11 <- psi11

  out <- list(
    target = target,
    target_resp = target_resp,
    p_true = p_true,
    p_resp = p_resp,
    n_doses = n_doses,
    n_trials = n_trials,
    lambda1 = trials$boundary$lambda1,
    lambda2 = trials$boundary$lambda2,
    eta1 = trials$boundary$eta1,
    sel_percent = sel_percent,
    percent_no_obd = mean(is.na(selection$obd)) * 100,
    percent_early_stop = mean(trials$early_stop) * 100,
    n_pts_dose = n_pts_dose,
    n_tox_dose = n_tox_dose,
    n_resp_dose = n_resp_dose,
    total_n_pts = sum(n_pts_dose),
    total_n_tox = sum(n_tox_dose),
    total_n_resp = sum(n_resp_dose),
    stop_reasons = stop_reasons,
    est_method = est_method,
    obd_method = obd_method,
    settings = settings,
    call = call_expr
  )
  if (keep_trials) {
    out$trials <- trials
    out$obd <- selection$obd
  }

  structure(out, class = "boin_et_oc")
}
