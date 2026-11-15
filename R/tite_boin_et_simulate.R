#' Simulate TITE-BOIN-ET Trials
#'
#' @description
#'   Run the TITE-BOIN-ET dose-finding algorithm of Takeda, Morita and Taguri
#'   (2020) over many simulated trials, in which patients arrive over time and
#'   the DLT and the response are assessed over windows, and return the
#'   patient, DLT and response counts at every dose. This is the simulation
#'   engine behind \code{\link{sim_tite_boin_et}}, exposed for users who want
#'   the raw trial data.
#'
#' @inheritParams boin_et_simulate
#'
#' @param window
#'   Numeric scalar. Length of the DLT assessment window.
#'
#' @param accrual_rate
#'   Numeric scalar. Average number of patients arriving per unit of time, in
#'   the same unit as \code{window}. The \code{accrual} argument of the
#'   'boinet' package is the average time between arrivals, its inverse.
#'
#' @param resp_window
#'   Numeric scalar. Length of the response assessment window. Defaults to
#'   \code{window}.
#'
#' @param accrual
#'   Character scalar. Distribution of the times between arrivals:
#'   \code{"uniform"} (the default, as in 'boinet') on
#'   \code{(0, 2 / accrual_rate)}, \code{"exponential"} for a Poisson process,
#'   or \code{"fixed"} for arrivals exactly \code{1 / accrual_rate} apart.
#'
#' @param event_time
#'   Character scalar. Distribution of the times to DLT and to response:
#'   \code{"weibull"} (the default), with the share of the events in the
#'   second half of the window set by \code{late_fraction} and
#'   \code{resp_late_fraction}, or \code{"uniform"} on
#'   \code{(0, window / p)} for an event of probability \code{p}. With
#'   \code{"weibull"} the probabilities must be below 1.
#'
#' @param late_fraction
#'   Numeric scalar strictly between 0 and 1. Share of the DLTs that occur in
#'   the second half of the DLT window under the Weibull distribution.
#'   Defaults to 0.5.
#'
#' @param resp_late_fraction
#'   Numeric scalar strictly between 0 and 1. The same for the responses.
#'   Defaults to 0.5.
#'
#' @param resp_cor
#'   Numeric scalar strictly between -1 and 1. Correlation of the Gaussian
#'   copula that links the time to DLT and the time to response of a patient.
#'   Defaults to 0.2, the default of the 'boinet' package.
#'
#' @return
#'   An object of class \code{boin_et_trials}, as returned by
#'   \code{\link{boin_et_simulate}}, with \code{design} set to
#'   \code{"TITE-BOIN-ET"} and one more component
#'   \item{duration}{Numeric vector of the time of the last decision of each
#'     trial, which is the trial duration when all cohorts are enrolled.}
#'   The DLTs and responses are those observed at the last decision.
#'
#' @details
#'   The first patient of a cohort arrives at the time of the previous
#'   decision, and the next decision is due when the next patient arrives;
#'   after the last cohort it is due once every patient has been followed over
#'   both windows. While fewer than half of the patients at the current dose
#'   have finished both assessments, the decision waits for the next
#'   assessment to end. The DLT and response rates of a dose then count a
#'   patient still under assessment by the fraction of the window followed so
#'   far (Takeda et al., 2020), and the decision rule of BOIN-ET is applied to
#'   these rates. The admissibility of the doses, which counts the patients
#'   still under assessment as patients without the event, and the moves
#'   between doses follow \code{\link{boin_et_simulate}}, including the
#'   changes from 'boinet' unless \code{boinet_compatible} is \code{TRUE}.
#'   These are the rules of \code{boinet::tite.boinet()}.
#'
#'   The two event times of a patient come from their marginal distributions
#'   linked by a Gaussian copula, so that a DLT and a response occur within
#'   their windows with the probabilities \code{p_true} and \code{p_resp}.
#'   Each patient draws two uniform variates, so the trials agree with those
#'   of 'boinet' in distribution, not one by one.
#'
#' @references
#'   Takeda, K., Morita, S. and Taguri, M. (2020). TITE-BOIN-ET: Time-to-event
#'   Bayesian optimal interval design to accelerate dose-finding based on both
#'   efficacy and toxicity outcomes. Pharmaceutical Statistics, 19,
#'   335-349. \doi{10.1002/pst.1995}
#'
#'   Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
#'   Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
#'   and Toxicity Outcomes. The R Journal, 17(2), 15-34.
#'   \doi{10.32614/RJ-2025-012}
#'
#' @examples
#' trials <- tite_boin_et_simulate(
#'   target = 0.30, target_resp = 0.60,
#'   p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
#'   p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
#'   n_cohort = 12, cohort_size = 3, window = 30, resp_window = 45,
#'   accrual_rate = 1 / 10, n_trials = 200, seed = 123
#' )
#' trials
#' summary(trials$duration)
#'
#' @seealso \code{\link{sim_tite_boin_et}}, \code{\link{boin_et_select_obd}}
#' @export
tite_boin_et_simulate <- function(target, target_resp, p_true, p_resp,
                                  n_cohort, cohort_size, window, accrual_rate,
                                  resp_window = window,
                                  accrual = c("uniform", "exponential",
                                              "fixed"),
                                  event_time = c("weibull", "uniform"),
                                  late_fraction = 0.5,
                                  resp_late_fraction = 0.5, resp_cor = 0.2,
                                  n_trials = 10000, start_dose = 1,
                                  n_earlystop = NULL, p_saf = NULL,
                                  p_tox = NULL, p_saf_resp = NULL,
                                  cutoff_eli = 0.95, cutoff_eli_resp = 0.99,
                                  p_eli = NULL, p_eli_resp = NULL,
                                  boinet_compatible = FALSE,
                                  prior = rep(1 / 6, 6), seed = 123) {

  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  check_flag(boinet_compatible, "boinet_compatible")

  check_p_true(p_true)
  n_doses <- length(p_true)
  if (!is.numeric(p_resp) || length(p_resp) != n_doses ||
      any(!is.finite(p_resp)) || any(p_resp < 0) || any(p_resp > 1)) {
    stop("'p_resp' must hold one response probability between 0 and 1 ",
         "per dose in 'p_true'", call. = FALSE)
  }
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(start_dose, "start_dose", 1L)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of doses in 'p_true'",
         call. = FALSE)
  }
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  check_scalar_prob(resp_cor, "resp_cor", lower = -1, upper = 1)
  boin_et_check_timing(window, resp_window, accrual_rate, late_fraction,
                       resp_late_fraction)
  par_t <- boin_et_event_params(p_true, window, late_fraction, event_time,
                                "p_true")
  par_e <- boin_et_event_params(p_resp, resp_window, resp_late_fraction,
                                event_time, "p_resp")

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)
  if (is.null(n_earlystop)) n_earlystop <- max_total_pts
  check_count(n_earlystop, "n_earlystop", 1L)

  bd <- boin_et_boundary(target, target_resp, p_saf = p_saf, p_tox = p_tox,
                         p_saf_resp = p_saf_resp, prior = prior)
  eli <- boin_et_check_eli(p_eli, p_eli_resp, bd)
  adm <- boin_et_admissible_table(eli$p_eli, eli$p_eli_resp, cutoff_eli,
                                  cutoff_eli_resp, max_total_pts)

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- tite_boin_et_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    p_resp = as.numeric(p_resp),
    resp_cor = as.numeric(resp_cor),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    n_earlystop = as.integer(n_earlystop),
    tox_ok = as.integer(adm$tox_ok),
    resp_ok = as.integer(adm$resp_ok),
    max_n = as.integer(max_total_pts),
    boinet_compatible = as.integer(boinet_compatible),
    lambda1 = bd$lambda1,
    lambda2 = bd$lambda2,
    eta1 = bd$eta1,
    window = as.numeric(window),
    resp_window = as.numeric(resp_window),
    accrual_rate = as.numeric(accrual_rate),
    arrival = match(accrual, c("uniform", "exponential", "fixed")) - 1L,
    event_time = match(event_time, c("weibull", "uniform")) - 1L,
    shape_t = par_t$shape,
    scale_t = par_t$scale,
    shape_e = par_e$shape,
    scale_e = par_e$scale
  )

  dose_names <- paste0("DL", seq_len(n_doses))
  colnames(res$n_pts) <- dose_names
  colnames(res$n_tox) <- dose_names
  colnames(res$n_resp) <- dose_names

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$n_tox,
      n_resp = res$n_resp,
      cohorts_used = res$cohorts_used,
      stop_reason = boin_et_stop_reasons()[res$stop_code + 1L],
      early_stop = res$stop_code <= 1L,
      duration = res$duration,
      design = "TITE-BOIN-ET",
      boundary = bd,
      settings = list(
        target = bd$target,
        target_resp = bd$target_resp,
        p_true = p_true,
        p_resp = p_resp,
        n_cohort = n_cohort,
        cohort_size = cohort_size_vec,
        window = window,
        resp_window = resp_window,
        accrual_rate = accrual_rate,
        accrual = accrual,
        event_time = event_time,
        late_fraction = late_fraction,
        resp_late_fraction = resp_late_fraction,
        n_trials = n_trials,
        start_dose = start_dose,
        n_earlystop = n_earlystop,
        p_saf = bd$p_saf,
        p_tox = bd$p_tox,
        p_saf_resp = bd$p_saf_resp,
        cutoff_eli = cutoff_eli,
        cutoff_eli_resp = cutoff_eli_resp,
        p_eli = eli$p_eli,
        p_eli_resp = eli$p_eli_resp,
        resp_cor = resp_cor,
        boinet_compatible = boinet_compatible,
        prior = prior,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "boin_et_trials"
  )
}
