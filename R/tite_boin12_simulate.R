#' Simulate TITE-BOIN12 Trials
#'
#' @description
#'   Simulate trials of the TITE-BOIN12 design of Zhou et al. (2022), the
#'   time-to-event version of BOIN12 in which patients arrive over time, the
#'   DLT and the response are assessed over windows and the decisions use the
#'   outcomes still pending through the approximated likelihood. Return the
#'   patients, DLTs and responses at every dose. This is the simulation engine
#'   behind \code{\link{sim_tite_boin12}}.
#'
#' @inheritParams boin12_simulate
#'
#' @param window
#'   Numeric scalar. Length of the DLT assessment window.
#'
#' @param accrual_rate
#'   Numeric scalar. Average number of patients arriving per unit of time, in
#'   the same unit as \code{window}.
#'
#' @param resp_window
#'   Numeric scalar. Length of the response assessment window. Defaults to
#'   \code{window}.
#'
#' @param accrual
#'   Character scalar. Distribution of the times between arrivals:
#'   \code{"uniform"} (the default) on \code{(0, 2 / accrual_rate)},
#'   \code{"exponential"} for a Poisson process, or \code{"fixed"} for
#'   arrivals exactly \code{1 / accrual_rate} apart.
#'
#' @param event_time
#'   Character scalar. Distribution of the times to DLT and to response:
#'   \code{"weibull"} (the default), with the share of the events in the
#'   second half of the window set by \code{late_fraction} and
#'   \code{resp_late_fraction}, or \code{"uniform"} over the window. With
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
#' @param min_completed_share
#'   Numeric scalar from 0 to 1. Share of the patients at the current dose who
#'   must have finished both assessments before the next decision; the
#'   decision, and the accrual, wait until then. Defaults to 0.5, the rule of
#'   Zhou et al. (2022), who suspend the accrual when more than half of the
#'   patients at the current dose have a pending outcome. With 1 every
#'   decision waits for all the outcomes, as in BOIN12.
#'
#' @param pending_estimate
#'   Character scalar. Value that replaces a pending outcome Y in the
#'   quasi-binomial likelihood: \code{"conditional"} (the default), the
#'   probability of the event given that it has not occurred during the time
#'   followed, \code{pi_hat * (1 - w) / (1 - pi_hat * w)}, with \code{w} the
#'   fraction of the window followed, as in Equation (11) of Zhou et al.
#'   (2022); or \code{"marginal"}, the estimate \code{pi_hat} itself, which
#'   reproduces the desirabilities of Table S.3 of their Supplement. See
#'   Details.
#'
#' @param resp_cor
#'   Numeric scalar between -1 and 1. Correlation between the DLT and the
#'   response of a patient, see \code{resp_cor_type}; with
#'   \code{"copula"} it also links the time to DLT and the time to response.
#'   Defaults to 0.
#'
#' @return
#'   An object of class \code{boin12_trials}, as returned by
#'   \code{\link{boin12_simulate}}, with \code{design} set to
#'   \code{"TITE-BOIN12"} and one more component
#'   \item{duration}{Numeric vector of the duration of each trial: the time at
#'     which every patient has finished both assessments, or the time of the
#'     decision that stopped the trial because no dose was admissible or, with
#'     \code{no_candidate = "stop"}, because no dose allowed for the next
#'     cohort was admissible.}
#'   The DLTs and responses are the final outcomes of all the patients.
#'
#' @details
#'   The first patient of a cohort arrives at the time of the previous
#'   decision, and the next decision is due when the next patient arrives.
#'   While fewer than a share \code{min_completed_share} of the patients at
#'   the current dose have finished both assessments, the decision waits for
#'   the next assessment to end. After the last cohort, and when the trial
#'   stops because \code{n_earlystop} is reached, the trial waits until every
#'   patient has finished both assessments, and the OBD is selected from the
#'   complete data, as in the example of Zhou et al. (2022).
#'
#'   At a decision, an outcome is observed once its event has occurred or its
#'   window has ended, and pending otherwise. For each dose, with \code{w} the
#'   fraction of the window followed by a pending patient, the DLT
#'   probability is estimated by \code{pi_hat} = DLTs observed / (DLTs
#'   observed + patients who completed the window without DLT + sum of
#'   \code{w}) (Equation 14 of the Supplement), and the response probability
#'   in the same way. The quasi-events of BOIN12 count each patient with a
#'   pending outcome by the expected utility, the pending outcome being
#'   replaced as set by \code{pending_estimate} and the two outcomes taken as
#'   independent for this purpose only. The desirability is that of BOIN12
#'   with these quasi-events and the number of patients treated, and the
#'   rules of \code{\link{boin12_simulate}} are applied with \code{pi_hat} of
#'   the DLTs as the DLT rate of the current dose. With every outcome
#'   observed, these are exactly the quantities of BOIN12, so that
#'   \code{min_completed_share = 1} with \code{accrual = "fixed"} gives the
#'   trials of \code{\link{boin12_simulate}} with the same seed.
#'
#'   Two points are not settled by Zhou et al. (2022) and are decided here as
#'   follows. Equation (11) replaces a pending outcome by its conditional
#'   expectation, but the desirabilities of their Table S.3 follow from
#'   \code{pi_hat} itself; both are offered by \code{pending_estimate}. The
#'   admissibility of a dose, which they update at every decision, uses the
#'   posterior of the approximated likelihood, Beta(1 + DLTs observed, 1 +
#'   patients who completed the window without DLT + sum of \code{w}), and
#'   the same for the responses, so that a response that is merely late does
#'   not exclude a dose. TITE-BOIN-ET (\code{\link{tite_boin_et_simulate}})
#'   instead counts the pending patients as patients without the event, as
#'   the 'boinet' package does.
#'
#'   Each patient draws two uniform variates as in
#'   \code{\link{boin12_simulate}}; the time to DLT and the time to response
#'   are quantiles of their distributions at these variates, so that a DLT
#'   and a response occur within their windows with the probabilities
#'   \code{p_true} and \code{p_resp}. The times between arrivals draw further
#'   variates unless \code{accrual = "fixed"}, and with
#'   \code{tie_break = "random"} a decision among tied doses draws one more
#'   variate after those of the cohort, as in \code{\link{boin12_simulate}}.
#'
#' @references
#'   Zhou, Y., Lin, R., Lee, J. J., Li, D., Wang, L., Li, R. and Yuan, Y.
#'   (2022). TITE-BOIN12: A Bayesian phase I/II trial design to find the
#'   optimal biological dose with late-onset toxicity and efficacy. Statistics
#'   in Medicine, 41(11), 1918-1931. \doi{10.1002/sim.9337}
#'
#'   Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
#'   optimal interval phase I/II trial design for utility-based dose finding
#'   in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
#'   1393-1402. \doi{10.1200/PO.20.00257}
#'
#' @examples
#' # Scenario 1 of Zhou et al. (2022), windows of 45 and 60 days and three
#' # patients a month
#' trials <- tite_boin12_simulate(
#'   target = 0.35, target_resp = 0.25,
#'   p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'   p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
#'   n_cohort = 12, cohort_size = 3, window = 45, resp_window = 60,
#'   accrual_rate = 3 / 30, n_trials = 200, seed = 123
#' )
#' trials
#' summary(trials$duration)
#'
#' @seealso \code{\link{sim_tite_boin12}}, \code{\link{boin12_simulate}},
#'   \code{\link{boin12_select_obd}}
#'
#' @export
tite_boin12_simulate <- function(target, target_resp, p_true, p_resp,
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
                                 n_trials = 10000, start_dose = 1,
                                 n_earlystop = NULL,
                                 n_earlystop_rule = c("with_stay", "simple"),
                                 p_saf = NULL, p_tox = NULL,
                                 cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
                                 psi00 = 40, psi11 = 60,
                                 benchmark_weight = 0.5,
                                 utility_prior = c(1, 1), n_star = 6,
                                 n_explore = 9,
                                 tie_break = c("higher", "random"),
                                 no_candidate = c("nearest", "stop"),
                                 resp_cor = 0,
                                 resp_cor_type = c("copula", "binary"),
                                 seed = 123) {

  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  pending_estimate <- match.arg(pending_estimate)
  n_earlystop_rule <- match.arg(n_earlystop_rule)
  tie_break <- match.arg(tie_break)
  no_candidate <- match.arg(no_candidate)
  resp_cor_type <- match.arg(resp_cor_type)

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
  check_count(n_star, "n_star", 1L)
  check_count(n_explore, "n_explore", 1L)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  check_scalar_prob(resp_cor, "resp_cor", lower = -1, upper = 1)
  if (!is.numeric(min_completed_share) || length(min_completed_share) != 1L ||
      !is.finite(min_completed_share) || min_completed_share < 0 ||
      min_completed_share > 1) {
    stop("'min_completed_share' must be a single number from 0 to 1",
         call. = FALSE)
  }
  boin_et_check_timing(window, resp_window, accrual_rate, late_fraction,
                       resp_late_fraction)
  par_t <- boin_et_event_params(p_true, window, late_fraction, event_time,
                                "p_true")
  par_e <- boin_et_event_params(p_resp, resp_window, resp_late_fraction,
                                event_time, "p_resp")
  des <- boin12_check_design(target, target_resp, p_saf = p_saf,
                             p_tox = p_tox, psi00 = psi00, psi11 = psi11,
                             benchmark_weight = benchmark_weight,
                             utility_prior = utility_prior)

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)
  if (!is.null(n_earlystop)) check_count(n_earlystop, "n_earlystop", 1L)

  if (resp_cor_type == "binary") {
    link <- boin_et_resp_link(p_true, p_resp, resp_cor)
  } else {
    link <- list(q0 = rep(0, n_doses), q1 = rep(0, n_doses))
  }

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- tite_boin12_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    p_resp = as.numeric(p_resp),
    q0 = as.numeric(link$q0),
    q1 = as.numeric(link$q1),
    cor_type = if (resp_cor_type == "binary") 1L else 0L,
    resp_cor = as.numeric(resp_cor),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    n_earlystop = if (is.null(n_earlystop)) 0L else as.integer(n_earlystop),
    early_stop_simple = identical(n_earlystop_rule, "simple"),
    n_star = as.integer(n_star),
    n_explore = as.integer(n_explore),
    tie_break = if (tie_break == "random") 1L else 0L,
    no_candidate = if (no_candidate == "stop") 1L else 0L,
    lambda_e = des$lambda_e,
    lambda_d = des$lambda_d,
    psi00 = as.numeric(psi00),
    psi11 = as.numeric(psi11),
    alpha = des$utility_prior[1L],
    beta = des$utility_prior[2L],
    u_benchmark = des$u_benchmark / 100,
    target = as.numeric(target),
    target_resp = as.numeric(target_resp),
    cutoff_eli = as.numeric(cutoff_eli),
    cutoff_eli_resp = as.numeric(cutoff_eli_resp),
    pending = if (pending_estimate == "marginal") 1L else 0L,
    window = as.numeric(window),
    resp_window = as.numeric(resp_window),
    accrual_rate = as.numeric(accrual_rate),
    arrival = match(accrual, c("uniform", "exponential", "fixed")) - 1L,
    event_time = match(event_time, c("weibull", "uniform")) - 1L,
    shape_t = par_t$shape,
    scale_t = par_t$scale,
    shape_e = par_e$shape,
    scale_e = par_e$scale,
    min_completed_share = as.numeric(min_completed_share)
  )

  dose_names <- paste0("DL", seq_len(n_doses))
  for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
    colnames(res[[nm]]) <- dose_names
  }

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$n_tox,
      n_resp = res$n_resp,
      n_tox_resp = res$n_tox_resp,
      cohorts_used = res$cohorts_used,
      stop_reason = boin12_stop_reasons()[res$stop_code + 1L],
      early_stop = res$stop_code %in% c(0L, 3L),
      duration = res$duration,
      design = "TITE-BOIN12",
      boundary = list(lambda_e = des$lambda_e, lambda_d = des$lambda_d,
                      u_low = des$u_low, u_benchmark = des$u_benchmark),
      settings = list(
        target = target,
        target_resp = target_resp,
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
        min_completed_share = min_completed_share,
        pending_estimate = pending_estimate,
        n_trials = n_trials,
        start_dose = start_dose,
        n_earlystop = n_earlystop,
        n_earlystop_rule = n_earlystop_rule,
        p_saf = des$p_saf,
        p_tox = des$p_tox,
        cutoff_eli = cutoff_eli,
        cutoff_eli_resp = cutoff_eli_resp,
        psi00 = psi00,
        psi11 = psi11,
        benchmark_weight = benchmark_weight,
        utility_prior = des$utility_prior,
        n_star = n_star,
        n_explore = n_explore,
        tie_break = tie_break,
        no_candidate = no_candidate,
        resp_cor = resp_cor,
        resp_cor_type = resp_cor_type,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "boin12_trials"
  )
}
