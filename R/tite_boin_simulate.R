#' Simulate TITE-BOIN Trials
#'
#' @description
#'   Run the time-to-event BOIN (TITE-BOIN) design over many simulated trials,
#'   in which patients arrive over time and new patients are treated while the
#'   DLT assessment of earlier patients is still pending. Return the patient and
#'   DLT counts at every dose together with the duration of each trial and the
#'   time spent with accrual suspended.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param p_true
#'   Numeric vector. True DLT probability at each dose level, in increasing dose
#'   order. The probability of a DLT within the assessment window.
#'
#' @param n_cohort
#'   Integer scalar. Number of cohorts in a trial.
#'
#' @param cohort_size
#'   Integer scalar or vector. Number of patients per cohort. A scalar is used for
#'   every cohort. A vector shorter than \code{n_cohort} is padded with its last
#'   element and a longer one is truncated.
#'
#' @param window
#'   Numeric scalar. Length of the DLT assessment window.
#'
#' @param accrual_rate
#'   Numeric scalar. Average number of patients arriving per unit of time, in the
#'   same unit as \code{window}.
#'
#' @param method
#'   Character scalar, \code{"imputation"} (the default) or \code{"ess"}. How the
#'   pending patients enter the estimate of the DLT probability. See
#'   \code{\link{tite_boin_decision_table}}.
#'
#' @param accrual
#'   Character scalar. Distribution of the times between arrivals:
#'   \code{"exponential"} (the default) for a Poisson process,
#'   \code{"uniform"} for a uniform distribution on
#'   \code{(0, 2 / accrual_rate)}, as in the code of Lin and Yuan (2020), or
#'   \code{"fixed"} for arrivals exactly \code{1 / accrual_rate} apart.
#'
#' @param dlt_time
#'   Character scalar. Distribution of the time to DLT within the window:
#'   \code{"weibull"} (the default), with the share of DLTs in the second half
#'   of the window set by \code{late_fraction}, or \code{"uniform"}.
#'
#' @param late_fraction
#'   Numeric scalar strictly between 0 and 1. Share of the DLTs that occur in the
#'   second half of the assessment window under the Weibull distribution.
#'   Defaults to 0.5, the setting of Yuan et al. (2018) and Lin and Yuan (2020).
#'
#' @param prior_weights
#'   Numeric vector of three non-negative values. Prior probabilities that a DLT
#'   occurs in the first, second and last third of the assessment window, used
#'   to weight the follow-up of the pending patients in STFT (Yuan et al., 2018,
#'   Supplementary Appendix D) and ESS (Lin and Yuan, 2020, Supplementary S1).
#'   They are rescaled to sum to one. Defaults to equal weights, which give the
#'   unweighted STFT. The decision table does not change with the weights; only
#'   the value compared with its boundaries does. This describes the design and
#'   is separate from \code{dlt_time} and \code{late_fraction}, which describe
#'   how the simulated DLTs actually occur.
#'
#' @param max_pending_ratio
#'   Numeric scalar or \code{NULL}. Suspension rule on the proportion of pending
#'   patients. See \code{\link{tite_boin_decision_table}}.
#'
#' @param min_completed
#'   Integer scalar or \code{NULL}. Number of completed assessments required
#'   for escalation. See \code{\link{tite_boin_decision_table}}.
#'
#' @param n_trials
#'   Integer scalar. Number of trials to simulate. Defaults to 10000.
#'
#' @param start_dose
#'   Integer scalar. Dose level for the first cohort. Defaults to 1.
#'
#' @param n_earlystop
#'   Integer scalar. The trial stops once this many patients have been treated at
#'   the current dose and the design would stay there. Defaults to 18.
#'
#' @param p_saf
#'   Numeric scalar. Highest DLT probability deemed subtherapeutic.
#'   Defaults to \code{0.6 * target}.
#'
#' @param p_tox
#'   Numeric scalar. Lowest DLT probability deemed overly toxic.
#'   Defaults to \code{1.4 * target}.
#'
#' @param cutoff_eli
#'   Numeric scalar. Posterior probability cutoff for dose elimination.
#'   Defaults to 0.95.
#'
#' @param extrasafe
#'   Logical scalar. Apply the stricter safety stopping rule at the lowest dose.
#'   Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is
#'   relaxed for the safety stopping rule. Defaults to 0.05.
#'
#' @param n_earlystop_rule
#'   Character scalar, either \code{"with_stay"} or \code{"simple"}. See
#'   \code{\link{boin_simulate}}.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. When \code{NULL} the current
#'   state of the random number generator is used. The state of the calling
#'   session is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{tite_boin_trials}, which is a list with components
#'   \item{n_pts}{Integer matrix, \code{n_trials} by number of doses, of patients treated.}
#'   \item{n_tox}{Integer matrix of the same shape, of DLTs, counted once every
#'     patient has completed the assessment.}
#'   \item{eliminated}{Logical matrix of the same shape, marking doses eliminated during the trial.}
#'   \item{cohorts_used}{Integer vector of cohorts treated in each trial.}
#'   \item{stop_reason}{Character vector giving why each trial ended, with the
#'     values of \code{\link{boin_simulate}}.}
#'   \item{duration}{Numeric vector. Time from the first arrival until every
#'     enrolled patient has completed the assessment.}
#'   \item{n_suspensions}{Integer vector. Number of decisions at which accrual
#'     had to be suspended.}
#'   \item{time_suspended}{Numeric vector. Total time with accrual suspended.}
#'   \item{boundary}{The \code{boin_boundary} object used.}
#'   \item{settings}{List of the design parameters.}
#'
#' @details
#'   Patients arrive one at a time. The dose for a cohort is decided when its
#'   first patient arrives, from the data available at that moment, and the other
#'   patients of the cohort receive the same dose. The decision follows
#'   \code{\link{tite_boin_decision_table}}. A patient with a DLT completes the
#'   assessment when the DLT occurs, and a patient without one at the end of the
#'   window. When the decision is to suspend accrual, the arriving patient waits
#'   until the next pending patient at the current dose completes the assessment,
#'   and the decision is taken again at that moment. Later patients arrive after
#'   the one who waited.
#'
#'   At every decision the rules are applied in the order used by
#'   \code{\link{boin_simulate}}: elimination at the current dose and the extra
#'   safety rule, early stopping on \code{n_earlystop}, then the dose
#'   transition. Because DLTs can be observed after a dose has been left, the
#'   elimination rule is also checked at the other doses; a dose that meets it is
#'   eliminated together with all higher doses, and the trial moves below it.
#'   Elimination counts every treated patient, the pending ones as without DLT.
#'
#'   One uniform variate is drawn from R's random number stream per patient, in
#'   enrollment order and in the cohort blocks of \code{\link{boin_simulate}}.
#'   It decides whether the patient has a DLT and, under the Weibull
#'   distribution, gives the time to DLT through the inverse distribution
#'   function, which falls inside the window exactly when there is a DLT. The
#'   times between arrivals come from a second generator (xoshiro256**) seeded
#'   with \code{seed}, or with one draw from R's stream when \code{seed} is
#'   \code{NULL}. As a result, whenever no patient is pending at any decision,
#'   for example with \code{accrual = "fixed"} and \code{window} shorter than
#'   \code{1 / accrual_rate}, the trials are identical to those of
#'   \code{\link{boin_simulate}} with the same seed.
#'
#'   The follow-up of a pending patient counts towards STFT as the share of the
#'   window already covered, or, with \code{prior_weights}, as the prior
#'   probability that a DLT would have occurred by then. Both articles evaluate
#'   their designs with equal weights, the default.
#'
#'   The titration phase and the \code{stay_on_1_of_3} option of
#'   \code{\link{boin_simulate}} are not available, as neither article defines
#'   them for the time-to-event design.
#'
#' @references
#'   Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K. E. (2018). Time-to-Event
#'   Bayesian Optimal Interval Design to Accelerate Phase I Trials. Clinical
#'   Cancer Research, 24(20), 4921-4930.
#'
#'   Lin, R. and Yuan, Y. (2020). Time-to-Event Model-Assisted Designs for
#'   Dose-Finding Trials with Delayed Toxicity. Biostatistics, 21(4), 807-824.
#'
#' @examples
#' # A three month window with two patients a month
#' trials <- tite_boin_simulate(
#'   target = 0.30,
#'   p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   window = 3,
#'   accrual_rate = 2,
#'   n_trials = 200,
#'   seed = 123
#' )
#' trials
#'
#' summary(trials$duration)
#'
#' @seealso \code{\link{tite_boin_decision_table}}, \code{\link{boin_simulate}}
#'
#' @export
tite_boin_simulate <- function(target, p_true, n_cohort, cohort_size,
                               window, accrual_rate,
                               method = c("imputation", "ess"),
                               accrual = c("exponential", "uniform", "fixed"),
                               dlt_time = c("weibull", "uniform"),
                               late_fraction = 0.5,
                               prior_weights = c(1, 1, 1) / 3,
                               max_pending_ratio = NULL, min_completed = NULL,
                               n_trials = 10000, start_dose = 1, n_earlystop = 18,
                               p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                               extrasafe = FALSE, offset = 0.05,
                               n_earlystop_rule = c("with_stay", "simple"),
                               seed = 123) {

  method <- match.arg(method)
  accrual <- match.arg(accrual)
  dlt_time <- match.arg(dlt_time)
  n_earlystop_rule <- match.arg(n_earlystop_rule)
  rules <- tite_rules(method, max_pending_ratio, min_completed)

  check_p_true(p_true)
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(n_earlystop, "n_earlystop", 1L)
  check_count(start_dose, "start_dose", 1L)
  check_flag(extrasafe, "extrasafe")
  for (arg in c("window", "accrual_rate")) {
    value <- get(arg)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) {
      stop("'", arg, "' must be a single positive number", call. = FALSE)
    }
  }
  check_scalar_prob(late_fraction, "late_fraction")
  if (!is.numeric(prior_weights) || length(prior_weights) != 3L ||
      any(!is.finite(prior_weights)) || any(prior_weights < 0) ||
      sum(prior_weights) <= 0) {
    stop("'prior_weights' must be three non-negative numbers with a positive sum",
         call. = FALSE)
  }
  prior_weights <- prior_weights / sum(prior_weights)
  # Equal weights take the unweighted computation of STFT in the engine.
  weighted <- any(abs(prior_weights - 1 / 3) > 1e-12)
  if (dlt_time == "weibull" && any(p_true >= 1)) {
    stop("'p_true' must be below 1 when 'dlt_time' is \"weibull\"", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         !is.finite(seed) || seed != round(seed))) {
    stop("'seed' must be a single whole number or NULL", call. = FALSE)
  }
  if (n_earlystop <= 6) {
    warning("'n_earlystop' is low; values between 9 and 18 are recommended",
            call. = FALSE)
  }

  n_doses <- length(p_true)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of doses in 'p_true'", call. = FALSE)
  }

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)

  bound <- boin_boundary(target, max_total_pts, p_saf = p_saf, p_tox = p_tox,
                         cutoff_eli = cutoff_eli, extrasafe = extrasafe,
                         offset = offset)

  # The engine encodes a missing elimination boundary as zero.
  b_elim <- bound$b_elim
  b_elim[is.na(b_elim)] <- 0L

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
    stream_seed <- seed
  } else {
    stream_seed <- floor(stats::runif(1L) * 2147483647)
  }

  res <- tite_boin_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    n_earlystop = as.integer(n_earlystop),
    early_stop_simple = identical(n_earlystop_rule, "simple"),
    extrasafe = extrasafe,
    target = as.numeric(target),
    cutoff_eli = as.numeric(cutoff_eli),
    offset = as.numeric(offset),
    b_esc = bound$b_esc,
    b_deesc = bound$b_deesc,
    b_elim = as.integer(b_elim),
    max_total_pts = as.integer(max_total_pts),
    method = if (method == "imputation") 0L else 1L,
    lambda_e = bound$lambda_e,
    lambda_d = bound$lambda_d,
    max_pending_ratio = as.numeric(rules$max_pending_ratio),
    min_completed = rules$min_completed,
    window = as.numeric(window),
    accrual = match(accrual, c("exponential", "uniform", "fixed")) - 1L,
    accrual_rate = as.numeric(accrual_rate),
    dlt_time = match(dlt_time, c("weibull", "uniform")) - 1L,
    late_fraction = as.numeric(late_fraction),
    weighted = weighted,
    prior_weights = as.numeric(prior_weights),
    stream_seed = as.numeric(stream_seed)
  )

  reasons <- c("lowest_dose_eliminated", "lowest_dose_too_toxic",
               "n_earlystop", "max_sample_size", "max_cohorts")

  dose_names <- paste0("DL", seq_len(n_doses))
  colnames(res$n_pts) <- dose_names
  colnames(res$n_tox) <- dose_names
  colnames(res$eliminated) <- dose_names

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$n_tox,
      eliminated = res$eliminated,
      cohorts_used = res$cohorts_used,
      stop_reason = reasons[res$stop_code + 1L],
      duration = res$duration,
      n_suspensions = res$n_suspensions,
      time_suspended = res$time_suspended,
      boundary = bound,
      settings = list(
        target = target,
        p_true = p_true,
        n_cohort = n_cohort,
        cohort_size = cohort_size_vec,
        window = window,
        accrual_rate = accrual_rate,
        method = method,
        accrual = accrual,
        dlt_time = dlt_time,
        late_fraction = late_fraction,
        prior_weights = prior_weights,
        max_pending_ratio = rules$max_pending_ratio,
        min_completed = rules$min_completed,
        n_trials = n_trials,
        start_dose = start_dose,
        n_earlystop = n_earlystop,
        n_earlystop_rule = n_earlystop_rule,
        p_saf = bound$p_saf,
        p_tox = bound$p_tox,
        cutoff_eli = cutoff_eli,
        extrasafe = extrasafe,
        offset = offset,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "tite_boin_trials"
  )
}
