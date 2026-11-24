#' Simulate TITE-gBOIN-ET Trials
#'
#' @description
#'   Run the TITE-gBOIN-ET dose-finding algorithm of Takeda, Yamaguchi,
#'   Taguri and Morita (2023) over many simulated trials, for toxicity and
#'   efficacy outcomes with several categories that each carry a score, in
#'   which patients arrive over time and the outcomes are assessed over
#'   windows, and return the patients, their categories and the total scores
#'   at every dose. This is the simulation engine behind
#'   \code{\link{sim_tite_gboin_et}}, exposed for users who want the raw
#'   trial data.
#'
#' @inheritParams gboin_et_simulate
#'
#' @param window
#'   Numeric scalar. Length of the toxicity assessment window.
#'
#' @param accrual_rate
#'   Numeric scalar. Average number of patients arriving per unit of time, in
#'   the same unit as \code{window}. The \code{accrual} argument of the
#'   'boinet' package is the average time between arrivals, its inverse.
#'
#' @param resp_window
#'   Numeric scalar. Length of the efficacy assessment window. Defaults to
#'   \code{window}.
#'
#' @param accrual
#'   Character scalar. Distribution of the times between arrivals:
#'   \code{"uniform"} (the default, as in 'boinet') on
#'   \code{(0, 2 / accrual_rate)}, \code{"exponential"} for a Poisson process,
#'   or \code{"fixed"} for arrivals exactly \code{1 / accrual_rate} apart.
#'
#' @param event_time
#'   Character scalar. Distribution of the times to the toxicity event and to
#'   the efficacy event: \code{"weibull"} (the default), with the share of
#'   the events in the second half of the window set by \code{late_fraction}
#'   and \code{resp_late_fraction}, or \code{"uniform"} on
#'   \code{(0, window / p)} for an event of probability \code{p}. With
#'   \code{"weibull"} the first category of both outcomes, the one without
#'   an event, must have a positive probability at every dose.
#'
#' @param late_fraction
#'   Numeric scalar strictly between 0 and 1. Share of the toxicity events
#'   that occur in the second half of the toxicity window under the Weibull
#'   distribution. Defaults to 0.5.
#'
#' @param resp_late_fraction
#'   Numeric scalar strictly between 0 and 1. The same for the efficacy
#'   events. Defaults to 0.5.
#'
#' @param resp_cor
#'   Numeric scalar strictly between -1 and 1. Correlation of the Gaussian
#'   copula that links the time to the toxicity event and the time to the
#'   efficacy event of a patient. Defaults to 0.2, the default of the
#'   'boinet' package.
#'
#' @param min_completed_share
#'   Numeric scalar from 0 to 1. Share of the patients at the current dose who
#'   must have finished both assessments before the next decision; the
#'   decision waits until then. Defaults to 0.5, the rule of Takeda et al.
#'   (2023) and of the 'boinet' package.
#'
#' @return
#'   An object of class \code{boin_et_trials}, as returned by
#'   \code{\link{gboin_et_simulate}}, with \code{design} set to
#'   \code{"TITE-gBOIN-ET"} and one more component
#'   \item{duration}{Numeric vector of the time of the last decision of each
#'     trial, which is the trial duration when all cohorts are enrolled.}
#'   The categories and the total scores are those observed at the last
#'   decision; a patient still under assessment is counted in the first
#'   category.
#'
#' @details
#'   The TITE-gBOIN-ET design is the gBOIN-ET design
#'   (\code{\link{gboin_et_simulate}}) for outcomes that are observed late:
#'   the next cohort is enrolled while earlier patients are still assessed.
#'   The first patient of a cohort arrives at the time of the previous
#'   decision, and the next decision is due when the next patient arrives;
#'   after the last cohort it is due once every patient has been followed
#'   over both windows. While fewer than a share \code{min_completed_share}
#'   of the patients at the current dose have finished both assessments, the
#'   decision waits for the next assessment to end.
#'
#'   At a decision, the mean toxicity score of a dose is the total of the
#'   normalized scores observed so far divided by the effective number of
#'   patients: the patients whose toxicity assessment has ended plus, for
#'   every patient still under assessment, the fraction of the window
#'   followed so far (Takeda et al., 2023). The mean efficacy score is
#'   formed in the same way, and the decision rule of gBOIN-ET is applied to
#'   the two means of the current dose. The admissibility of a dose uses the
#'   total scores observed so far and the number of patients treated at the
#'   dose, so that a patient still under assessment counts as a patient
#'   without an event, and the moves between doses follow
#'   \code{\link{boin_et_simulate}}, including the changes from 'boinet'
#'   unless \code{boinet_compatible} is \code{TRUE}. These are the rules of
#'   \code{boinet::tite.gboinet()}.
#'
#'   A patient has a toxicity event when the time to the event falls within
#'   the window, which happens with the probability of the categories above
#'   the first, and an efficacy event likewise. The two event times come
#'   from their marginal distributions linked by a Gaussian copula, and the
#'   category of a patient with an event is drawn from the probabilities of
#'   the categories above the first, independently of everything else; it is
#'   observed at the time of the event, and a patient without an event is
#'   known to have none when the window ends. Each patient draws one uniform
#'   variate for the toxicity event, one for the efficacy event, and one
#'   more for each outcome with more than two categories. With two
#'   categories of scores 0 and 1 for both outcomes the trials are, for the
#'   same seed, those of \code{\link{tite_boin_et_simulate}}; with arrivals
#'   a fixed time apart that is longer than both windows, they are those of
#'   \code{\link{gboin_et_simulate}}. The 'boinet' package draws its random
#'   numbers differently, so the trials agree with it in distribution, not
#'   one by one.
#'
#'   The program of the supporting information of Takeda et al. (2023)
#'   differs from these rules in two points. Where one dose among the
#'   current dose and its neighbors has the highest mean efficacy score, the
#'   program passes that dose to \code{sample()}, which then draws from all
#'   doses up to it, so that the trial can move to any lower dose; and the
#'   program excludes doses with the effective numbers of patients in the
#'   place of the numbers of patients. Neither is described in the article,
#'   and 'boinet' has neither. The first changes the published operating
#'   characteristics substantially; see
#'   \code{inst/reproduce/reproduce_Takeda2023.R}.
#'
#'   The mean score of a dose is computed as the total of the scores as
#'   given, divided by the largest score times the effective number of
#'   patients, as in \code{\link{gboin_et_simulate}}; when no patient of the
#'   dose is still under assessment, doses whose mean scores are equal as
#'   fractions tie on every platform.
#'
#'   Possible values of \code{stop_reason} are those of
#'   \code{\link{boin_et_simulate}}.
#'
#' @references
#'   Takeda, K., Yamaguchi, Y., Taguri, M. and Morita, S. (2023).
#'   TITE-gBOIN-ET: Time-to-event generalized Bayesian optimal interval
#'   design to accelerate dose-finding accounting for ordinal graded efficacy
#'   and toxicity outcomes. Biometrical Journal, 65(7), 2200265.
#'   \doi{10.1002/bimj.202200265}
#'
#'   Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
#'   Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
#'   and Toxicity Outcomes. The R Journal, 17(2), 15-34.
#'   \doi{10.32614/RJ-2025-012}
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
#' trials <- tite_gboin_et_simulate(
#'   target = 0.47 / 1.5, target_resp = 1.75 / 3,
#'   p_grade = p_grade, p_resp_grade = p_resp_grade,
#'   score = c(0, 0.5, 1, 1.5), score_resp = c(0, 0.25, 1, 3),
#'   n_cohort = 12, cohort_size = 3, window = 30, resp_window = 60,
#'   accrual_rate = 1 / 10, n_trials = 200, seed = 123
#' )
#' trials
#' summary(trials$duration)
#'
#' @seealso \code{\link{sim_tite_gboin_et}},
#'   \code{\link{gboin_et_select_obd}}, \code{\link{gboin_et_simulate}},
#'   \code{\link{tite_boin_et_simulate}}
#'
#' @export
tite_gboin_et_simulate <- function(target, target_resp, p_grade,
                                   p_resp_grade, score, score_resp, n_cohort,
                                   cohort_size, window, accrual_rate,
                                   resp_window = window,
                                   accrual = c("uniform", "exponential",
                                               "fixed"),
                                   event_time = c("weibull", "uniform"),
                                   late_fraction = 0.5,
                                   resp_late_fraction = 0.5, resp_cor = 0.2,
                                   min_completed_share = 0.5,
                                   n_trials = 10000, start_dose = 1,
                                   n_earlystop = NULL, p_saf = NULL,
                                   p_tox = NULL, p_saf_resp = NULL,
                                   lambda1 = NULL, lambda2 = NULL,
                                   eta1 = NULL, cutoff_eli = 0.95,
                                   cutoff_eli_resp = 0.99, p_eli = NULL,
                                   p_eli_resp = NULL,
                                   boinet_compatible = FALSE,
                                   prior = rep(1 / 6, 6), seed = 123) {

  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  check_flag(boinet_compatible, "boinet_compatible")

  tox <- gboin_check_grades(p_grade, score)
  n_doses <- ncol(p_grade)
  resp <- gboin_et_check_resp(p_resp_grade, score_resp, n_doses)
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(start_dose, "start_dose", 1L)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of doses in 'p_grade'",
         call. = FALSE)
  }
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
  if (event_time == "weibull" &&
      (any(p_grade[1L, ] <= 0) || any(p_resp_grade[1L, ] <= 0))) {
    stop("the first rows of 'p_grade' and 'p_resp_grade' must be positive ",
         "with event_time = \"weibull\"", call. = FALSE)
  }
  # Probabilities of an event: those of the categories above the first.
  p_event <- tox$cum_prob[2L, ]
  p_event_resp <- resp$cum_prob[2L, ]
  par_t <- boin_et_event_params(p_event, window, late_fraction, event_time,
                                "p_grade")
  par_e <- boin_et_event_params(p_event_resp, resp_window,
                                resp_late_fraction, event_time,
                                "p_resp_grade")

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)
  if (is.null(n_earlystop)) n_earlystop <- max_total_pts
  check_count(n_earlystop, "n_earlystop", 1L)

  bd <- gboin_et_boundary_arg(target, target_resp, p_saf, p_tox, p_saf_resp,
                              prior, lambda1, lambda2, eta1)
  eli <- boin_et_check_eli(p_eli, p_eli_resp, bd)

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- tite_gboin_et_simulate_cpp(
    n_trials = as.integer(n_trials),
    cum_tox = tox$cum_prob,
    cum_resp = resp$cum_prob,
    score_tox = as.numeric(score),
    score_resp = as.numeric(score_resp),
    resp_cor = as.numeric(resp_cor),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    n_earlystop = as.integer(n_earlystop),
    lambda1 = as.numeric(bd$lambda1),
    lambda2 = as.numeric(bd$lambda2),
    eta1 = as.numeric(bd$eta1),
    p_eli = as.numeric(eli$p_eli),
    p_eli_resp = as.numeric(eli$p_eli_resp),
    cutoff_eli = as.numeric(cutoff_eli),
    cutoff_eli_resp = as.numeric(cutoff_eli_resp),
    boinet_compatible = as.integer(boinet_compatible),
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

  reasons <- boin_et_stop_reasons()

  dose_names <- paste0("DL", seq_len(n_doses))
  # Row names are used when every row has one; rbind() of named and unnamed
  # vectors leaves some of them empty.
  grade_names <- rownames(p_grade)
  if (is.null(grade_names) || any(!nzchar(grade_names))) {
    grade_names <- paste0("G", seq_len(nrow(p_grade)))
  }
  resp_names <- rownames(p_resp_grade)
  if (is.null(resp_names) || any(!nzchar(resp_names))) {
    resp_names <- paste0("R", seq_len(nrow(p_resp_grade)))
  }
  colnames(res$n_pts) <- dose_names
  colnames(res$tox_sum) <- dose_names
  colnames(res$resp_sum) <- dose_names
  dimnames(res$n_tox_cat) <- list(NULL, dose_names, grade_names)
  dimnames(res$n_resp_cat) <- list(NULL, dose_names, resp_names)

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$tox_sum,
      n_resp = res$resp_sum,
      n_grade = res$n_tox_cat,
      n_resp_grade = res$n_resp_cat,
      cohorts_used = res$cohorts_used,
      stop_reason = reasons[res$stop_code + 1L],
      early_stop = res$stop_code <= 1L,
      duration = res$duration,
      design = "TITE-gBOIN-ET",
      boundary = bd,
      settings = list(
        target = bd$target,
        target_resp = bd$target_resp,
        p_true = tox$p_true,
        p_resp = resp$p_resp,
        p_grade = p_grade,
        p_resp_grade = p_resp_grade,
        score = score,
        score_resp = score_resp,
        score_norm = tox$score_norm,
        score_resp_norm = resp$score_norm,
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
