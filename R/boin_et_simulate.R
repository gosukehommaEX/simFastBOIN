#' Simulate BOIN-ET Trials
#'
#' @description
#'   Run the BOIN-ET dose-finding algorithm of Takeda, Taguri and Morita
#'   (2018) over many simulated trials and return the patient, DLT and
#'   response counts at every dose. This is the simulation engine behind
#'   \code{\link{sim_boin_et}}, exposed for users who want the raw trial data.
#'
#' @inheritParams boin_et_decision_table
#'
#' @param p_true
#'   Numeric vector. True DLT probability at each dose level, in increasing
#'   dose order.
#'
#' @param p_resp
#'   Numeric vector of the same length as \code{p_true}. True response
#'   probability at each dose level.
#'
#' @param n_cohort
#'   Integer scalar. Number of cohorts in a trial.
#'
#' @param cohort_size
#'   Integer scalar or vector. Number of patients per cohort. A scalar is used
#'   for every cohort. A vector shorter than \code{n_cohort} is padded with its
#'   last element and a longer one is truncated.
#'
#' @param n_trials
#'   Integer scalar. Number of trials to simulate. Defaults to 10000.
#'
#' @param start_dose
#'   Integer scalar. Dose level for the first cohort. Defaults to 1.
#'
#' @param n_earlystop
#'   Integer scalar or \code{NULL}. The trial stops once this many patients
#'   have been treated at any dose. Defaults to \code{NULL}, the maximum sample
#'   size, which is the default of the 'boinet' package.
#'
#' @param resp_cor
#'   Numeric scalar strictly between -1 and 1. Correlation between the DLT and
#'   the response of a patient. Defaults to 0.2, the default of the 'boinet'
#'   package.
#'
#' @param resp_cor_type
#'   Character scalar. How \code{resp_cor} links the two outcomes:
#'   \code{"copula"} (the default) uses it as the correlation of a Gaussian
#'   copula, which is how the 'boinet' package generates its data;
#'   \code{"binary"} uses it as the correlation of the two binary outcomes, as
#'   in the simulations of Takeda et al. (2018).
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. When \code{NULL} the current
#'   state of the random number generator is used. The state of the calling
#'   session is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{boin_et_trials}, which is a list with components
#'   \item{n_pts}{Integer matrix, \code{n_trials} by number of doses, of
#'     patients treated.}
#'   \item{n_tox}{Integer matrix of the same shape, of DLTs observed.}
#'   \item{n_resp}{Integer matrix of the same shape, of responses observed.}
#'   \item{cohorts_used}{Integer vector of cohorts enrolled in each trial.}
#'   \item{stop_reason}{Character vector giving why each trial ended.}
#'   \item{early_stop}{Logical vector, \code{TRUE} for the trials stopped
#'     early, which select no OBD.}
#'   \item{boundary}{The \code{\link{boin_et_boundary}} object used.}
#'   \item{settings}{List of the design parameters.}
#'
#' @details
#'   Every cohort is treated at the current dose and its outcomes are observed
#'   before the next decision, which follows \code{\link{boin_et_decision_table}}.
#'   After each cohort the admissibility of every dose is re-derived from its
#'   data. The trial stops early when no dose is admissible, or when the dose
#'   should be lowered and no admissible dose is found at or below the
#'   intended one. When the intended dose is not admissible, an escalation
#'   moves to the lowest admissible dose above it (the next cohort stays when
#'   there is none), and a stay or a de-escalation moves to the highest
#'   admissible dose below it. An escalation from the highest dose stays
#'   there. The trial also stops, without being counted as an early stop,
#'   once a dose has \code{n_earlystop} patients. These are the rules of
#'   \code{boinet::boinet()}.
#'
#'   Each patient draws two uniform variates in enrollment order, the first for
#'   the DLT and the second for the response, and a tie under the rule for the
#'   neighboring doses draws one more. The 'boinet' package draws its random
#'   numbers differently (event times from a copula, and arrival times), so the
#'   trials agree with it in distribution, not one by one.
#'
#'   Possible values of \code{stop_reason} are \code{"no_admissible_dose"},
#'   \code{"no_lower_admissible_dose"}, \code{"n_earlystop"} and
#'   \code{"max_cohorts"}; the first two are early stops.
#'
#' @references
#'   Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
#'   interval design for dose finding based on both efficacy and toxicity
#'   outcomes. Pharmaceutical Statistics, 17(4), 383-395.
#'   \doi{10.1002/pst.1864}
#'
#' @examples
#' trials <- boin_et_simulate(
#'   target = 0.30, target_resp = 0.60,
#'   p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
#'   p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
#'   n_cohort = 12, cohort_size = 3, n_trials = 200, seed = 123
#' )
#' trials
#'
#' head(trials$n_resp)
#' table(trials$stop_reason)
#'
#' @seealso \code{\link{sim_boin_et}}, \code{\link{boin_et_select_obd}}
#'
#' @export
boin_et_simulate <- function(target, target_resp, p_true, p_resp, n_cohort,
                             cohort_size, n_trials = 10000, start_dose = 1,
                             n_earlystop = NULL, p_saf = NULL, p_tox = NULL,
                             p_saf_resp = NULL, cutoff_eli = 0.95,
                             cutoff_eli_resp = 0.99, p_eli = NULL,
                             p_eli_resp = NULL, resp_cor = 0.2,
                             resp_cor_type = c("copula", "binary"),
                             prior = rep(1 / 6, 6), seed = 123) {

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
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  check_scalar_prob(resp_cor, "resp_cor", lower = -1, upper = 1)

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)
  if (is.null(n_earlystop)) n_earlystop <- max_total_pts
  check_count(n_earlystop, "n_earlystop", 1L)

  bd <- boin_et_boundary(target, target_resp, p_saf = p_saf, p_tox = p_tox,
                         p_saf_resp = p_saf_resp, prior = prior)
  codes <- boin_et_decision_codes(bd$lambda1, bd$lambda2, bd$eta1,
                                  max_total_pts)
  eli <- boin_et_check_eli(p_eli, p_eli_resp, bd)
  adm <- boin_et_admissible_table(eli$p_eli, eli$p_eli_resp, cutoff_eli,
                                  cutoff_eli_resp, max_total_pts)

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

  res <- boin_et_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    p_resp = as.numeric(p_resp),
    q0 = as.numeric(link$q0),
    q1 = as.numeric(link$q1),
    cor_type = if (resp_cor_type == "binary") 1L else 0L,
    resp_cor = as.numeric(resp_cor),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    n_earlystop = as.integer(n_earlystop),
    dec = as.integer(codes),
    tox_ok = as.integer(adm$tox_ok),
    resp_ok = as.integer(adm$resp_ok),
    max_n = as.integer(max_total_pts)
  )

  reasons <- c("no_admissible_dose", "no_lower_admissible_dose",
               "n_earlystop", "max_cohorts")

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
      stop_reason = reasons[res$stop_code + 1L],
      early_stop = res$stop_code <= 1L,
      boundary = bd,
      settings = list(
        target = bd$target,
        target_resp = bd$target_resp,
        p_true = p_true,
        p_resp = p_resp,
        n_cohort = n_cohort,
        cohort_size = cohort_size_vec,
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
        resp_cor_type = resp_cor_type,
        prior = prior,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "boin_et_trials"
  )
}
