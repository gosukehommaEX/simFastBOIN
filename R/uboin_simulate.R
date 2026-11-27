#' Simulate U-BOIN Trials
#'
#' @description
#'   Simulate trials of the utility-based BOIN (U-BOIN) phase I/II design of
#'   Zhou, Lee and Yuan (2019), which looks for the optimal biological dose
#'   (OBD) from a binary toxicity and a binary efficacy outcome observed before
#'   the next cohort is treated. This is the simulation engine behind
#'   \code{\link{sim_uboin}}, exposed for users who want the raw trial data.
#'
#' @param tox_upper
#'   Numeric scalar. Highest acceptable DLT probability, for example 0.30.
#'
#' @param resp_lower
#'   Numeric scalar. Lowest acceptable response probability, for example 0.20.
#'
#' @param p_true
#'   Numeric vector or \code{NULL}. True DLT probability of each dose.
#'
#' @param p_resp
#'   Numeric vector or \code{NULL}. True response probability of each dose.
#'
#' @param n_cohort
#'   Integer scalar. Number of cohorts in a trial; the maximum sample size is
#'   \code{n_cohort * cohort_size}.
#'
#' @param cohort_size
#'   Integer scalar. Number of patients per cohort.
#'
#' @param n_trials
#'   Integer scalar. Number of trials to simulate. Defaults to 10000.
#'
#' @param target
#'   Numeric scalar or \code{NULL}. Target DLT probability of the BOIN design in
#'   stage I. Defaults to \code{NULL}, which uses \code{tox_upper - 0.05}.
#'
#' @param resp_assoc
#'   Numeric scalar. Association parameter of the Gumbel model that joins the
#'   DLT and the response of a patient; 0 makes them independent. Defaults to
#'   0.
#'
#' @param p_outcome
#'   Numeric matrix or \code{NULL}. Probabilities of the four outcomes of a
#'   patient at each dose, one row per dose, in the columns (no response, DLT),
#'   (no response, no DLT), (response, DLT) and (response, no DLT). When given,
#'   \code{p_true}, \code{p_resp} and \code{resp_assoc} must be \code{NULL} or
#'   are not used. Defaults to \code{NULL}.
#'
#' @param utility
#'   Numeric vector of length four. Utility of the four outcomes, in the order
#'   of \code{p_outcome}. Defaults to \code{c(0, 30, 50, 100)}, example 1 of the
#'   article.
#'
#' @param n_stage1
#'   Integer scalar. Stage II starts once a dose has this many patients (s1 of
#'   the article). Defaults to 12.
#'
#' @param n_stop
#'   Integer scalar or \code{NULL}. The trial ends once a dose has this many
#'   patients (s2 of the article). Defaults to \code{NULL}, the maximum sample
#'   size.
#'
#' @param start_dose
#'   Integer scalar. Dose level of the first cohort. Defaults to 1.
#'
#' @param p_saf
#'   Numeric scalar. Highest DLT probability deemed subtherapeutic, for the
#'   boundaries of stage I. Defaults to \code{0.6 * target}.
#'
#' @param p_tox
#'   Numeric scalar. Lowest DLT probability deemed overly toxic, for the
#'   boundaries of stage I. Defaults to \code{1.4 * target}.
#'
#' @param cutoff_eli
#'   Numeric scalar. A dose is toxic when the posterior probability that its
#'   DLT probability exceeds \code{tox_upper} is above this value. Defaults to
#'   0.95.
#'
#' @param cutoff_eli_resp
#'   Numeric scalar. A dose is futile when the posterior probability that its
#'   response probability is below \code{resp_lower} is above this value.
#'   Defaults to 0.90.
#'
#' @param prior_tox
#'   Numeric vector of length two. Beta prior of the DLT probability in the
#'   toxicity rule, which eliminates doses in stage I as well as in stage
#'   II. Defaults to \code{c(1, 1)}, the uniform prior of the BOIN design
#'   that stage I of the article uses.
#'
#' @param prior_resp
#'   Numeric vector of length two. Beta prior of the response probability in
#'   the futility rule. Defaults to \code{c(1, 1)}.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. The state of the calling
#'   session is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{uboin_trials}, which is a list with components
#'   \item{n_outcome}{Integer array with dimensions trial by dose by outcome:
#'     number of patients with each outcome.}
#'   \item{n_pts, n_tox, n_resp}{Integer matrices with one row per trial and one
#'     column per dose: patients treated, DLTs and responses.}
#'   \item{eliminated}{Logical matrix: doses eliminated in stage I.}
#'   \item{cohorts_used}{Integer vector: cohorts treated in each trial.}
#'   \item{stage2}{Logical vector: whether each trial reached stage II.}
#'   \item{stop_reason}{Character vector: why each trial ended.}
#'   \item{boundary}{List of the boundaries of stage I.}
#'   \item{design}{\code{"U-BOIN"}.}
#'   \item{settings}{List of the design parameters.}
#'
#' @details
#'   Stage I is the BOIN design with the target \code{target}. After each
#'   cohort, a dose that is toxic and has at least three patients is
#'   eliminated together with every higher dose, and the trial stops when the
#'   lowest dose is eliminated; the next cohort escalates (unless the next dose
#'   is eliminated), stays or de-escalates by the boundaries. Once a dose has
#'   \code{n_stage1} patients, stage II starts with the next cohort:
#'   \itemize{
#'     \item (B1) If the DLT rate observed at the highest dose tried is at most
#'       the escalation boundary, that dose is not the highest dose and the dose
#'       above it is not eliminated, the next cohort is treated at the dose
#'       above it.
#'     \item (B2) Otherwise the next cohort is treated at the admissible dose
#'       with the largest posterior mean utility, the lowest among ties. A dose
#'       is admissible when it has been tried, is not eliminated, and is
#'       neither toxic nor futile. The trial stops, without an OBD, when no dose
#'       is admissible.
#'   }
#'   The trial ends at the maximum sample size or when a dose has
#'   \code{n_stop} patients. The posterior mean utility of a dose is
#'   \code{sum(utility * (1/4 + n_k)) / (1 + n)}, the utilities weighted by
#'   the posterior means of the outcome probabilities under the
#'   Dirichlet(1/4, 1/4, 1/4, 1/4) prior of the article. Utilities are
#'   compared after rounding to ten decimals.
#'
#'   The article does not give the priors of the toxicity and futility rules.
#'   With Beta(1, 1) priors the simulations of the article are reproduced
#'   (\code{inst/reproduce/reproduce_Zhou2019.R}); the decision tables of its
#'   supporting information are reproduced with \code{prior_resp = c(0.05,
#'   0.05)}. The article does not say either whether the admissible doses of
#'   rule B2 include doses not yet tried; they do not here, which the
#'   simulations of the article support. The stage I elimination uses the
#'   toxicity rule, as the article describes. A dose that becomes toxic in
#'   stage II does not exclude the doses above it.
#'
#'   One uniform random variate is drawn per patient, in enrollment order, and
#'   the outcome is the first whose cumulative probability exceeds it.
#'
#'   Only the pick-the-winner allocation of the article is provided; the
#'   adaptive and equal randomization and the version for efficacy observed
#'   late (Section 2.5 of the article) are not.
#'
#'   Values of \code{stop_reason} are \code{"max_sample_size"},
#'   \code{"n_stop"}, \code{"lowest_dose_eliminated"} and
#'   \code{"no_admissible_dose"}.
#'
#' @references
#'   Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian optimal
#'   interval (U-BOIN) phase I/II design to identify the optimal biological
#'   dose for targeted and immune therapies. Statistics in Medicine, 38(28),
#'   S5299-S5316. \doi{10.1002/sim.8361}
#'
#' @examples
#' trials <- uboin_simulate(
#'   tox_upper = 0.30,
#'   resp_lower = 0.20,
#'   p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
#'   p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65),
#'   n_cohort = 18,
#'   cohort_size = 3,
#'   resp_assoc = 0.2,
#'   n_trials = 200,
#'   seed = 123
#' )
#' trials
#'
#' @seealso \code{\link{sim_uboin}}, \code{\link{uboin_select_obd}},
#'   \code{\link{uboin_utility_table}}
#'
#' @export
uboin_simulate <- function(tox_upper, resp_lower, p_true = NULL,
                           p_resp = NULL, n_cohort, cohort_size,
                           n_trials = 10000, target = NULL, resp_assoc = 0,
                           p_outcome = NULL, utility = c(0, 30, 50, 100),
                           n_stage1 = 12, n_stop = NULL, start_dose = 1,
                           p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                           cutoff_eli_resp = 0.90, prior_tox = c(1, 1),
                           prior_resp = c(1, 1), seed = 123) {

  if (is.null(p_outcome)) {
    if (is.null(p_true) || is.null(p_resp)) {
      stop("give 'p_true' and 'p_resp', or 'p_outcome'", call. = FALSE)
    }
    check_p_true(p_true)
    if (!is.numeric(p_resp) || length(p_resp) != length(p_true) ||
        any(!is.finite(p_resp)) || any(p_resp < 0) || any(p_resp > 1)) {
      stop("'p_resp' must hold one response probability between 0 and 1 ",
           "per dose in 'p_true'", call. = FALSE)
    }
    p_outcome <- uboin_outcome_probs(p_true, p_resp, resp_assoc)
  } else {
    if (!is.null(p_true) || !is.null(p_resp)) {
      stop("give 'p_true' and 'p_resp', or 'p_outcome', not both",
           call. = FALSE)
    }
    p_outcome <- uboin_check_outcome(p_outcome)
    p_true <- p_outcome[, 1L] + p_outcome[, 3L]
    p_resp <- p_outcome[, 3L] + p_outcome[, 4L]
  }
  n_doses <- nrow(p_outcome)

  check_count(n_cohort, "n_cohort", 1L)
  check_count(cohort_size, "cohort_size", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(n_stage1, "n_stage1", 1L)
  check_count(start_dose, "start_dose", 1L)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of doses", call. = FALSE)
  }
  max_total_pts <- as.integer(n_cohort * cohort_size)
  if (is.null(n_stop)) n_stop <- max_total_pts
  check_count(n_stop, "n_stop", 1L)
  if (n_stop < n_stage1) {
    stop("'n_stop' must not be smaller than 'n_stage1'", call. = FALSE)
  }

  des <- uboin_check_design(tox_upper, resp_lower, target, p_saf, p_tox,
                            utility, cutoff_eli, cutoff_eli_resp, prior_tox,
                            prior_resp)
  tab <- uboin_tables(des, max_total_pts)

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- uboin_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_outcome = p_outcome,
    n_cohort = as.integer(n_cohort),
    cohort_size = as.integer(cohort_size),
    start_dose = as.integer(start_dose),
    n_stage1 = as.integer(n_stage1),
    n_stop = as.integer(n_stop),
    utility = des$utility,
    prior_a = 0.25,
    max_n = max_total_pts,
    elim = tab$elim,
    toxic = tab$toxic,
    futile = tab$futile,
    b_esc = tab$b_esc,
    b_deesc = tab$b_deesc
  )

  dose_names <- paste0("DL", seq_len(n_doses))
  n_outcome <- res$n_outcome
  dimnames(n_outcome) <- list(NULL, dose_names, uboin_outcome_names())
  n_pts <- apply(n_outcome, c(1L, 2L), sum)
  n_tox <- n_outcome[, , 1L] + n_outcome[, , 3L]
  n_resp <- n_outcome[, , 3L] + n_outcome[, , 4L]
  shape <- function(x) {
    x <- matrix(as.integer(x), nrow = n_trials, ncol = n_doses)
    colnames(x) <- dose_names
    x
  }
  eliminated <- matrix(res$eliminated, nrow = n_trials, ncol = n_doses)
  colnames(eliminated) <- dose_names

  structure(
    list(
      n_outcome = n_outcome,
      n_pts = shape(n_pts),
      n_tox = shape(n_tox),
      n_resp = shape(n_resp),
      eliminated = eliminated,
      cohorts_used = res$cohorts_used,
      stage2 = res$stage2,
      stop_reason = uboin_stop_reasons()[res$stop_code + 1L],
      boundary = list(lambda_e = des$lambda_e, lambda_d = des$lambda_d,
                      b_esc = tab$b_esc, b_deesc = tab$b_deesc,
                      b_elim = tab$b_elim),
      design = "U-BOIN",
      settings = list(
        tox_upper = des$tox_upper, resp_lower = des$resp_lower,
        target = des$target, p_true = unname(p_true),
        p_resp = unname(p_resp), p_outcome = p_outcome,
        utility = des$utility, n_cohort = n_cohort,
        cohort_size = as.integer(cohort_size), n_trials = n_trials,
        n_stage1 = n_stage1, n_stop = n_stop, start_dose = start_dose,
        p_saf = des$p_saf, p_tox = des$p_tox, cutoff_eli = cutoff_eli,
        cutoff_eli_resp = cutoff_eli_resp, prior_tox = des$prior_tox,
        prior_resp = des$prior_resp, max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "uboin_trials"
  )
}
