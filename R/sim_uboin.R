#' Operating Characteristics of the U-BOIN Design
#'
#' @description
#'   Simulate a trial run with the utility-based BOIN (U-BOIN) phase I/II design
#'   of Zhou, Lee and Yuan (2019) many times under one scenario of true DLT and
#'   response probabilities, and summarize how often each dose is selected as
#'   the optimal biological dose (OBD), how many patients are treated at each
#'   dose and how many DLTs and responses are observed.
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
#'   toxicity rule. Defaults to \code{c(1, 1)}.
#'
#' @param prior_resp
#'   Numeric vector of length two. Beta prior of the response probability in
#'   the futility rule. Defaults to \code{c(1, 1)}.
#'
#' @param keep_trials
#'   Logical scalar. Keep the trial by trial data in the result.
#'   Defaults to \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Report progress while the simulation runs.
#'   Defaults to \code{FALSE}.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. The state of the calling
#'   session is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{uboin_oc}, which is a list with components
#'   \item{sel_percent}{Percentage of trials selecting each dose as the OBD.}
#'   \item{percent_no_obd}{Percentage of trials without an OBD, which the
#'     article calls the percentage of early stopping.}
#'   \item{n_pts_dose, n_tox_dose, n_resp_dose}{Average numbers of patients,
#'     DLTs and responses at each dose.}
#'   \item{total_n_pts, total_n_tox, total_n_resp}{Their totals.}
#'   \item{true_utility}{True mean utility of each dose,
#'     \code{p_outcome \%*\% utility}.}
#'   \item{percent_stage2}{Percentage of trials that reached stage II.}
#'   \item{stop_reason_percent}{Percentage of trials by reason for ending.}
#'   \item{trials, obd}{The trials and their OBDs when \code{keep_trials} is
#'     \code{TRUE}.}
#'   together with the design parameters and the call.
#'
#' @details
#'   The trials are simulated with \code{\link{uboin_simulate}}, whose help page
#'   describes the rules, and the OBD is selected with
#'   \code{\link{uboin_select_obd}}: the admissible dose with the largest
#'   posterior mean utility at the end of the trial. A trial stopped because
#'   the lowest dose was eliminated or no dose was admissible selects none.
#'
#'   With \code{resp_assoc = 0.2}, cohorts of 3, \code{n_cohort = 18},
#'   \code{tox_upper = 0.30}, \code{resp_lower = 0.20} and the other defaults,
#'   the results agree with Table 3 of the article and Table S9 of its
#'   supporting information within simulation error
#'   (\code{inst/reproduce/reproduce_Zhou2019.R}).
#'
#' @references
#'   Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian optimal
#'   interval (U-BOIN) phase I/II design to identify the optimal biological
#'   dose for targeted and immune therapies. Statistics in Medicine, 38(28),
#'   S5299-S5316. \doi{10.1002/sim.8361}
#'
#' @examples
#' # Scenario 1 of Simulation A of the article
#' oc <- sim_uboin(
#'   tox_upper = 0.30,
#'   resp_lower = 0.20,
#'   p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
#'   p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65),
#'   n_cohort = 18,
#'   cohort_size = 3,
#'   resp_assoc = 0.2,
#'   n_trials = 500,
#'   seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_uboin_multi}}, \code{\link{uboin_simulate}},
#'   \code{\link{uboin_select_obd}}, \code{\link{uboin_utility_table}}
#'
#' @export
sim_uboin <- function(tox_upper, resp_lower, p_true = NULL, p_resp = NULL,
                      n_cohort, cohort_size, n_trials = 10000, target = NULL,
                      resp_assoc = 0, p_outcome = NULL,
                      utility = c(0, 30, 50, 100), n_stage1 = 12,
                      n_stop = NULL, start_dose = 1, p_saf = NULL,
                      p_tox = NULL, cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
                      prior_tox = c(1, 1), prior_resp = c(1, 1),
                      keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- uboin_simulate(
    tox_upper = tox_upper, resp_lower = resp_lower, p_true = p_true,
    p_resp = p_resp, n_cohort = n_cohort, cohort_size = cohort_size,
    n_trials = n_trials, target = target, resp_assoc = resp_assoc,
    p_outcome = p_outcome, utility = utility, n_stage1 = n_stage1,
    n_stop = n_stop, start_dose = start_dose, p_saf = p_saf, p_tox = p_tox,
    cutoff_eli = cutoff_eli, cutoff_eli_resp = cutoff_eli_resp,
    prior_tox = prior_tox, prior_resp = prior_resp, seed = seed
  )

  if (verbose) message("Selecting the OBD ...")

  selection <- uboin_select_obd(
    trials$n_outcome, tox_upper = tox_upper, resp_lower = resp_lower,
    eliminated = trials$eliminated, utility = utility,
    cutoff_eli = cutoff_eli, cutoff_eli_resp = cutoff_eli_resp,
    prior_tox = prior_tox, prior_resp = prior_resp
  )

  # A trial stopped because the lowest dose was eliminated or no dose was
  # admissible selects no dose; the selection from its data finds none as
  # well, and this makes the rule explicit.
  obd <- selection$obd
  obd[trials$stop_reason %in% c("lowest_dose_eliminated",
                                "no_admissible_dose")] <- NA_integer_

  if (verbose) message("Done.")

  uboin_summarize(trials, obd, keep_trials, call_expr)
}
