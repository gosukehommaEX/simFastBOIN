#' Operating Characteristics of the U-BOIN Design Across Scenarios
#'
#' @description
#'   Run \code{\link{sim_uboin}} under several scenarios of true DLT and
#'   response probabilities with the same design, and collect the main
#'   operating characteristics of every scenario in one table.
#'
#' @param tox_upper
#'   Numeric scalar. Highest acceptable DLT probability, for example 0.30.
#'
#' @param resp_lower
#'   Numeric scalar. Lowest acceptable response probability, for example 0.20.
#'
#' @param scenarios
#'   List of scenarios, each a list with elements \code{p_true} and
#'   \code{p_resp} and optionally \code{name}. Names of the list are used as
#'   scenario names. All scenarios must have the same number of doses.
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
#' @param keep_trials
#'   Logical scalar. Keep the trial by trial data in the result of every
#'   scenario. Defaults to \code{FALSE}.
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
#'   An object of class \code{uboin_oc_multi}, which is a list with components
#'   \item{results}{Named list of \code{uboin_oc} objects, one per scenario.}
#'   \item{summary_table}{Data frame with the true rates, the selection
#'     percentages and the average numbers of patients, DLTs and responses of
#'     every scenario and dose.}
#'   \item{scenario_names}{Character vector of scenario names.}
#'   \item{n_doses}{Number of doses.}
#'   \item{call}{The matched call.}
#'
#' @details
#'   Every scenario is simulated with the same seed, so that the scenarios
#'   differ by their probabilities and not by the random numbers. The outcome
#'   probabilities of a scenario are given by its marginal rates and
#'   \code{resp_assoc}; to give them directly, run \code{\link{sim_uboin}} with
#'   \code{p_outcome} for each scenario.
#'
#' @references
#'   Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian optimal
#'   interval (U-BOIN) phase I/II design to identify the optimal biological
#'   dose for targeted and immune therapies. Statistics in Medicine, 38(28),
#'   S5299-S5316. \doi{10.1002/sim.8361}
#'
#' @examples
#' scenarios <- list(
#'   "Scenario 1" = list(p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
#'                       p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65)),
#'   "Scenario 3" = list(p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
#'                       p_resp = c(0.08, 0.46, 0.25, 0.20, 0.10))
#' )
#'
#' oc <- sim_uboin_multi(
#'   tox_upper = 0.30,
#'   resp_lower = 0.20,
#'   scenarios = scenarios,
#'   n_cohort = 18,
#'   cohort_size = 3,
#'   resp_assoc = 0.2,
#'   n_trials = 200,
#'   seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_uboin}}
#'
#' @export
sim_uboin_multi <- function(tox_upper, resp_lower, scenarios, n_cohort,
                            cohort_size, n_trials = 10000, target = NULL,
                            resp_assoc = 0, utility = c(0, 30, 50, 100),
                            n_stage1 = 12, n_stop = NULL, start_dose = 1,
                            p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                            cutoff_eli_resp = 0.90, prior_tox = c(1, 1),
                            prior_resp = c(1, 1), keep_trials = FALSE,
                            verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  check_flag(verbose, "verbose")

  scenarios <- normalize_backfill_scenarios(scenarios)
  scenario_names <- vapply(scenarios, function(z) z$name, character(1))
  n_doses <- length(scenarios[[1L]]$p_true)

  results <- vector("list", length(scenarios))
  names(results) <- scenario_names

  for (i in seq_along(scenarios)) {
    if (verbose) {
      message("Scenario ", i, " of ", length(scenarios), ": ",
              scenario_names[i])
    }
    results[[i]] <- sim_uboin(
      tox_upper = tox_upper, resp_lower = resp_lower,
      p_true = scenarios[[i]]$p_true, p_resp = scenarios[[i]]$p_resp,
      n_cohort = n_cohort, cohort_size = cohort_size, n_trials = n_trials,
      target = target, resp_assoc = resp_assoc, utility = utility,
      n_stage1 = n_stage1, n_stop = n_stop, start_dose = start_dose,
      p_saf = p_saf, p_tox = p_tox, cutoff_eli = cutoff_eli,
      cutoff_eli_resp = cutoff_eli_resp, prior_tox = prior_tox,
      prior_resp = prior_resp, keep_trials = keep_trials, verbose = FALSE,
      seed = seed
    )
  }

  structure(
    list(
      results = results,
      summary_table = boin_et_multi_table(results, scenario_names, n_doses),
      scenario_names = scenario_names,
      n_doses = n_doses,
      call = call_expr
    ),
    class = "uboin_oc_multi"
  )
}
