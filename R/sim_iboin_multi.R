#' Operating Characteristics of the iBOIN Design in Several Scenarios
#'
#' @description
#'   Run \code{\link{sim_iboin}} for several true dose-toxicity curves with
#'   the same design and the same prior information, and collect the results
#'   in one table.
#'
#' @inheritParams sim_iboin
#'
#' @param scenarios
#'   A list of scenarios, each a numeric vector of true DLT probabilities or
#'   a list with the components \code{name} and \code{p_true}, as for
#'   \code{\link{sim_boin_multi}}. Every scenario must have the same number
#'   of doses.
#'
#' @return
#'   An object of class \code{boin_oc_multi}, which is a list with components
#'   \item{results}{The result of \code{\link{sim_iboin}} for every scenario.}
#'   \item{summary_table}{Data frame with the true DLT rates, the selection
#'     percentages and the average numbers of patients and DLTs of every
#'     scenario.}
#'   \item{scenario_names}{The names of the scenarios.}
#'   \item{n_doses}{The number of doses.}
#'   \item{call}{The call.}
#'
#' @details
#'   Every scenario is simulated with the same seed, as in
#'   \code{\link{sim_boin_multi}}.
#'
#' @references
#'   Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
#'   Incorporating historical information to improve phase I clinical trials.
#'   Pharmaceutical Statistics, 20, 1017-1034. \doi{10.1002/pst.2121}
#'
#' @examples
#' # The prior MTD is dose 3; it is the true MTD in the first scenario only
#' res <- sim_iboin_multi(
#'   target = 0.30,
#'   scenarios = list(
#'     "MTD at dose 3" = c(0.08, 0.15, 0.31, 0.45, 0.55),
#'     "MTD at dose 4" = c(0.09, 0.12, 0.15, 0.30, 0.45)
#'   ),
#'   skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
#'   prior_ess = 3,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_earlystop = 100,
#'   n_trials = 300,
#'   seed = 123
#' )
#' res
#'
#' @seealso \code{\link{sim_iboin}}, \code{\link{sim_boin_multi}}
#'
#' @export
sim_iboin_multi <- function(target, scenarios, skeleton, prior_ess, n_cohort,
                            cohort_size, n_trials = 10000, start_dose = 1,
                            n_earlystop = 18, p_saf = NULL, p_tox = NULL,
                            cutoff_eli = 0.95, extrasafe = FALSE,
                            offset = 0.05, robust = FALSE,
                            select_with_prior = TRUE, min_mtd_sample = 1,
                            overdose_cutoff = NULL,
                            n_earlystop_rule = c("with_stay", "simple"),
                            keep_trials = FALSE, verbose = FALSE,
                            seed = 123) {

  call_expr <- match.call()
  n_earlystop_rule <- match.arg(n_earlystop_rule)

  scenarios <- normalize_scenarios(scenarios)
  scenario_names <- vapply(scenarios, function(z) z$name, character(1))
  n_scenarios <- length(scenarios)
  n_doses <- length(scenarios[[1L]]$p_true)

  results <- vector("list", n_scenarios)
  names(results) <- scenario_names

  for (i in seq_len(n_scenarios)) {
    if (verbose) {
      message("Scenario ", i, " of ", n_scenarios, ": ", scenario_names[i])
    }
    results[[i]] <- sim_iboin(
      target = target, p_true = scenarios[[i]]$p_true, skeleton = skeleton,
      prior_ess = prior_ess, n_cohort = n_cohort, cohort_size = cohort_size,
      n_trials = n_trials, start_dose = start_dose,
      n_earlystop = n_earlystop, p_saf = p_saf, p_tox = p_tox,
      cutoff_eli = cutoff_eli, extrasafe = extrasafe, offset = offset,
      robust = robust, select_with_prior = select_with_prior,
      min_mtd_sample = min_mtd_sample, overdose_cutoff = overdose_cutoff,
      n_earlystop_rule = n_earlystop_rule, keep_trials = keep_trials,
      verbose = FALSE, seed = seed
    )
  }

  summary_table <- oc_multi_table(results, scenario_names, n_doses)

  structure(
    list(
      results = results,
      summary_table = summary_table,
      scenario_names = scenario_names,
      n_doses = n_doses,
      call = call_expr
    ),
    class = "boin_oc_multi"
  )
}
