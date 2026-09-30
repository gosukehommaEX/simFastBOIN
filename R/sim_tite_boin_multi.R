#' Operating Characteristics of a TITE-BOIN Design Across Several Scenarios
#'
#' @description
#'   Run \code{\link{sim_tite_boin}} under a list of dose-toxicity scenarios that
#'   share the same design, and collect the results into one table.
#'
#' @inheritParams sim_tite_boin
#'
#' @param scenarios
#'   A list of scenarios. Each element is a list with a \code{name} and a
#'   \code{p_true} vector, or a plain numeric vector of true DLT probabilities in
#'   which case the list names are used as scenario names. Every scenario must
#'   have the same number of doses.
#'
#' @param keep_trials
#'   Logical scalar. Keep the trial by trial data of every scenario.
#'   Defaults to \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Report progress while the simulations run.
#'   Defaults to \code{FALSE}.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed, applied to every scenario, so
#'   that a scenario simulated here matches the same scenario simulated on its own
#'   with \code{\link{sim_tite_boin}}. Defaults to 123.
#'
#' @return
#'   An object of class \code{c("tite_boin_oc_multi", "boin_oc_multi")}, which
#'   is a list with components
#'   \item{results}{Named list of \code{tite_boin_oc} objects, one per scenario.}
#'   \item{summary_table}{Data frame collecting all scenarios, with the four rows
#'     of \code{\link{sim_boin_multi}} and two more per scenario for the average
#'     trial duration and the percentage of trials with accrual suspended, both
#'     in the last column.}
#'   \item{scenario_names}{Character vector of scenario names.}
#'   \item{n_doses}{Number of dose levels.}
#'   \item{call}{The matched call.}
#'
#' @examples
#' scenarios <- list(
#'   list(name = "MTD at dose 3", p_true = c(0.05, 0.15, 0.30, 0.45, 0.60)),
#'   list(name = "All toxic",     p_true = c(0.35, 0.45, 0.55, 0.65, 0.75))
#' )
#'
#' oc <- sim_tite_boin_multi(
#'   target = 0.30,
#'   scenarios = scenarios,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   window = 3,
#'   accrual_rate = 2,
#'   n_trials = 200,
#'   seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_tite_boin}}, \code{\link{sim_boin_multi}}
#'
#' @export
sim_tite_boin_multi <- function(target, scenarios, n_cohort, cohort_size,
                                window, accrual_rate,
                                method = c("imputation", "ess"),
                                accrual = c("exponential", "uniform", "fixed"),
                                dlt_time = c("weibull", "uniform"),
                                late_fraction = 0.5, prior_weights = c(1, 1, 1) / 3,
                                max_pending_ratio = NULL, min_completed = NULL,
                                n_trials = 10000, start_dose = 1,
                                n_earlystop = 18, p_saf = NULL, p_tox = NULL,
                                cutoff_eli = 0.95, extrasafe = FALSE,
                                offset = 0.05, bound_mtd = FALSE,
                                mtd_max_estimate = NULL, min_mtd_sample = 1,
                                overdose_cutoff = NULL,
                                n_earlystop_rule = c("with_stay", "simple"),
                                keep_trials = FALSE, verbose = FALSE,
                                seed = 123) {

  call_expr <- match.call()
  method <- match.arg(method)
  accrual <- match.arg(accrual)
  dlt_time <- match.arg(dlt_time)
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
    results[[i]] <- sim_tite_boin(
      target = target, p_true = scenarios[[i]]$p_true, n_cohort = n_cohort,
      cohort_size = cohort_size, window = window, accrual_rate = accrual_rate,
      method = method, accrual = accrual, dlt_time = dlt_time,
      late_fraction = late_fraction, prior_weights = prior_weights,
      max_pending_ratio = max_pending_ratio,
      min_completed = min_completed, n_trials = n_trials,
      start_dose = start_dose, n_earlystop = n_earlystop, p_saf = p_saf,
      p_tox = p_tox, cutoff_eli = cutoff_eli, extrasafe = extrasafe,
      offset = offset, bound_mtd = bound_mtd,
      mtd_max_estimate = mtd_max_estimate, min_mtd_sample = min_mtd_sample,
      overdose_cutoff = overdose_cutoff, n_earlystop_rule = n_earlystop_rule,
      keep_trials = keep_trials, verbose = FALSE, seed = seed
    )
  }

  summary_table <- tite_oc_multi_table(results, scenario_names, n_doses)

  structure(
    list(
      results = results,
      summary_table = summary_table,
      scenario_names = scenario_names,
      n_doses = n_doses,
      call = call_expr
    ),
    class = c("tite_boin_oc_multi", "boin_oc_multi")
  )
}
