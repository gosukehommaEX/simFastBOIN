#' Operating Characteristics of a BE-BOIN Design Across Several Scenarios
#'
#' @description
#'   Run \code{\link{sim_be_boin}} under a list of scenarios of DLT and response
#'   probabilities that share the same design, and collect the results into one
#'   table.
#'
#' @inheritParams sim_be_boin
#' @inheritParams sim_bf_boin_multi
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed, applied to every scenario, so
#'   that a scenario simulated here matches the same scenario simulated on its
#'   own with \code{\link{sim_be_boin}}. Defaults to 123.
#'
#' @return
#'   An object of class \code{c("be_boin_oc_multi", "backfill_oc_multi",
#'   "boin_oc_multi")}, with the components described in
#'   \code{\link{sim_bf_boin_multi}}; \code{results} holds \code{be_boin_oc}
#'   objects.
#'
#' @examples
#' scenarios <- list(
#'   list(name = "MTD at dose 3", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
#'        p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58)),
#'   list(name = "Plateau", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
#'        p_resp = c(0.30, 0.32, 0.35, 0.36, 0.36))
#' )
#'
#' oc <- sim_be_boin_multi(
#'   target = 0.25,
#'   scenarios = scenarios,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   window = 3,
#'   accrual_rate = 2,
#'   n_earlystop = 9,
#'   n_trials = 200,
#'   seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_be_boin}}, \code{\link{sim_bf_boin_multi}}
#'
#' @export
sim_be_boin_multi <- function(target, scenarios, n_cohort, cohort_size,
                              window, accrual_rate, n_cap = 12,
                              backfill_dose = c("highest", "lowest"),
                              conflict_dose = c("lowest", "highest"),
                              no_slot = c("wait", "leave"),
                              accrual = c("exponential", "uniform", "fixed"),
                              dlt_time = c("weibull", "uniform"),
                              late_fraction = 0.5, resp_window = window,
                              resp_late_fraction = 0.5, resp_cor = 0,
                              max_pending_ratio = 0.49, min_completed = 0,
                              min_follow_up = 0.25,
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
  backfill_dose <- match.arg(backfill_dose)
  conflict_dose <- match.arg(conflict_dose)
  no_slot <- match.arg(no_slot)
  accrual <- match.arg(accrual)
  dlt_time <- match.arg(dlt_time)
  n_earlystop_rule <- match.arg(n_earlystop_rule)

  scenarios <- normalize_backfill_scenarios(scenarios)
  scenario_names <- vapply(scenarios, function(z) z$name, character(1))
  n_doses <- length(scenarios[[1L]]$p_true)

  results <- vector("list", length(scenarios))
  names(results) <- scenario_names

  for (i in seq_along(scenarios)) {
    if (verbose) {
      message("Scenario ", i, " of ", length(scenarios), ": ", scenario_names[i])
    }
    results[[i]] <- sim_be_boin(
      target = target, p_true = scenarios[[i]]$p_true,
      p_resp = scenarios[[i]]$p_resp, n_cohort = n_cohort,
      cohort_size = cohort_size, window = window, accrual_rate = accrual_rate,
      n_cap = n_cap, backfill_dose = backfill_dose,
      conflict_dose = conflict_dose, no_slot = no_slot, accrual = accrual,
      dlt_time = dlt_time, late_fraction = late_fraction,
      resp_window = resp_window, resp_late_fraction = resp_late_fraction,
      resp_cor = resp_cor, max_pending_ratio = max_pending_ratio,
      min_completed = min_completed, min_follow_up = min_follow_up,
      n_trials = n_trials, start_dose = start_dose,
      n_earlystop = n_earlystop, p_saf = p_saf, p_tox = p_tox,
      cutoff_eli = cutoff_eli, extrasafe = extrasafe, offset = offset,
      bound_mtd = bound_mtd, mtd_max_estimate = mtd_max_estimate,
      min_mtd_sample = min_mtd_sample, overdose_cutoff = overdose_cutoff,
      n_earlystop_rule = n_earlystop_rule, keep_trials = keep_trials,
      verbose = FALSE, seed = seed
    )
  }

  structure(
    list(
      results = results,
      summary_table = backfill_oc_multi_table(results, scenario_names, n_doses),
      scenario_names = scenario_names,
      n_doses = n_doses,
      call = call_expr
    ),
    class = c("be_boin_oc_multi", "backfill_oc_multi", "boin_oc_multi")
  )
}
