#' Operating Characteristics of the TITE-BOIN12 Design Across Scenarios
#'
#' @description
#'   Run \code{\link{sim_tite_boin12}} for several scenarios of DLT and
#'   response probabilities with the same design and collect the results in
#'   one table.
#'
#' @inheritParams sim_tite_boin12
#' @inheritParams sim_boin12_multi
#'
#' @return
#'   An object of class \code{boin12_oc_multi}, as returned by
#'   \code{\link{sim_boin12_multi}}. The average trial duration of each
#'   scenario is the component \code{duration} of its element of
#'   \code{results}.
#'
#' @details
#'   Every scenario is simulated with the same \code{seed}, so that a scenario
#'   simulated here matches the same scenario simulated on its own with
#'   \code{\link{sim_tite_boin12}}.
#'
#' @examples
#' # Scenarios 1 and 2 of Zhou et al. (2022)
#' scenarios <- list(
#'   "Scenario 1" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                       p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
#'   "Scenario 2" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                       p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
#' )
#' oc <- sim_tite_boin12_multi(
#'   target = 0.35, target_resp = 0.25, scenarios = scenarios,
#'   n_cohort = 12, cohort_size = 3, window = 45, resp_window = 60,
#'   accrual_rate = 3 / 30, n_trials = 200, seed = 1
#' )
#' oc
#' sapply(oc$results, function(r) r$duration)
#'
#' @seealso \code{\link{sim_tite_boin12}}, \code{\link{sim_boin12_multi}}
#'
#' @export
sim_tite_boin12_multi <- function(target, target_resp, scenarios, n_cohort,
                                  cohort_size, window, accrual_rate,
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
                                  est_method = c("posterior_mean", "observed",
                                                 "tail_probability"),
                                  keep_trials = FALSE, verbose = FALSE,
                                  seed = 123) {

  call_expr <- match.call()
  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  pending_estimate <- match.arg(pending_estimate)
  n_earlystop_rule <- match.arg(n_earlystop_rule)
  tie_break <- match.arg(tie_break)
  no_candidate <- match.arg(no_candidate)
  resp_cor_type <- match.arg(resp_cor_type)
  est_method <- match.arg(est_method)
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
    results[[i]] <- sim_tite_boin12(
      target = target, target_resp = target_resp,
      p_true = scenarios[[i]]$p_true, p_resp = scenarios[[i]]$p_resp,
      n_cohort = n_cohort, cohort_size = cohort_size, window = window,
      accrual_rate = accrual_rate, resp_window = resp_window,
      accrual = accrual, event_time = event_time,
      late_fraction = late_fraction, resp_late_fraction = resp_late_fraction,
      min_completed_share = min_completed_share,
      pending_estimate = pending_estimate, n_trials = n_trials,
      start_dose = start_dose, n_earlystop = n_earlystop,
      n_earlystop_rule = n_earlystop_rule, p_saf = p_saf, p_tox = p_tox,
      cutoff_eli = cutoff_eli, cutoff_eli_resp = cutoff_eli_resp,
      psi00 = psi00, psi11 = psi11, benchmark_weight = benchmark_weight,
      utility_prior = utility_prior, n_star = n_star, n_explore = n_explore,
      tie_break = tie_break, no_candidate = no_candidate,
      resp_cor = resp_cor, resp_cor_type = resp_cor_type,
      est_method = est_method, keep_trials = keep_trials, verbose = FALSE,
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
    class = "boin12_oc_multi"
  )
}
