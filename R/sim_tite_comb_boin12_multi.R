#' Operating Characteristics of the TITE-Comb-BOIN12 Design Across Scenarios
#'
#' @description
#'   Run \code{\link{sim_tite_comb_boin12}} for several scenarios of DLT and
#'   response probabilities with the same design and collect the results in
#'   one table.
#'
#' @inheritParams sim_tite_comb_boin12
#' @inheritParams sim_comb_boin12_multi
#'
#' @return
#'   An object of class \code{comb_boin12_oc_multi}, as returned by
#'   \code{\link{sim_comb_boin12_multi}}, whose summary table also has the
#'   average trial duration of every scenario (\code{Duration}).
#'
#' @details
#'   Every scenario is simulated with the same \code{seed}, so that a scenario
#'   simulated here matches the same scenario simulated on its own with
#'   \code{\link{sim_tite_comb_boin12}}.
#'
#' @examples
#' # Scenarios 16 and 17 of Lu et al. (2025), windows of 1 and 2 months
#' scenarios <- list(
#'   "Scenario 16" = list(
#'     p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2),
#'     p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)),
#'   "Scenario 17" = list(
#'     p_true = matrix(c(0.05, 0.10, 0.15, 0.20, 0.40, 0.50), nrow = 2),
#'     p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.55), nrow = 2))
#' )
#' oc <- sim_tite_comb_boin12_multi(
#'   target = 0.33, target_resp = 0.2, scenarios = scenarios,
#'   n_cohort = 14, cohort_size = 3, window = 1, resp_window = 2,
#'   accrual_rate = 2, event_time = "uniform", n_trials = 100,
#'   psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11, benchmark_weight = 0.6,
#'   extrasafe = TRUE, offset = 0.10, est_method = "observed", seed = 1
#' )
#' oc
#'
#' @seealso \code{\link{sim_tite_comb_boin12}}
#'
#' @export
sim_tite_comb_boin12_multi <- function(target, target_resp, scenarios,
                                       n_cohort, cohort_size, window,
                                       accrual_rate, resp_window = window,
                                       accrual = c("uniform", "exponential",
                                                   "fixed"),
                                       event_time = c("weibull", "uniform"),
                                       late_fraction = 0.5,
                                       resp_late_fraction = 0.5,
                                       min_completed_share = 0.5,
                                       pending_estimate = c("conditional",
                                                            "marginal"),
                                       pending_use = c("all", "utility"),
                                       n_trials = 10000,
                                       start_dose = c(1, 1),
                                       n_earlystop = NULL,
                                       n_earlystop_rule = c("with_stay",
                                                            "simple"),
                                       p_saf = NULL, p_tox = NULL,
                                       cutoff_eli = 0.95,
                                       cutoff_eli_resp = 0.90,
                                       extrasafe = FALSE, offset = 0.05,
                                       psi00 = 40, psi11 = 60,
                                       benchmark_weight = 0.5,
                                       utility_prior = c(1, 1), n_star = 6,
                                       tie_break = c("higher", "random"),
                                       no_candidate = c("nearest", "stop"),
                                       resp_cor = 0,
                                       resp_cor_type = c("copula", "binary"),
                                       est_method = c("logistic",
                                                      "observed"),
                                       obd_bound = c("rectangle", "contour"),
                                       dose_a = NULL, dose_b = NULL,
                                       logistic_scale = c(2.5, 2.5),
                                       n_draws = 4000, keep_trials = FALSE,
                                       verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  pending_estimate <- match.arg(pending_estimate)
  pending_use <- match.arg(pending_use)
  n_earlystop_rule <- match.arg(n_earlystop_rule)
  tie_break <- match.arg(tie_break)
  no_candidate <- match.arg(no_candidate)
  resp_cor_type <- match.arg(resp_cor_type)
  est_method <- match.arg(est_method)
  obd_bound <- match.arg(obd_bound)
  check_flag(verbose, "verbose")

  scenarios <- normalize_comb_resp_scenarios(scenarios)
  scenario_names <- vapply(scenarios, function(z) z$name, character(1))

  results <- vector("list", length(scenarios))
  names(results) <- scenario_names

  for (i in seq_along(scenarios)) {
    if (verbose) {
      message("Scenario ", i, " of ", length(scenarios), ": ",
              scenario_names[i])
    }
    results[[i]] <- sim_tite_comb_boin12(
      target = target, target_resp = target_resp,
      p_true = scenarios[[i]]$p_true, p_resp = scenarios[[i]]$p_resp,
      n_cohort = n_cohort, cohort_size = cohort_size, window = window,
      accrual_rate = accrual_rate, resp_window = resp_window,
      accrual = accrual, event_time = event_time,
      late_fraction = late_fraction, resp_late_fraction = resp_late_fraction,
      min_completed_share = min_completed_share,
      pending_estimate = pending_estimate, pending_use = pending_use,
      n_trials = n_trials, start_dose = start_dose,
      n_earlystop = n_earlystop,
      n_earlystop_rule = n_earlystop_rule, p_saf = p_saf, p_tox = p_tox,
      cutoff_eli = cutoff_eli, cutoff_eli_resp = cutoff_eli_resp,
      extrasafe = extrasafe, offset = offset, psi00 = psi00, psi11 = psi11,
      benchmark_weight = benchmark_weight, utility_prior = utility_prior,
      n_star = n_star, tie_break = tie_break, no_candidate = no_candidate,
      resp_cor = resp_cor, resp_cor_type = resp_cor_type,
      est_method = est_method, obd_bound = obd_bound, dose_a = dose_a,
      dose_b = dose_b, logistic_scale = logistic_scale, n_draws = n_draws,
      keep_trials = keep_trials, verbose = FALSE, seed = seed
    )
  }

  structure(
    list(
      results = results,
      summary_table = comb_boin12_multi_table(results, scenario_names),
      scenario_names = scenario_names,
      n_rows = nrow(scenarios[[1L]]$p_true),
      n_cols = ncol(scenarios[[1L]]$p_true),
      call = call_expr
    ),
    class = "comb_boin12_oc_multi"
  )
}
