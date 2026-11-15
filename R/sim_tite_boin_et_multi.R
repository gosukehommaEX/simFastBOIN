#' Operating Characteristics of the TITE-BOIN-ET Design Across Scenarios
#'
#' @description
#'   Run \code{\link{sim_tite_boin_et}} for several scenarios of DLT and
#'   response probabilities with the same design and collect the results in
#'   one table.
#'
#' @inheritParams sim_tite_boin_et
#' @inheritParams sim_boin_et_multi
#'
#' @return
#'   An object of class \code{boin_et_oc_multi}, as returned by
#'   \code{\link{sim_boin_et_multi}}, whose \code{results} are the
#'   \code{boin_et_oc} objects of \code{\link{sim_tite_boin_et}}.
#'
#' @details
#'   Every scenario is simulated with the same \code{seed}.
#'
#' @examples
#' scenarios <- list(
#'   "Scenario 1" = list(p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
#'                       p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70)),
#'   "Scenario 7" = list(p_true = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
#'                       p_resp = c(0.30, 0.70, 0.60, 0.50, 0.40, 0.30))
#' )
#' res <- sim_tite_boin_et_multi(
#'   target = 0.30, target_resp = 0.60, scenarios = scenarios,
#'   n_cohort = 12, cohort_size = 3, window = 30, resp_window = 45,
#'   accrual_rate = 1 / 10, n_trials = 300, seed = 123
#' )
#' res
#'
#' @seealso \code{\link{sim_tite_boin_et}}
#' @export
sim_tite_boin_et_multi <- function(target, target_resp, scenarios, n_cohort,
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
                                   cutoff_eli = 0.95, cutoff_eli_resp = 0.99,
                                   p_eli = NULL, p_eli_resp = NULL,
                                   boinet_compatible = FALSE,
                                   prior = rep(1 / 6, 6),
                                   est_method = c("observed", "multi_iso",
                                                  "fp_logistic"),
                                   obd_method = c("max_resp",
                                                  "utility_weighted",
                                                  "utility_truncated_linear",
                                                  "utility_scoring"),
                                   w1 = 0.33, w2 = 1.09, t_low = NULL,
                                   t_upp = NULL, e_low = NULL, e_upp = NULL,
                                   psi00 = 40, psi11 = 60,
                                   keep_trials = FALSE, verbose = FALSE,
                                   seed = 123) {

  call_expr <- match.call()
  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  est_method <- match.arg(est_method)
  obd_method <- match.arg(obd_method)
  check_flag(verbose, "verbose")

  scenarios <- normalize_backfill_scenarios(scenarios)
  scenario_names <- vapply(scenarios, function(z) z$name, character(1))
  n_doses <- length(scenarios[[1L]]$p_true)

  results <- vector("list", length(scenarios))
  names(results) <- scenario_names

  for (i in seq_along(scenarios)) {
    if (verbose) {
      message("Scenario ", i, " of ", length(scenarios), ": ", scenario_names[i])
    }
    results[[i]] <- sim_tite_boin_et(
      target = target, target_resp = target_resp,
      p_true = scenarios[[i]]$p_true, p_resp = scenarios[[i]]$p_resp,
      n_cohort = n_cohort, cohort_size = cohort_size, window = window,
      accrual_rate = accrual_rate, resp_window = resp_window,
      accrual = accrual, event_time = event_time,
      late_fraction = late_fraction, resp_late_fraction = resp_late_fraction,
      resp_cor = resp_cor, min_completed_share = min_completed_share,
      n_trials = n_trials, start_dose = start_dose,
      n_earlystop = n_earlystop, p_saf = p_saf, p_tox = p_tox,
      p_saf_resp = p_saf_resp, cutoff_eli = cutoff_eli,
      cutoff_eli_resp = cutoff_eli_resp, p_eli = p_eli,
      p_eli_resp = p_eli_resp, boinet_compatible = boinet_compatible,
      prior = prior, est_method = est_method, obd_method = obd_method,
      w1 = w1, w2 = w2, t_low = t_low, t_upp = t_upp, e_low = e_low,
      e_upp = e_upp, psi00 = psi00, psi11 = psi11,
      keep_trials = keep_trials, verbose = FALSE, seed = seed
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
    class = "boin_et_oc_multi"
  )
}
