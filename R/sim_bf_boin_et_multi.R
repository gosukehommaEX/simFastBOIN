#' Operating Characteristics of the BF-BOIN-ET Design Across Scenarios
#'
#' @description
#'   Run \code{\link{sim_bf_boin_et}} for several scenarios of DLT and
#'   response probabilities with the same design and collect the results in
#'   one table.
#'
#' @inheritParams sim_bf_boin_et
#' @inheritParams sim_boin_et_multi
#'
#' @return
#'   An object of class \code{bf_boin_et_oc_multi}, which is a list with
#'   components
#'   \item{results}{Named list of \code{bf_boin_et_oc} objects, one per
#'     scenario.}
#'   \item{summary_table}{Data frame with seven rows per scenario, as shown
#'     by \code{\link{print.bf_boin_et_oc_multi}}.}
#'   \item{scenario_names}{Character vector of scenario names.}
#'   \item{n_doses}{Number of doses.}
#'   \item{call}{The matched call.}
#'
#' @details
#'   Every scenario is simulated with the same \code{seed}.
#'
#' @examples
#' scenarios <- list(
#'   "Example" = list(p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
#'                    p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70)),
#'   "Plateau" = list(p_true = c(0.05, 0.10, 0.15, 0.30, 0.45),
#'                    p_resp = c(0.10, 0.30, 0.45, 0.50, 0.52))
#' )
#' res <- sim_bf_boin_et_multi(
#'   target = 0.30, target_resp = 0.25, scenarios = scenarios,
#'   n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
#'   n_trials = 200, seed = 123
#' )
#' res
#'
#' @seealso \code{\link{sim_bf_boin_et}}
#' @export
sim_bf_boin_et_multi <- function(target, target_resp, scenarios, n_cohort,
                                 cohort_size, window, accrual_rate,
                                 resp_window = window,
                                 accrual = c("uniform", "exponential",
                                             "fixed"),
                                 event_time = c("weibull", "uniform"),
                                 late_fraction = 0.5,
                                 resp_late_fraction = 0.5, resp_cor = 0.2,
                                 n_trials = 10000, start_dose = 1,
                                 n_earlystop = 12, n_elimination = 6,
                                 suspend = 0, p_saf = NULL, p_tox = NULL,
                                 p_saf_resp = NULL, cutoff_eli = 0.95,
                                 cutoff_eli_resp = 0.90, psi00 = 40,
                                 psi01 = 0, psi10 = 100, psi11 = 60,
                                 backfill_extension = 100,
                                 keep_trials = FALSE, verbose = FALSE,
                                 seed = 123) {

  call_expr <- match.call()
  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
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
    results[[i]] <- sim_bf_boin_et(
      target = target, target_resp = target_resp,
      p_true = scenarios[[i]]$p_true, p_resp = scenarios[[i]]$p_resp,
      n_cohort = n_cohort, cohort_size = cohort_size, window = window,
      accrual_rate = accrual_rate, resp_window = resp_window,
      accrual = accrual, event_time = event_time,
      late_fraction = late_fraction, resp_late_fraction = resp_late_fraction,
      resp_cor = resp_cor, n_trials = n_trials, start_dose = start_dose,
      n_earlystop = n_earlystop, n_elimination = n_elimination,
      suspend = suspend, p_saf = p_saf, p_tox = p_tox,
      p_saf_resp = p_saf_resp, cutoff_eli = cutoff_eli,
      cutoff_eli_resp = cutoff_eli_resp, psi00 = psi00, psi01 = psi01,
      psi10 = psi10, psi11 = psi11,
      backfill_extension = backfill_extension, keep_trials = keep_trials,
      verbose = FALSE, seed = seed
    )
  }

  structure(
    list(
      results = results,
      summary_table = bf_boin_et_multi_table(results, scenario_names),
      scenario_names = scenario_names,
      n_doses = n_doses,
      call = call_expr
    ),
    class = "bf_boin_et_oc_multi"
  )
}
