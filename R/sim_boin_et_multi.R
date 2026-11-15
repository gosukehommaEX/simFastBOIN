#' Operating Characteristics of the BOIN-ET Design Across Scenarios
#'
#' @description
#'   Run \code{\link{sim_boin_et}} for several scenarios of DLT and response
#'   probabilities with the same design and collect the results in one
#'   table.
#'
#' @inheritParams sim_boin_et
#'
#' @param scenarios
#'   List of scenarios. Each element is a list with the components
#'   \code{p_true} and \code{p_resp}, and optionally \code{name}. Names of the
#'   list are used when \code{name} is absent. All scenarios must have the same
#'   number of doses.
#'
#' @return
#'   An object of class \code{boin_et_oc_multi}, which is a list with
#'   components
#'   \item{results}{Named list of \code{boin_et_oc} objects, one per scenario.}
#'   \item{summary_table}{Data frame with six rows per scenario, as shown by
#'     \code{\link{print.boin_et_oc_multi}}.}
#'   \item{scenario_names}{Character vector of scenario names.}
#'   \item{n_doses}{Number of doses.}
#'   \item{call}{The matched call.}
#'
#' @details
#'   Every scenario is simulated with the same \code{seed}.
#'
#' @examples
#' scenarios <- list(
#'   list(name = "Monotone",
#'        p_true = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.40),
#'        p_resp = c(0.05, 0.10, 0.25, 0.40, 0.65, 0.70)),
#'   list(name = "Peaked",
#'        p_true = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
#'        p_resp = c(0.30, 0.70, 0.60, 0.50, 0.40, 0.30))
#' )
#' sim_boin_et_multi(
#'   target = 0.30, target_resp = 0.60, scenarios = scenarios,
#'   n_cohort = 12, cohort_size = 3, n_trials = 500, seed = 123
#' )
#'
#' @seealso \code{\link{sim_boin_et}}
#'
#' @export
sim_boin_et_multi <- function(target, target_resp, scenarios, n_cohort,
                              cohort_size, n_trials = 10000, start_dose = 1,
                              n_earlystop = NULL, p_saf = NULL, p_tox = NULL,
                              p_saf_resp = NULL, cutoff_eli = 0.95,
                              cutoff_eli_resp = 0.99, p_eli = NULL,
                              p_eli_resp = NULL, resp_cor = 0.2,
                              resp_cor_type = c("copula", "binary"),
                              boinet_compatible = FALSE,
                              prior = rep(1 / 6, 6),
                              est_method = c("observed", "multi_iso",
                                             "fp_logistic"),
                              obd_method = c("max_resp", "utility_weighted",
                                             "utility_truncated_linear",
                                             "utility_scoring"),
                              w1 = 0.33, w2 = 1.09, t_low = NULL,
                              t_upp = NULL, e_low = NULL, e_upp = NULL,
                              psi00 = 40, psi11 = 60, keep_trials = FALSE,
                              verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  resp_cor_type <- match.arg(resp_cor_type)
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
    results[[i]] <- sim_boin_et(
      target = target, target_resp = target_resp,
      p_true = scenarios[[i]]$p_true, p_resp = scenarios[[i]]$p_resp,
      n_cohort = n_cohort, cohort_size = cohort_size, n_trials = n_trials,
      start_dose = start_dose, n_earlystop = n_earlystop, p_saf = p_saf,
      p_tox = p_tox, p_saf_resp = p_saf_resp, cutoff_eli = cutoff_eli,
      cutoff_eli_resp = cutoff_eli_resp, p_eli = p_eli,
      p_eli_resp = p_eli_resp, resp_cor = resp_cor,
      resp_cor_type = resp_cor_type, boinet_compatible = boinet_compatible,
      prior = prior, est_method = est_method,
      obd_method = obd_method, w1 = w1, w2 = w2, t_low = t_low, t_upp = t_upp,
      e_low = e_low, e_upp = e_upp, psi00 = psi00, psi11 = psi11,
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
