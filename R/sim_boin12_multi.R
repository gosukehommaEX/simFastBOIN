#' Operating Characteristics of the BOIN12 Design Across Scenarios
#'
#' @description
#'   Run \code{\link{sim_boin12}} for several scenarios of DLT and response
#'   probabilities with the same design and collect the results in one
#'   table.
#'
#' @inheritParams sim_boin12
#'
#' @param scenarios
#'   List of scenarios. Each element is a list with the components
#'   \code{p_true} and \code{p_resp}, and optionally \code{name}. Names of the
#'   list are used when \code{name} is absent. All scenarios must have the same
#'   number of doses.
#'
#' @return
#'   An object of class \code{boin12_oc_multi}, which is a list with
#'   components
#'   \item{results}{Named list of \code{boin12_oc} objects, one per scenario.}
#'   \item{summary_table}{Data frame with six rows per scenario, as shown by
#'     \code{\link{print.boin12_oc_multi}}.}
#'   \item{scenario_names}{Character vector of scenario names.}
#'   \item{n_doses}{Number of doses.}
#'   \item{call}{The matched call.}
#'
#' @details
#'   Every scenario is simulated with the same \code{seed}, so that a scenario
#'   simulated here matches the same scenario simulated on its own with
#'   \code{\link{sim_boin12}}.
#'
#' @examples
#' scenarios <- list(
#'   "Plateau" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                    p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
#'   "Increasing" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                       p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
#' )
#' oc <- sim_boin12_multi(
#'   target = 0.35, target_resp = 0.25, scenarios = scenarios,
#'   n_cohort = 12, cohort_size = 3, n_trials = 200, seed = 1
#' )
#' oc
#'
#' @seealso \code{\link{sim_boin12}}
#'
#' @export
sim_boin12_multi <- function(target, target_resp, scenarios, n_cohort,
                             cohort_size, n_trials = 10000, start_dose = 1,
                             n_earlystop = NULL,
                             n_earlystop_rule = c("with_stay", "simple"),
                             p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                             cutoff_eli_resp = 0.90, psi00 = 40, psi11 = 60,
                             benchmark_weight = 0.5, utility_prior = c(1, 1),
                             n_star = 6, n_explore = 9, resp_cor = 0,
                             resp_cor_type = c("copula", "binary"),
                             est_method = c("posterior_mean", "observed",
                                            "tail_probability"),
                             keep_trials = FALSE, verbose = FALSE,
                             seed = 123) {

  call_expr <- match.call()
  n_earlystop_rule <- match.arg(n_earlystop_rule)
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
    results[[i]] <- sim_boin12(
      target = target, target_resp = target_resp,
      p_true = scenarios[[i]]$p_true, p_resp = scenarios[[i]]$p_resp,
      n_cohort = n_cohort, cohort_size = cohort_size, n_trials = n_trials,
      start_dose = start_dose, n_earlystop = n_earlystop,
      n_earlystop_rule = n_earlystop_rule, p_saf = p_saf, p_tox = p_tox,
      cutoff_eli = cutoff_eli, cutoff_eli_resp = cutoff_eli_resp,
      psi00 = psi00, psi11 = psi11, benchmark_weight = benchmark_weight,
      utility_prior = utility_prior, n_star = n_star, n_explore = n_explore,
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
