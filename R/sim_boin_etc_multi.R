#' Operating Characteristics of the BOIN-ETC Design Across Scenarios
#'
#' @description
#'   Run \code{\link{sim_boin_etc}} for several scenarios of DLT and response
#'   probabilities with the same design and collect the results in one
#'   table.
#'
#' @inheritParams sim_boin_etc
#'
#' @param scenarios
#'   List of scenarios. Each element is a list with the matrices
#'   \code{p_true} and \code{p_resp} (one row per level of drug A, one column
#'   per level of drug B), and optionally \code{name}. Names of the list are
#'   used when \code{name} is absent. All scenarios must have the same
#'   numbers of levels.
#'
#' @return
#'   An object of class \code{boin_etc_oc_multi}, which is a list with
#'   components
#'   \item{results}{Named list of \code{boin_etc_oc} objects, one per
#'     scenario.}
#'   \item{summary_table}{Data frame with one row per scenario, as shown by
#'     \code{\link{print.boin_etc_oc_multi}}.}
#'   \item{scenario_names}{Character vector of scenario names.}
#'   \item{n_rows, n_cols}{Numbers of levels of drug A and drug B.}
#'   \item{call}{The matched call.}
#'
#' @details
#'   Every scenario is simulated with the same \code{seed}, so that a scenario
#'   simulated here matches the same scenario simulated on its own with
#'   \code{\link{sim_boin_etc}}.
#'
#' @examples
#' # Scenarios 1 and 2 of Kakizume et al. (2024)
#' scenarios <- list(
#'   "Scenario 1" = list(
#'     p_true = matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2),
#'     p_resp = matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)),
#'   "Scenario 2" = list(
#'     p_true = matrix(c(0.05, 0.15, 0.15, 0.35, 0.30, 0.55), nrow = 2),
#'     p_resp = matrix(c(0.25, 0.60, 0.60, 0.65, 0.65, 0.70), nrow = 2))
#' )
#' res <- sim_boin_etc_multi(target = 0.30, target_resp = 0.60,
#'                           scenarios = scenarios, n_cohort = 8,
#'                           cohort_size = 3, n_trials = 200)
#' res
#'
#' @seealso \code{\link{sim_boin_etc}}, \code{\link{print.boin_etc_oc_multi}}
#'
#' @export
sim_boin_etc_multi <- function(target, target_resp, scenarios, n_cohort,
                               cohort_size, n_trials = 10000,
                               subtrials = c("etc1", "etc2", "etc3"),
                               n_earlystop = 6,
                               etc1_start = c("diagonal", "below"),
                               rule4 = c("paper", "authors_code"),
                               p_saf = NULL, p_tox = NULL,
                               p_saf_resp = NULL, prior = rep(1 / 6, 6),
                               cutoff_eli = 0.95, p_eli = NULL, resp_cor = 0,
                               resp_cor_type = c("copula", "binary"),
                               psi00 = 40, psi11 = 60,
                               utility_prior = c(1, 1),
                               odc_tie = c("lowest", "all"),
                               keep_trials = FALSE, verbose = FALSE,
                               seed = 123) {

  call_expr <- match.call()
  subtrials <- match.arg(subtrials)
  etc1_start <- match.arg(etc1_start)
  rule4 <- match.arg(rule4)
  resp_cor_type <- match.arg(resp_cor_type)
  odc_tie <- match.arg(odc_tie)
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
    results[[i]] <- sim_boin_etc(
      target = target, target_resp = target_resp,
      p_true = scenarios[[i]]$p_true, p_resp = scenarios[[i]]$p_resp,
      n_cohort = n_cohort, cohort_size = cohort_size, n_trials = n_trials,
      subtrials = subtrials, n_earlystop = n_earlystop,
      etc1_start = etc1_start, rule4 = rule4, p_saf = p_saf, p_tox = p_tox,
      p_saf_resp = p_saf_resp, prior = prior, cutoff_eli = cutoff_eli,
      p_eli = p_eli, resp_cor = resp_cor, resp_cor_type = resp_cor_type,
      psi00 = psi00, psi11 = psi11, utility_prior = utility_prior,
      odc_tie = odc_tie, keep_trials = keep_trials, verbose = FALSE,
      seed = seed
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
    class = "boin_etc_oc_multi"
  )
}
