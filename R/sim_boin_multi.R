#' Operating Characteristics Across Several Scenarios
#'
#' @description
#'   Run \code{\link{sim_boin}} under a list of dose-toxicity scenarios that share
#'   the same design, and collect the results into one table.
#'
#' @param target
#'   Numeric scalar. Target DLT probability.
#'
#' @param scenarios
#'   A list of scenarios. Each element is a list with a \code{name} and a
#'   \code{p_true} vector, or a plain numeric vector of true DLT probabilities in
#'   which case the list names are used as scenario names. Every scenario must
#'   have the same number of doses.
#'
#' @param n_cohort
#'   Integer scalar. Number of cohorts in a trial.
#'
#' @param cohort_size
#'   Integer scalar or vector. Number of patients per cohort.
#'
#' @param n_trials
#'   Integer scalar. Number of trials per scenario. Defaults to 10000.
#'
#' @param start_dose
#'   Integer scalar. Dose level for the first cohort. Defaults to 1. It is
#'   ignored when \code{titration} is \code{TRUE}.
#'
#' @param n_earlystop
#'   Integer scalar. Early stopping sample size at the current dose. Defaults to 18.
#'
#' @param p_saf
#'   Numeric scalar. Highest DLT probability deemed subtherapeutic.
#'   Defaults to \code{0.6 * target}.
#'
#' @param p_tox
#'   Numeric scalar. Lowest DLT probability deemed overly toxic.
#'   Defaults to \code{1.4 * target}.
#'
#' @param cutoff_eli
#'   Numeric scalar. Posterior probability cutoff for dose elimination.
#'   Defaults to 0.95.
#'
#' @param extrasafe
#'   Logical scalar. Apply the stricter safety stopping rule at the lowest dose.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Relaxation of \code{cutoff_eli} for the
#'   safety stopping rule. Defaults to 0.05.
#'
#' @param titration
#'   Logical scalar. Start with single patient cohorts until the first DLT.
#'
#' @param stay_on_1_of_3
#'   Logical scalar. Make one DLT out of three patients a stay rather than a
#'   de-escalation. See \code{\link{boin_boundary}}.
#'
#' @param bound_mtd
#'   Logical scalar. Bound the isotonic estimate at the selected MTD by the
#'   de-escalation boundary.
#'
#' @param mtd_max_estimate
#'   Numeric scalar or \code{NULL}. Largest isotonic estimate a dose may have and
#'   still be selected as the MTD. See \code{\link{sim_boin}}.
#'
#' @param overdose_cutoff
#'   Numeric scalar or \code{NULL}. Doses whose true DLT probability exceeds this
#'   value count as overdoses. Defaults to \code{NULL}, which uses \code{target}.
#'
#' @param min_mtd_sample
#'   Integer scalar. Smallest number of patients for a dose to be eligible.
#'
#' @param n_earlystop_rule
#'   Character scalar, either \code{"with_stay"} or \code{"simple"}.
#'
#' @param keep_trials
#'   Logical scalar. Keep the trial by trial data of every scenario.
#'
#' @param verbose
#'   Logical scalar. Report progress while the simulations run.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed, applied to every scenario, so
#'   that a scenario simulated here matches the same scenario simulated on its own
#'   with \code{\link{sim_boin}}.
#'
#' @return
#'   An object of class \code{boin_oc_multi}, which is a list with components
#'   \item{results}{Named list of \code{boin_oc} objects, one per scenario.}
#'   \item{summary_table}{Data frame collecting all scenarios.}
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
#' oc <- sim_boin_multi(
#'   target = 0.30,
#'   scenarios = scenarios,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_boin}}
#'
#' @export
sim_boin_multi <- function(target, scenarios, n_cohort, cohort_size,
                           n_trials = 10000, start_dose = 1, n_earlystop = 18,
                           p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                           extrasafe = FALSE, offset = 0.05, titration = FALSE,
                           stay_on_1_of_3 = FALSE, bound_mtd = FALSE,
                           mtd_max_estimate = NULL, min_mtd_sample = 1,
                           overdose_cutoff = NULL,
                           n_earlystop_rule = c("with_stay", "simple"),
                           keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  n_earlystop_rule <- match.arg(n_earlystop_rule)

  scenarios <- normalise_scenarios(scenarios)
  scenario_names <- vapply(scenarios, function(z) z$name, character(1))
  n_scenarios <- length(scenarios)
  n_doses <- length(scenarios[[1L]]$p_true)

  results <- vector("list", n_scenarios)
  names(results) <- scenario_names

  for (i in seq_len(n_scenarios)) {
    if (verbose) {
      message("Scenario ", i, " of ", n_scenarios, ": ", scenario_names[i])
    }
    results[[i]] <- sim_boin(
      target = target, p_true = scenarios[[i]]$p_true, n_cohort = n_cohort,
      cohort_size = cohort_size, n_trials = n_trials, start_dose = start_dose,
      n_earlystop = n_earlystop, p_saf = p_saf, p_tox = p_tox,
      cutoff_eli = cutoff_eli, extrasafe = extrasafe, offset = offset,
      titration = titration, stay_on_1_of_3 = stay_on_1_of_3,
      bound_mtd = bound_mtd, mtd_max_estimate = mtd_max_estimate,
      min_mtd_sample = min_mtd_sample, overdose_cutoff = overdose_cutoff,
      n_earlystop_rule = n_earlystop_rule,
      keep_trials = keep_trials, verbose = FALSE, seed = seed
    )
  }

  dose_names <- paste0("DL", seq_len(n_doses))
  item_labels <- c("True DLT rate (%)", "MTD selected (%)",
                   "Patients treated", "Patients with DLT")

  blocks <- lapply(seq_len(n_scenarios), function(i) {
    res <- results[[i]]
    values <- rbind(res$p_true * 100, res$sel_percent,
                    res$n_pts_dose, res$n_tox_dose)
    dimnames(values) <- list(NULL, dose_names)
    block <- data.frame(
      Scenario = c(scenario_names[i], "", "", ""),
      Item = item_labels,
      stringsAsFactors = FALSE
    )
    block <- cbind(block, as.data.frame(round(values, 1)))
    block[["Total / No MTD"]] <- c(
      NA_real_,
      round(res$percent_no_mtd, 1),
      round(res$total_n_pts, 1),
      round(res$total_n_tox, 1)
    )
    block
  })

  summary_table <- do.call(rbind, blocks)
  rownames(summary_table) <- NULL

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
