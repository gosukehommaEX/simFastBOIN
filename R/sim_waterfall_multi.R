#' Operating Characteristics of the Waterfall Design Across Scenarios
#'
#' @description
#'   Run \code{\link{sim_waterfall}} under several matrices of true DLT
#'   probabilities with the same design, and collect the main operating
#'   characteristics of every scenario in one table.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param scenarios
#'   List of scenarios. Each element is either a matrix of true DLT
#'   probabilities, with one row per dose level of drug A and one column per
#'   dose level of drug B, or a list with elements \code{name} and
#'   \code{p_true}. Names of the list are used as scenario names. All matrices
#'   must have the same dimensions.
#'
#' @param cohort_size
#'   Integer scalar. Number of patients per cohort.
#'
#' @param n_cohort
#'   Integer vector with one element per dose level of drug A, or \code{NULL}.
#'   Maximum number of cohorts of the first, second, ... subtrial. Defaults to
#'   \code{NULL}, which gives every subtrial about four patients per dose of
#'   the subtrial (see Details).
#'
#' @param n_trials
#'   Integer scalar. Number of trials to simulate. Defaults to 10000.
#'
#' @param n_earlystop
#'   Integer scalar. A subtrial ends once this many patients have been treated
#'   at the dose chosen for its next cohort. Defaults to 12.
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
#'   Numeric scalar. Posterior probability cutoff for eliminating a combination.
#'   Defaults to 0.95.
#'
#' @param extrasafe
#'   Logical scalar. Apply the stricter safety stopping rule at the lowest dose
#'   of every subtrial. Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is
#'   relaxed for the safety stopping rule. Defaults to 0.05.
#'
#' @param titration
#'   Logical scalar. Start the first subtrial with one patient per dose, up to
#'   the first DLT. Ignored when \code{cohort_size} is one.
#'   Defaults to \code{FALSE}.
#'
#' @param bound_mtd
#'   Logical scalar. Require the estimate at the candidate MTD of a subtrial to
#'   be below the de-escalation boundary. Defaults to \code{FALSE}.
#'
#' @param boin_compatible
#'   Logical scalar. Drop the data of a subtrial that ends without a candidate
#'   MTD, as \code{BOIN::get.oc.comb()} does (see Details).
#'   Defaults to \code{FALSE}.
#'
#' @param keep_trials
#'   Logical scalar. Keep the full trial by trial data in the result.
#'   Defaults to \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Report progress while the simulation runs.
#'   Defaults to \code{FALSE}.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. The state of the calling session
#'   is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{waterfall_oc_multi}, which is a list with components
#'   \item{results}{Named list of \code{waterfall_oc} objects, one per scenario.}
#'   \item{summary_table}{Data frame with one row per scenario: the true MTD contour, the percentage of trials selecting exactly that contour, the percentages of patients treated on and above it, the percentage of trials without any MTD, and the average numbers of patients and DLTs.}
#'   \item{scenario_names}{Character vector of scenario names.}
#'   \item{n_rows, n_cols}{Numbers of dose levels of drug A and drug B.}
#'   \item{call}{The matched call.}
#'
#' @details
#'   Every scenario is simulated with the same seed, so that the scenarios
#'   differ by their DLT probabilities and not by the random numbers.
#'
#' @references
#'   Zhang, L. and Yuan, Y. (2016). A Practical Bayesian Design to Identify the
#'   Maximum Tolerated Dose Contour for Drug Combination Trials. Statistics in
#'   Medicine, 35(27), 4924-4936.
#'
#' @examples
#' scenarios <- list(
#'   "Contour A1B3, A2B2" = matrix(c(0.03, 0.10, 0.28,
#'                                   0.10, 0.30, 0.50), nrow = 2, byrow = TRUE),
#'   "Contour A1B1" = matrix(c(0.30, 0.40, 0.50,
#'                             0.42, 0.49, 0.55), nrow = 2, byrow = TRUE)
#' )
#'
#' oc <- sim_waterfall_multi(
#'   target = 0.30,
#'   scenarios = scenarios,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#' oc
#'
#' # The full results of one scenario
#' oc$results[["Contour A1B1"]]
#'
#' @seealso \code{\link{sim_waterfall}}
#'
#' @export
sim_waterfall_multi <- function(target, scenarios, cohort_size,
                                n_cohort = NULL, n_trials = 10000,
                                n_earlystop = 12, p_saf = NULL, p_tox = NULL,
                                cutoff_eli = 0.95, extrasafe = FALSE,
                                offset = 0.05, titration = FALSE,
                                bound_mtd = FALSE, boin_compatible = FALSE,
                                keep_trials = FALSE, verbose = FALSE,
                                seed = 123) {

  call_expr <- match.call()

  # normalize_comb_scenarios() also checks that the scenarios have the same
  # dimensions and unique names.
  scenarios <- normalize_comb_scenarios(scenarios)
  scenario_names <- vapply(scenarios, function(z) z$name, character(1))
  n_scenarios <- length(scenarios)
  first_dim <- dim(scenarios[[1L]]$p_true)

  results <- vector("list", n_scenarios)
  names(results) <- scenario_names

  for (i in seq_len(n_scenarios)) {
    if (verbose) {
      message("Scenario ", i, " of ", n_scenarios, ": ", scenario_names[i])
    }
    results[[i]] <- sim_waterfall(
      target = target, p_true = scenarios[[i]]$p_true,
      cohort_size = cohort_size, n_cohort = n_cohort, n_trials = n_trials,
      n_earlystop = n_earlystop, p_saf = p_saf, p_tox = p_tox,
      cutoff_eli = cutoff_eli, extrasafe = extrasafe, offset = offset,
      titration = titration, bound_mtd = bound_mtd,
      boin_compatible = boin_compatible, keep_trials = keep_trials,
      verbose = FALSE, seed = seed
    )
  }

  structure(
    list(
      results = results,
      summary_table = waterfall_oc_multi_table(results, scenario_names),
      scenario_names = scenario_names,
      n_rows = first_dim[1L],
      n_cols = first_dim[2L],
      call = call_expr
    ),
    class = "waterfall_oc_multi"
  )
}
