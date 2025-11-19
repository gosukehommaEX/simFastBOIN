#' Run BOIN Simulation
#'
#' @description
#'   Execute multiple BOIN trial simulations and return aggregated operating
#'   characteristics. This is the main entry point for running simulation studies
#'   to evaluate BOIN design performance.
#'
#' @param n_trials
#'   Numeric. Number of trials to simulate. Default is 10000.
#'   Typically 1000-10000 for detailed operating characteristics.
#'
#' @param target
#'   Numeric. The target toxicity probability (e.g., 0.30 for 30%).
#'
#' @param p_true
#'   Numeric vector. True toxicity probabilities for each dose.
#'   Length determines number of doses evaluated.
#'
#' @param n_doses
#'   Numeric. Number of doses evaluated.
#'
#' @param n_cohort
#'   Numeric. Maximum number of cohorts per trial.
#'
#' @param cohort_size
#'   Numeric vector or scalar specifying patients per cohort.
#'   If vector (e.g., c(4, 3, 3)), each element specifies size for corresponding cohort.
#'   If scalar, all cohorts use the same size.
#'
#' @param lambda_e
#'   Numeric. Escalation boundary from `get_boin_boundary()`.
#'
#' @param lambda_d
#'   Numeric. De-escalation boundary from `get_boin_boundary()`.
#'
#' @param decision_table
#'   Character matrix. Decision table from `get_boin_decision()`.
#'
#' @param stopping_boundaries
#'   Character matrix. Trial stopping rule table from `get_boin_stopping_boundaries()`.
#'
#' @param n_earlystop
#'   Numeric. Sample size triggering early stopping. Default is 18.
#'
#' @param min_mtd_sample
#'   Numeric. Minimum sample size for MTD consideration. Default is 6.
#'
#' @param seed
#'   Numeric. Random seed for reproducibility. Default is 123.
#'
#' @return
#'   A list containing:
#'   \item{detailed_results}{List of results from each individual trial}
#'   \item{summary}{Aggregated summary statistics from `summarize_simulation()`}
#'
#' @details
#'   This function orchestrates the entire simulation workflow:
#'   1. Sets random seed for reproducibility
#'   2. Runs n_trials independent trial simulations
#'   3. Aggregates results into summary statistics
#'   4. Returns both detailed and summary results
#'
#'   The detailed results allow for further custom analysis if needed.
#'   Summary statistics include MTD selection rates, patient allocation,
#'   and safety metrics needed for protocol development.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' # Example: Simulation for 3-dose study with 30% target DLT rate
#' target <- 0.30
#' p_true <- c(0.10, 0.25, 0.40)
#' boin_bound <- get_boin_boundary(target)
#' decision_table <- get_boin_decision(target, boin_bound$lambda_e,
#'                                 boin_bound$lambda_d, 18, 0.95)
#' stopping_boundaries <- get_boin_stopping_boundaries(target, 18, 0.90)
#'
#' result <- sim_boin_multiple_trial(
#'   n_trials = 1000,
#'   target = target,
#'   p_true = p_true,
#'   n_doses = 3,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   decision_table = decision_table,
#'   stopping_boundaries = stopping_boundaries,
#'   seed = 123
#' )
#' print_simulation_summary(result$summary, n_doses = 3)
#' }
#'
#' @export
sim_boin_multiple_trial <- function(
    n_trials = 10000,
    target,
    p_true,
    n_doses,
    n_cohort,
    cohort_size,
    decision_table,
    stopping_boundaries,
    n_earlystop = 18,
    min_mtd_sample = 1,
    seed = 123
) {

  set.seed(seed)

  cat("========================================\n")
  cat("Starting BOIN Simulation\n")
  cat("Number of trials:", n_trials, "\n")
  cat("Target DLT rate:", target * 100, "%\n")
  cat("Number of doses:", n_doses, "\n")
  cat("========================================\n\n")

  # Use list for storing results
  simulation_results <- vector("list", n_trials)

  # Determine progress reporting interval
  progress_interval <- max(1, n_trials %/% 10)  # Report 10 times during simulation

  # Main simulation loop
  for (trial in seq_len(n_trials)) {

    # Progress reporting (optimized to occur less frequently)
    if (trial %% progress_interval == 0) {
      cat("Progress: ", trial, " / ", n_trials, " trials completed\n", sep = "")
    }

    simulation_results[[trial]] <- sim_boin_one_trial(
      target = target,
      p_true = p_true,
      n_doses = n_doses,
      n_cohort = n_cohort,
      cohort_size = cohort_size,
      decision_table = decision_table,
      stopping_boundaries = stopping_boundaries,
      n_earlystop = n_earlystop,
      min_mtd_sample = min_mtd_sample
    )
  }

  cat("\nSimulation completed!\n\n")

  # Compute summary statistics
  summary_result <- summarize_simulation_boin(simulation_results, n_doses)

  return(list(
    detailed_results = simulation_results,
    summary = summary_result
  ))
}
