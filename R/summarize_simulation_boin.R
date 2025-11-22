#' Summarize BOIN Simulation Results
#'
#' @description
#'   Aggregate results from multiple trial simulations into comprehensive
#'   summary statistics. Uses vectorized operations for efficient computation.
#'   Returns an S3 object of class "boin_summary" which has a custom print method.
#'
#' @param simulation_results
#'   List. Results from `simulate_one_trial()`, one element per simulation.
#'   Each element should be a list containing:
#'   - n_pts: numeric vector of patients at each dose
#'   - n_tox: numeric vector of DLTs at each dose
#'   - mtd: numeric value (dose level) or NA
#'
#' @param n_doses
#'   Numeric. Number of doses evaluated in the simulations.
#'
#' @return
#'   An object of class "boin_summary" (a list) containing:
#'   \item{mtd_selection_percent}{Numeric vector. Percentage of trials selecting
#'                                 each dose as MTD}
#'   \item{avg_n_pts}{Numeric vector. Average number of patients treated at each dose}
#'   \item{avg_n_tox}{Numeric vector. Average number of DLTs at each dose}
#'   \item{percent_no_mtd}{Numeric. Percentage of trials where no MTD was selected}
#'   \item{avg_total_n_pts}{Numeric. Average total patients across all doses}
#'   \item{avg_total_n_tox}{Numeric. Average total DLTs across all doses}
#'
#' @details
#'   This function is typically called after running multiple trial simulations
#'   via `sim_boin()`. It computes operating characteristics including MTD
#'   selection rates, patient allocation, and safety metrics.
#'
#'   Operating characteristics are key for evaluating trial design performance
#'   and are presented in trial protocol applications.
#'
#'   The returned object is of class "boin_summary" which has a custom print method.
#'   Simply calling `print(result)` or typing `result` will display a formatted
#'   summary table.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' # This function is typically used internally by sim_boin()
#' # For demonstration of aggregated results:
#' target <- 0.30
#' p_true <- c(0.10, 0.25, 0.40)
#' boin_bound <- get_boin_boundary(target)
#' decision_table <- get_boin_decision(target, boin_bound$lambda_e,
#'                                 boin_bound$lambda_d, 18, 0.95)
#' stopping_boundaries <- get_boin_stopping_boundaries(target, 18, 0.90)
#' result <- sim_boin(
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
#'
#' # S3 method automatically formats output
#' print(result$summary)
#'
#' # Or simply
#' result$summary
#'
#' # For knitr::kable format
#' print(result$summary, kable_output = TRUE)
#' }
#'
#' @export
summarize_simulation_boin <- function(simulation_results, n_doses) {

  n_trials <- length(simulation_results)

  # Extract data directly into matrices
  # Use vapply for type-safe extraction
  n_pts_all <- vapply(simulation_results, function(x) x$n_pts, numeric(n_doses))
  n_tox_all <- vapply(simulation_results, function(x) x$n_tox, numeric(n_doses))
  mtd_vector <- vapply(simulation_results, function(x) x$mtd, numeric(1), USE.NAMES = FALSE)

  # Transpose matrices to have trials in rows, doses in columns
  n_pts_all <- t(n_pts_all)
  n_tox_all <- t(n_tox_all)

  # Compute MTD selection matrix vectorized
  # Initialize as numeric, then identify which doses were selected
  mtd_selected <- matrix(0, nrow = n_trials, ncol = n_doses)
  valid_mtd <- !is.na(mtd_vector)
  mtd_selected[cbind(which(valid_mtd), mtd_vector[valid_mtd])] <- 1

  # Compute all statistics in vectorized form
  summary <- list(
    # MTD selection percentage
    mtd_selection_percent = colMeans(mtd_selected) * 100,

    # Average number of patients per dose
    avg_n_pts = colMeans(n_pts_all),

    # Average number of DLTs per dose
    avg_n_tox = colMeans(n_tox_all),

    # Percentage with no MTD selected
    percent_no_mtd = (1 - mean(valid_mtd)) * 100,

    # Average total patients across all doses
    avg_total_n_pts = mean(rowSums(n_pts_all)),

    # Average total DLTs across all doses
    avg_total_n_tox = mean(rowSums(n_tox_all))
  )

  # Set S3 class
  class(summary) <- c("boin_summary", "list")

  return(summary)
}
