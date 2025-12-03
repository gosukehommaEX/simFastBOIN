#' Summarize BOIN Simulation Results
#'
#' @description
#'   Aggregate results from multiple trial simulations into comprehensive
#'   summary statistics. Uses vectorized operations for efficient computation.
#'   Returns an S3 object of class "boin_summary" which has a custom print method.
#'
#' @param simulation_results List. Results from \code{simulate_one_trial()}, one element per simulation.
#'   Each element should be a list containing:
#'   \itemize{
#'     \item n_pts: numeric vector of patients at each dose
#'     \item n_tox: numeric vector of DLTs at each dose
#'     \item mtd: numeric value (dose level) or NA
#'   }
#' @param p_true Numeric vector. True toxicity probabilities for each dose.
#'   Length determines number of doses evaluated.
#'
#' @return An object of class "boin_summary" (a list) containing:
#'   \item{p_true}{Numeric vector. True toxicity probabilities}
#'   \item{mtd_selection_percent}{Numeric vector (length n_doses + 1). Percentage of trials
#'     selecting each dose as MTD. The last element is the percentage with no MTD selected.}
#'   \item{avg_n_pts}{Numeric vector. Average number of patients treated at each dose}
#'   \item{avg_n_tox}{Numeric vector. Average number of DLTs at each dose}
#'   \item{avg_total_n_pts}{Numeric. Average total patients across all doses}
#'   \item{avg_total_n_tox}{Numeric. Average total DLTs across all doses}
#'
#' @details
#'   This function is typically called after running multiple trial simulations
#'   via \code{sim_boin()}. It computes operating characteristics including MTD
#'   selection rates, patient allocation, and safety metrics.
#'
#'   Operating characteristics are key for evaluating trial design performance
#'   and are presented in trial protocol applications.
#'
#'   The returned object is of class "boin_summary" which has a custom print method.
#'   Simply calling \code{print(result)} or typing \code{result} will display a formatted
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
#'
#' result <- sim_boin(
#'   n_trials = 1000,
#'   target = target,
#'   p_true = p_true,
#'   n_cohort = 10,
#'   cohort_size = 3,
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
summarize_simulation_boin <- function(simulation_results, p_true) {

  n_trials <- length(simulation_results)
  n_doses <- length(p_true)

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

  # Compute MTD selection percentages including "No MTD" as the last column
  mtd_selection_by_dose <- colMeans(mtd_selected) * 100
  percent_no_mtd <- (1 - mean(valid_mtd)) * 100

  # Combine into single vector: [DL1, DL2, ..., DLn, No MTD]
  mtd_selection_percent <- c(mtd_selection_by_dose, percent_no_mtd)

  # Compute all statistics in vectorized form
  summary <- list(
    # True toxicity probabilities
    p_true = p_true,

    # MTD selection percentage (including No MTD as last element)
    mtd_selection_percent = mtd_selection_percent,

    # Average number of patients per dose
    avg_n_pts = colMeans(n_pts_all),

    # Average number of DLTs per dose
    avg_n_tox = colMeans(n_tox_all),

    # Average total patients across all doses
    avg_total_n_pts = mean(rowSums(n_pts_all)),

    # Average total DLTs across all doses
    avg_total_n_tox = mean(rowSums(n_tox_all))
  )

  # Set S3 class
  class(summary) <- c("boin_summary", "list")

  return(summary)
}
