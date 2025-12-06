#' Run BOIN Simulation for Multiple Scenarios
#'
#' @description
#'   Execute BOIN trial simulations across multiple dose-toxicity scenarios and
#'   compute operating characteristics for each scenario. Results are aggregated
#'   into a combined summary table suitable for protocol development and reports.
#'
#' @param scenarios
#'   List of scenario definitions. Each element must be a list containing:
#'   \describe{
#'     \item{\code{name}}{Character. Scenario identifier (e.g., "Scenario 1")}
#'     \item{\code{p_true}}{Numeric vector. True toxicity probabilities for each dose}
#'   }
#'
#' @param target
#'   Numeric. Target toxicity probability (e.g., 0.30 for 30%).
#'
#' @param p_saf
#'   Numeric. Highest toxicity probability deemed acceptable for safety.
#'   Default is 0.6 * target. Used with p_tox for safety/efficacy dose identification.
#'
#' @param p_tox
#'   Numeric. Lowest toxicity probability deemed unacceptable for toxicity.
#'   Default is 1.4 * target. Used with p_saf for safety/efficacy dose identification.
#'
#' @param n_trials
#'   Numeric. Number of trials to simulate per scenario. Default is 10000.
#'
#' @param n_cohort
#'   Numeric. Maximum number of cohorts per trial.
#'
#' @param cohort_size
#'   Numeric vector or scalar. Patients per cohort.
#'
#' @param n_earlystop
#'   Numeric. Sample size triggering early stopping. Default is 18.
#'
#' @param cutoff_eli
#'   Numeric. Cutoff probability for dose elimination. Default is 0.95.
#'
#' @param extrasafe
#'   Logical. Apply extra safety stopping rule at lowest dose. Default is FALSE.
#'
#' @param offset
#'   Numeric. Offset for safety cutoff when extrasafe = TRUE. Default is 0.05.
#'
#' @param n_earlystop_rule
#'   Character. Early stopping rule: "with_stay" or "simple". Default is "with_stay".
#'
#' @param titration
#'   Logical. Perform accelerated dose titration phase. Default is FALSE.
#'
#' @param min_mtd_sample
#'   Numeric. Minimum patients required for MTD consideration. Default is 1.
#'
#' @param boundMTD
#'   Logical. Impose constraint that MTD's isotonic estimate <= lambda_d. Default is FALSE.
#'
#' @param return_details
#'   Logical. If TRUE, return detailed trial-level results. Default is FALSE.
#'
#' @param verbose
#'   Logical. If TRUE, print progress messages to console. If FALSE, run silently. Default is FALSE.
#'
#' @param seed
#'   Numeric. Random seed for reproducibility. Default is 123.
#'
#' @return
#'   A list with class "boin_multi_summary" containing:
#'   \item{results_by_scenario}{List of simulation results for each scenario}
#'   \item{combined_summary_df}{Data frame with aggregated operating characteristics}
#'   \item{scenario_names}{Character vector of scenario identifiers}
#'   \item{n_doses}{Numeric. Number of doses evaluated}
#'   \item{call}{The function call}
#'
#' @details
#'   Progress messages are printed to console only if verbose = TRUE.
#'   For large-scale studies (50+ scenarios), computation time may be substantial.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' scenarios <- list(
#'   list(name = "Scenario 1", p_true = c(0.05, 0.10, 0.20, 0.30, 0.45)),
#'   list(name = "Scenario 2", p_true = c(0.10, 0.15, 0.30, 0.45, 0.60))
#' )
#'
#' # Silent mode (default)
#' result <- sim_boin_multi(
#'   scenarios = scenarios,
#'   target = 0.30,
#'   n_trials = 10000,
#'   n_cohort = 48,
#'   cohort_size = 3,
#'   seed = 123
#' )
#'
#' # With progress messages
#' result_verbose <- sim_boin_multi(
#'   scenarios = scenarios,
#'   target = 0.30,
#'   n_trials = 10000,
#'   n_cohort = 48,
#'   cohort_size = 3,
#'   verbose = TRUE,
#'   seed = 123
#' )
#' }
#'
#' @importFrom utils head tail
#'
#' @export
sim_boin_multi <- function(scenarios,
                           target,
                           n_trials = 10000,
                           n_cohort,
                           cohort_size,
                           p_saf = NULL,
                           p_tox = NULL,
                           n_earlystop = 18,
                           cutoff_eli = 0.95,
                           extrasafe = FALSE,
                           offset = 0.05,
                           n_earlystop_rule = "with_stay",
                           titration = FALSE,
                           min_mtd_sample = 1,
                           boundMTD = FALSE,
                           return_details = FALSE,
                           verbose = FALSE,
                           seed = 123) {

  call_expr <- match.call()

  # Set default values for p_saf and p_tox if not provided
  if (is.null(p_saf)) {
    p_saf <- 0.6 * target
  }
  if (is.null(p_tox)) {
    p_tox <- 1.4 * target
  }

  # ===== INPUT VALIDATION =====
  if (!is.list(scenarios) || length(scenarios) == 0) {
    stop("'scenarios' must be a non-empty list of scenario definitions")
  }

  for (i in seq_along(scenarios)) {
    if (!is.list(scenarios[[i]]) ||
        !("name" %in% names(scenarios[[i]])) ||
        !("p_true" %in% names(scenarios[[i]]))) {
      stop("Each scenario must be a list with 'name' and 'p_true' elements")
    }
  }

  scenario_names <- sapply(scenarios, function(x) x$name)
  n_scenarios <- length(scenarios)
  n_doses <- length(scenarios[[1]]$p_true)

  for (i in seq_along(scenarios)) {
    if (length(scenarios[[i]]$p_true) != n_doses) {
      stop("All scenarios must have the same number of doses")
    }
  }

  # ===== SIMULATION LOOP =====
  if (verbose) {
    cat("========================================\n")
    cat("BOIN Multi-Scenario Simulation\n")
    cat("Number of scenarios:", n_scenarios, "\n")
    cat("Trials per scenario:", n_trials, "\n")
    cat("Number of doses:", n_doses, "\n")
    cat("========================================\n\n")
  }

  results_by_scenario <- list()

  for (i in seq_len(n_scenarios)) {
    scenario <- scenarios[[i]]

    if (verbose) {
      cat("Processing", scenario$name, "(", i, "of", n_scenarios, ")...\n")
    }

    result <- sim_boin(
      n_trials = n_trials,
      target = target,
      p_true = scenario$p_true,
      n_cohort = n_cohort,
      cohort_size = cohort_size,
      p_saf = p_saf,
      p_tox = p_tox,
      n_earlystop = n_earlystop,
      cutoff_eli = cutoff_eli,
      extrasafe = extrasafe,
      offset = offset,
      n_earlystop_rule = n_earlystop_rule,
      titration = titration,
      min_mtd_sample = min_mtd_sample,
      boundMTD = boundMTD,
      return_details = return_details,
      verbose = verbose,  # Pass verbose parameter to sim_boin
      seed = seed + i
    )

    results_by_scenario[[scenario$name]] <- result

    if (verbose) cat("\n")
  }

  # ===== AGGREGATE RESULTS INTO DATA FRAME =====
  if (verbose) cat("Aggregating results into summary table...\n")

  table_data <- data.frame()

  for (scenario_name in scenario_names) {
    res <- results_by_scenario[[scenario_name]]
    p_true_scenario <- res$summary$p_true
    dose_values <- round(p_true_scenario * 100, 1)

    # Row 1: True DLT Rate (%)
    row_data <- c(scenario_name, "True DLT Rate (%)", as.character(dose_values), "")
    table_data <- rbind(table_data, row_data)

    # Row 2: MTD Selected (%)
    mtd_sel_vec <- head(res$summary$mtd_selection_percent, -1)
    pct_no_mtd <- tail(res$summary$mtd_selection_percent, 1)
    mtd_values <- round(c(mtd_sel_vec, pct_no_mtd), 1)
    row_data <- c("", "MTD Selected (%)", as.character(mtd_values[seq_len(n_doses)]),
                  as.character(mtd_values[n_doses + 1]))
    table_data <- rbind(table_data, row_data)

    # Row 3: Number of Participants Treated
    pts_values <- round(c(res$summary$avg_n_pts, res$summary$avg_total_n_pts), 1)
    row_data <- c("", "Number of Participants Treated", as.character(pts_values[seq_len(n_doses)]),
                  as.character(pts_values[n_doses + 1]))
    table_data <- rbind(table_data, row_data)

    # Row 4: Number of Participants w/ DLTs
    tox_values <- round(c(res$summary$avg_n_tox, res$summary$avg_total_n_tox), 1)
    row_data <- c("", "Number of Participants w/ DLTs", as.character(tox_values[seq_len(n_doses)]),
                  as.character(tox_values[n_doses + 1]))
    table_data <- rbind(table_data, row_data)
  }

  dose_col_names <- paste0("DL", seq_len(n_doses))
  total_col_name <- "Total/No MTD"
  colnames(table_data) <- c("Scenario", "Item", dose_col_names, total_col_name)

  if (verbose) cat("  Completed.\n\n")

  # ===== CONSTRUCT RETURN OBJECT =====
  return_obj <- list(
    results_by_scenario = results_by_scenario,
    combined_summary_df = table_data,
    scenario_names = scenario_names,
    n_doses = n_doses,
    call = call_expr
  )

  class(return_obj) <- c("boin_multi_summary", "list")

  if (verbose) {
    cat("========================================\n")
    cat("Multi-Scenario Summary\n")
    cat("========================================\n\n")
    print(return_obj)
    cat("\n")
  }

  if (verbose) {
    return(return_obj)
  } else {
    invisible(return_obj)
  }
}
