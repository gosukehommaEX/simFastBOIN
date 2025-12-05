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
#'   Example: \code{list(name = "Scenario 1", p_true = c(0.05, 0.10, 0.15, 0.20))}
#'
#' @param target
#'   Numeric. Target toxicity probability (e.g., 0.30 for 30%).
#'
#' @param n_trials
#'   Numeric. Number of trials to simulate per scenario. Default is 10000.
#'   Larger values (e.g., 10000) yield more stable operating characteristics
#'   at the cost of increased computation time.
#'
#' @param n_cohort
#'   Numeric. Maximum number of cohorts per trial.
#'
#' @param cohort_size
#'   Numeric vector or scalar. Patients per cohort. If scalar, all cohorts
#'   use the same size. If vector, each element specifies size for corresponding cohort.
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
#'   Logical. If TRUE, return detailed trial-level results for each scenario.
#'   If FALSE (default), return summary statistics only for memory efficiency.
#'
#' @param seed
#'   Numeric. Random seed for reproducibility. Default is 123.
#'
#' @return
#'   A list with class "boin_multi_results" containing:
#'   \item{results_by_scenario}{List of length equal to number of scenarios.
#'     Each element contains the output from \code{sim_boin()} for that scenario.}
#'   \item{combined_summary_df}{Data frame with aggregated operating characteristics
#'     across all scenarios. Rows are organized as:
#'     True DLT Rate (%), MTD Selected (%), Number of Participants Treated,
#'     and Number of Participants w/ DLTs.
#'     Columns represent doses plus "N/S Total" for overall statistics.}
#'   \item{scenario_names}{Character vector of scenario identifiers.}
#'   \item{n_doses}{Numeric. Number of doses evaluated (inferred from first scenario).}
#'   \item{call}{The function call as entered by the user.}
#'
#' @details
#'   This function orchestrates simulations across multiple scenarios by:
#'   \enumerate{
#'     \item Validating scenario specifications
#'     \item Running \code{\link{sim_boin}} for each scenario independently
#'     \item Aggregating results into a unified data frame
#'     \item Organizing output for easy comparison across scenarios
#'   }
#'
#'   Operating characteristics computed for each scenario include:
#'   \describe{
#'     \item{True DLT Rate}{Actual dose-toxicity relationship}
#'     \item{MTD Selected}{Percentage selecting each dose as MTD (+ \% with no MTD)}
#'     \item{Number of Participants Treated}{Average enrollment per dose}
#'     \item{Number of Participants w/ DLTs}{Average DLT count per dose}
#'   }
#'
#'   Progress messages are printed to console for monitoring simulation status.
#'   For large-scale studies (50+ scenarios), computation time may be substantial.
#'   Use \code{return_details = FALSE} (default) for faster execution.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' # Define multiple scenarios with different dose-toxicity relationships
#' scenarios <- list(
#'   list(name = "Scenario 1: Linear Increase",
#'        p_true = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)),
#'   list(name = "Scenario 2: Steep Early",
#'        p_true = c(0.01, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.30)),
#'   list(name = "Scenario 3: Plateau",
#'        p_true = c(0.15, 0.20, 0.25, 0.30, 0.45, 0.50, 0.50, 0.50, 0.50))
#' )
#'
#' # Run simulations across all scenarios
#' result <- sim_boin_multi(
#'   scenarios = scenarios,
#'   target = 0.30,
#'   n_trials = 10000,
#'   n_cohort = 48,
#'   cohort_size = 3,
#'   seed = 123
#' )
#'
#' # View aggregated results (plain text)
#' print(result)
#'
#' # View as kable Markdown table
#' print(result, kable = TRUE, kable_format = "pipe")
#'
#' # Access scenario-specific results
#' result$results_by_scenario[["Scenario 1: Linear Increase"]]
#' }
#'
#' @import stats
#'
#' @export
sim_boin_multi <- function(scenarios,
                           target,
                           n_trials = 10000,
                           n_cohort,
                           cohort_size,
                           n_earlystop = 18,
                           cutoff_eli = 0.95,
                           extrasafe = FALSE,
                           offset = 0.05,
                           n_earlystop_rule = "with_stay",
                           titration = FALSE,
                           min_mtd_sample = 1,
                           boundMTD = FALSE,
                           return_details = FALSE,
                           seed = 123) {

  # Store the function call for output object
  call_expr <- match.call()

  # ===== INPUT VALIDATION =====
  # Validate scenarios parameter
  if (!is.list(scenarios) || length(scenarios) == 0) {
    stop("'scenarios' must be a non-empty list of scenario definitions")
  }

  # Check each scenario has required fields
  for (i in seq_along(scenarios)) {
    if (!is.list(scenarios[[i]]) ||
        !("name" %in% names(scenarios[[i]])) ||
        !("p_true" %in% names(scenarios[[i]]))) {
      stop("Each scenario must be a list with 'name' and 'p_true' elements")
    }
  }

  # Extract scenario names and dose-toxicity vectors
  scenario_names <- sapply(scenarios, function(x) x$name)
  n_scenarios <- length(scenarios)
  n_doses <- length(scenarios[[1]]$p_true)

  # Verify all scenarios have same number of doses
  for (i in seq_along(scenarios)) {
    if (length(scenarios[[i]]$p_true) != n_doses) {
      stop("All scenarios must have the same number of doses")
    }
  }

  # ===== SIMULATION LOOP =====
  # Print header
  cat("========================================\n")
  cat("BOIN Multi-Scenario Simulation\n")
  cat("Number of scenarios:", n_scenarios, "\n")
  cat("Trials per scenario:", n_trials, "\n")
  cat("Number of doses:", n_doses, "\n")
  cat("========================================\n\n")

  # Initialize list to store results for each scenario
  results_by_scenario <- list()

  # Execute sim_boin for each scenario
  for (i in seq_len(n_scenarios)) {
    scenario <- scenarios[[i]]

    cat("Processing", scenario$name, "(", i, "of", n_scenarios, ")...\n")

    # Run BOIN simulation for this scenario
    result <- sim_boin(
      n_trials = n_trials,
      target = target,
      p_true = scenario$p_true,
      n_cohort = n_cohort,
      cohort_size = cohort_size,
      n_earlystop = n_earlystop,
      cutoff_eli = cutoff_eli,
      extrasafe = extrasafe,
      offset = offset,
      n_earlystop_rule = n_earlystop_rule,
      titration = titration,
      min_mtd_sample = min_mtd_sample,
      boundMTD = boundMTD,
      return_details = return_details,
      seed = seed + i
    )

    # Store result with scenario name as key
    results_by_scenario[[scenario$name]] <- result

    cat("\n")
  }

  # ===== AGGREGATE RESULTS INTO DATA FRAME =====
  cat("Aggregating results into summary table...\n")

  # Initialize empty data frame for combined results
  table_data <- data.frame()

  # Process each scenario and add rows to table
  for (scenario_name in scenario_names) {
    res <- results_by_scenario[[scenario_name]]

    # Extract true p_true vector
    p_true_scenario <- res$summary$p_true

    # Create dose columns (DL1, DL2, ..., DLn)
    dose_values <- round(p_true_scenario * 100, 1)

    # Row 1: True DLT Rate (%)
    row_data <- c(scenario_name, "True DLT Rate (%)", as.character(dose_values), "")
    table_data <- rbind(table_data, row_data)

    # Row 2: MTD Selected (%)
    # Last element is percent with no MTD
    mtd_sel_vec <- head(res$summary$mtd_selection_percent, -1)
    pct_no_mtd <- tail(res$summary$mtd_selection_percent, 1)

    mtd_values <- round(c(mtd_sel_vec, pct_no_mtd), 1)
    row_data <- c("", "MTD Selected (%)", as.character(mtd_values[seq_len(n_doses)]), as.character(mtd_values[n_doses + 1]))
    table_data <- rbind(table_data, row_data)

    # Row 3: Number of Participants Treated
    pts_values <- round(c(res$summary$avg_n_pts, res$summary$avg_total_n_pts), 1)
    row_data <- c("", "Number of Participants Treated", as.character(pts_values[seq_len(n_doses)]), as.character(pts_values[n_doses + 1]))
    table_data <- rbind(table_data, row_data)

    # Row 4: Number of Participants w/ DLTs
    tox_values <- round(c(res$summary$avg_n_tox, res$summary$avg_total_n_tox), 1)
    row_data <- c("", "Number of Participants w/ DLTs", as.character(tox_values[seq_len(n_doses)]), as.character(tox_values[n_doses + 1]))
    table_data <- rbind(table_data, row_data)
  }

  # Construct column names: dose labels + total column
  dose_col_names <- paste0("DL", seq_len(n_doses))
  total_col_name <- "Total/No MTD"
  colnames(table_data) <- c("Scenario", "Item", dose_col_names, total_col_name)

  cat("  Completed.\n\n")

  # ===== CONSTRUCT RETURN OBJECT =====
  return_obj <- list(
    results_by_scenario = results_by_scenario,
    combined_summary_df = table_data,
    scenario_names = scenario_names,
    n_doses = n_doses,
    call = call_expr
  )

  class(return_obj) <- c("boin_multi_summary", "list")

  # Print final summary
  cat("========================================\n")
  cat("Multi-Scenario Summary\n")
  cat("========================================\n\n")
  print(return_obj)

  cat("\n")

  # Return invisibly to avoid duplication when printed
  invisible(return_obj)
}
