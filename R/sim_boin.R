#' Run BOIN Simulation
#'
#' @description
#'   Execute multiple BOIN trial simulations and return aggregated operating
#'   characteristics. This is the main entry point for running simulation studies
#'   to evaluate BOIN design performance across different dose-toxicity scenarios.
#'
#' @details
#'   This function orchestrates the entire simulation workflow:
#'   1. Sets random seed for reproducibility
#'   2. Runs n_trials independent trial simulations using `.sim_boin_one_trial()`
#'   3. Aggregates results into comprehensive summary statistics
#'   4. Returns both detailed trial-level results and summary statistics
#'
#'   The detailed results enable custom post-hoc analyses if needed.
#'   Summary statistics include MTD selection rates, patient allocation patterns,
#'   and safety metrics essential for protocol development and regulatory submissions.
#'
#'   Progress is reported at 10 intervals during execution to monitor simulation status.
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
#' @param cutoff_eli
#'   Numeric. Cutoff probability for dose elimination. Default is 0.95.
#'
#' @param seed
#'   Numeric. Random seed for reproducibility. Default is 123.
#'
#' @return
#'   A list containing:
#'   \item{detailed_results}{List of results from each individual trial}
#'   \item{summary}{Aggregated summary statistics from `.summarize_simulation()`}
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic 3-dose study with 30% target DLT rate
#' target <- 0.30
#' p_true <- c(0.10, 0.25, 0.40)
#' boin_bound <- get_boin_boundary(target)
#' decision_table <- get_boin_decision(
#'   target = target,
#'   lambda_e = boin_bound$lambda_e,
#'   lambda_d = boin_bound$lambda_d,
#'   max_sample_size = 18,
#'   cutoff_eli = 0.95
#' )
#' stopping_boundaries <- get_boin_stopping_boundaries(
#'   target = target,
#'   max_sample_size = 18,
#'   cutoff_stop = 0.90
#' )
#'
#' result <- sim_boin(
#'   n_trials = 1000,
#'   target = target,
#'   p_true = p_true,
#'   n_doses = 3,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   decision_table = decision_table,
#'   stopping_boundaries = stopping_boundaries,
#'   min_mtd_sample = 6,
#'   seed = 123
#' )
#'
#' # Display summary statistics
#' print(result$summary)
#'
#' # Example 2: 9-dose study with varying cohort sizes
#' target <- 0.30
#' p_true <- seq(0.05, 0.45, by = 0.05)
#' boin_bound <- get_boin_boundary(target)
#' decision_table <- get_boin_decision(
#'   target = target,
#'   lambda_e = boin_bound$lambda_e,
#'   lambda_d = boin_bound$lambda_d,
#'   max_sample_size = 18,
#'   cutoff_eli = 0.95
#' )
#' stopping_boundaries <- get_boin_stopping_boundaries(
#'   target = target,
#'   max_sample_size = 18,
#'   cutoff_stop = 0.90
#' )
#'
#' result <- sim_boin(
#'   n_trials = 10000,
#'   target = target,
#'   p_true = p_true,
#'   n_doses = 9,
#'   n_cohort = 48,
#'   cohort_size = 3,
#'   decision_table = decision_table,
#'   stopping_boundaries = stopping_boundaries,
#'   cutoff_eli = 0.95,
#'   min_mtd_sample = 6,
#'   seed = 123
#' )
#'
#' # Access detailed results for custom analysis
#' first_trial_result <- result$detailed_results[[1]]
#' print(first_trial_result$n_pts)  # Patients per dose in first trial
#' print(first_trial_result$mtd)    # Selected MTD in first trial
#' }
#'
#' @export
sim_boin <- function(
    n_trials = 10000,
    target,
    p_true,
    n_doses,
    n_cohort,
    cohort_size,
    decision_table,
    stopping_boundaries,
    n_earlystop = 18,
    min_mtd_sample = 6,
    cutoff_eli = 0.95,
    seed = 123
) {

  set.seed(seed)

  cat("========================================\n")
  cat("Starting BOIN Simulation\n")
  cat("Number of trials:", n_trials, "\n")
  cat("Target DLT rate:", target * 100, "%\n")
  cat("========================================\n\n")

  # Initialize list to store results from all trials
  simulation_results <- vector("list", n_trials)

  # Determine progress reporting interval
  progress_interval <- max(1, n_trials %/% 10)

  # Main simulation loop
  for (trial in seq_len(n_trials)) {

    # Progress reporting
    if (trial %% progress_interval == 0) {
      cat("Progress: ", trial, " / ", n_trials, " trials completed\n", sep = "")
    }

    simulation_results[[trial]] <- .sim_boin_one_trial(
      target = target,
      p_true = p_true,
      n_doses = n_doses,
      n_cohort = n_cohort,
      cohort_size = cohort_size,
      decision_table = decision_table,
      stopping_boundaries = stopping_boundaries,
      n_earlystop = n_earlystop,
      min_mtd_sample = min_mtd_sample,
      cutoff_eli = cutoff_eli
    )
  }

  cat("\nSimulation completed!\n\n")

  # Compute summary statistics
  summary_result <- .summarize_simulation(simulation_results, n_doses)

  return(list(
    detailed_results = simulation_results,
    summary = summary_result
  ))
}

#' Simulate One BOIN Trial
#'
#' @keywords internal
.sim_boin_one_trial <- function(
    target,
    p_true,
    n_doses,
    n_cohort,
    cohort_size,
    decision_table,
    stopping_boundaries,
    n_earlystop,
    min_mtd_sample,
    cutoff_eli
) {

  # Initialize patient and toxicity counts for each dose
  n_pts <- rep(0L, n_doses)
  n_tox <- rep(0L, n_doses)

  # Start at dose level 1
  current_dose <- 1L

  # Initialize stopping flag
  trial_stopped <- FALSE
  reason <- "Trial completed normally"

  # Initialize dose elimination status
  dose_eliminated <- rep(FALSE, n_doses)

  # Track number of cohorts completed
  cohorts_completed <- 0L

  # Cohort-by-cohort enrollment
  for (cohort in seq_len(n_cohort)) {

    # Stop if trial is stopped or all doses eliminated
    if (trial_stopped || all(dose_eliminated)) {
      if (all(dose_eliminated)) {
        reason <- "All doses eliminated"
      }
      break
    }

    # Determine cohort size for this cohort
    if (length(cohort_size) > 1) {
      cs <- cohort_size[cohort]
    } else {
      cs <- cohort_size
    }

    # Enroll patients at current dose
    n_pts[current_dose] <- n_pts[current_dose] + cs

    # Generate DLT outcomes
    tox_outcome <- rbinom(1, cs, p_true[current_dose])
    n_tox[current_dose] <- n_tox[current_dose] + tox_outcome

    cohorts_completed <- cohorts_completed + 1L

    # --- Early Stopping Check (at lowest dose) ---
    if (current_dose == 1) {
      n_current <- n_pts[1]
      y_current <- n_tox[1]

      if (n_current >= n_earlystop) {
        decision <- stopping_boundaries[y_current + 1, n_current]
        if (decision == "stop") {
          trial_stopped <- TRUE
          reason <- "Early stop at lowest dose"
          break
        }
      }
    }

    # --- Dose Elimination Check ---
    n_current <- n_pts[current_dose]
    y_current <- n_tox[current_dose]
    decision <- decision_table[y_current + 1, n_current]

    if (decision == "elim") {
      dose_eliminated[current_dose] <- TRUE

      # Eliminate all higher doses
      if (current_dose < n_doses) {
        dose_eliminated[(current_dose + 1):n_doses] <- TRUE
      }

      # Find next lower available dose
      if (current_dose > 1) {
        for (d in (current_dose - 1):1) {
          if (!dose_eliminated[d]) {
            current_dose <- d
            break
          }
        }
      }

      # If no lower dose available, stop trial
      if (dose_eliminated[current_dose]) {
        trial_stopped <- TRUE
        reason <- "No available dose"
        break
      }

      next  # Skip dose assignment for this cohort
    }

    # --- Dose Assignment for Next Cohort ---
    if (cohort < n_cohort && !trial_stopped) {

      if (decision == "stay") {
        # Stay at current dose
        current_dose <- current_dose

      } else if (decision == "escalate") {
        # Escalate to next higher dose if available
        if (current_dose < n_doses) {
          next_dose <- current_dose + 1L
          # Find next non-eliminated dose
          while (next_dose <= n_doses && dose_eliminated[next_dose]) {
            next_dose <- next_dose + 1L
          }
          if (next_dose <= n_doses) {
            current_dose <- next_dose
          }
        }

      } else if (decision == "deescalate") {
        # De-escalate to next lower dose if available
        if (current_dose > 1) {
          next_dose <- current_dose - 1L
          # Find next non-eliminated dose
          while (next_dose >= 1 && dose_eliminated[next_dose]) {
            next_dose <- next_dose - 1L
          }
          if (next_dose >= 1) {
            current_dose <- next_dose
          }
        }
      }
    }
  }

  # --- MTD Selection Phase ---
  mtd <- NA_integer_
  iso_est <- rep(NA_real_, n_doses)

  # Only consider doses with sufficient sample size and not eliminated
  eligible_doses <- which(n_pts >= min_mtd_sample & !dose_eliminated)

  if (length(eligible_doses) > 0) {

    # Isotonic regression to estimate dose-toxicity curve
    p_est <- n_tox / n_pts
    iso_est <- .isoreg_estimate(p_est)

    # Select dose closest to target among eligible doses
    distances <- abs(iso_est[eligible_doses] - target)
    best_idx <- which.min(distances)
    mtd <- eligible_doses[best_idx]
  }

  return(list(
    n_pts = n_pts,
    n_tox = n_tox,
    mtd = mtd,
    iso_est = iso_est,
    reason = reason,
    cohorts_completed = cohorts_completed
  ))
}

#' Perform Isotonic Regression for Dose-Toxicity Estimation
#'
#' @keywords internal
.isoreg_estimate <- function(p_est) {

  n_doses <- length(p_est)

  # Handle edge cases
  if (all(is.na(p_est))) {
    return(rep(NA_real_, n_doses))
  }

  # Replace NA with 0 for doses not yet evaluated
  p_est[is.na(p_est)] <- 0

  # Use built-in isotonic regression
  iso_fit <- stats::isoreg(seq_len(n_doses), p_est)

  # Return fitted values
  return(iso_fit$yf)
}

#' Summarize BOIN Simulation Results
#'
#' @keywords internal
.summarize_simulation <- function(simulation_results, n_doses) {

  n_trials <- length(simulation_results)

  # Initialize summary matrices
  mtd_count <- rep(0L, n_doses)
  total_n_pts <- rep(0, n_doses)
  total_n_tox <- rep(0, n_doses)
  no_mtd_count <- 0L

  # Aggregate results from all trials
  for (trial_result in simulation_results) {

    # Count MTD selections
    if (!is.na(trial_result$mtd)) {
      mtd_count[trial_result$mtd] <- mtd_count[trial_result$mtd] + 1L
    } else {
      no_mtd_count <- no_mtd_count + 1L
    }

    # Sum patient counts
    total_n_pts <- total_n_pts + trial_result$n_pts
    total_n_tox <- total_n_tox + trial_result$n_tox
  }

  # Calculate summary statistics
  mtd_selection_percent <- (mtd_count / n_trials) * 100
  avg_n_pts <- total_n_pts / n_trials
  avg_n_tox <- total_n_tox / n_trials
  percent_no_mtd <- (no_mtd_count / n_trials) * 100

  return(list(
    mtd_selection_percent = mtd_selection_percent,
    avg_n_pts = avg_n_pts,
    avg_n_tox = avg_n_tox,
    percent_no_mtd = percent_no_mtd,
    avg_total_n_pts = sum(avg_n_pts),
    avg_total_n_tox = sum(avg_n_tox)
  ))
}
