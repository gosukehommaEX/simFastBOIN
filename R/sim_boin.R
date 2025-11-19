#' Run BOIN Simulation
#'
#' @description
#'   Execute multiple BOIN trial simulations and return aggregated operating
#'   characteristics.
#'
#' @param n_trials
#'   Numeric. Number of trials to simulate. Default is 10000.
#'
#' @param target
#'   Numeric. The target toxicity probability (e.g., 0.30 for 30%).
#'
#' @param p_true
#'   Numeric vector. True toxicity probabilities for each dose.
#'
#' @param n_doses
#'   Numeric. Number of doses evaluated.
#'
#' @param n_cohort
#'   Numeric. Maximum number of cohorts per trial.
#'
#' @param cohort_size
#'   Numeric vector or scalar specifying patients per cohort.
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
#'   \item{summary}{Aggregated summary statistics}
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

  # Initialize
  n_pts <- rep(0, n_doses)
  n_tox <- rep(0, n_doses)
  current_dose <- 1
  eliminated_doses <- rep(FALSE, n_doses)

  # Pre-compute cohort size vector
  if (length(cohort_size) == 1) {
    cohort_size <- rep(cohort_size, n_cohort)
  } else if (length(cohort_size) < n_cohort) {
    cohort_size <- c(cohort_size, rep(cohort_size[length(cohort_size)], n_cohort - length(cohort_size)))
  }
  cohort_size <- cohort_size[1:n_cohort]

  # Main trial loop
  for (cohort in seq_len(n_cohort)) {

    # Early stopping check
    if (n_pts[current_dose] >= n_earlystop) {
      break
    }

    # Get cohort size
    current_cohort_size <- cohort_size[cohort]

    # Generate DLT data
    dlt_count <- rbinom(1, current_cohort_size, p_true[current_dose])

    # Update counts
    n_pts[current_dose] <- n_pts[current_dose] + current_cohort_size
    n_tox[current_dose] <- n_tox[current_dose] + dlt_count

    # Safety stopping rule at lowest dose
    if (current_dose == 1 && n_pts[1] >= 3) {
      n_pts_for_stop_table <- min(n_pts[1], ncol(stopping_boundaries))
      decision_stop <- stopping_boundaries[n_tox[1] + 1, n_pts_for_stop_table]

      if (!is.na(decision_stop) && decision_stop == "STOP") {
        return(list(
          n_pts = n_pts,
          n_tox = n_tox,
          mtd = NA,
          iso_est = rep(NA, n_doses),
          reason = "lowest_dose_too_toxic",
          cohorts_completed = cohort
        ))
      }
    }

    # Get dose decision
    if (n_pts[current_dose] <= ncol(decision_table)) {
      decision <- decision_table[n_tox[current_dose] + 1, n_pts[current_dose]]
    } else {
      decision <- NA
    }

    # Handle NA decisions
    if (is.na(decision)) {
      decision <- "S"
    }

    # Dose adjustment and elimination processing
    # Only "DE" triggers elimination
    if (decision == "DE") {
      eliminated_doses[current_dose:n_doses] <- TRUE

      if (current_dose > 1) {
        current_dose <- current_dose - 1
        while (current_dose > 1 && eliminated_doses[current_dose]) {
          current_dose <- current_dose - 1
        }
      } else {
        # Lowest dose eliminated - stop trial with no MTD
        return(list(
          n_pts = n_pts,
          n_tox = n_tox,
          mtd = NA,
          iso_est = rep(NA, n_doses),
          reason = "lowest_dose_eliminated",
          cohorts_completed = cohort
        ))
      }
    }

    # Update dose
    if (decision == "E") {
      if (current_dose < n_doses && !eliminated_doses[current_dose + 1]) {
        current_dose <- current_dose + 1
      }
    } else if (decision == "D") {
      if (current_dose > 1) {
        current_dose <- current_dose - 1
        while (current_dose > 1 && eliminated_doses[current_dose]) {
          current_dose <- current_dose - 1
        }
      }
    }
  }

  # MTD Selection Phase

  # Step 1: Compute isotonically-adjusted toxicity rates
  iso_est <- isotonic_regression(n_pts, n_tox, min_sample = min_mtd_sample)

  # Step 2: Set eliminated doses to NA
  iso_est[eliminated_doses] <- NA

  # Step 3: Compute distance from target for each dose
  diffs <- abs(iso_est - target)

  # Step 4: Check if any valid dose remains for MTD selection
  if (all(is.na(diffs))) {
    return(list(
      n_pts = n_pts,
      n_tox = n_tox,
      mtd = NA,
      iso_est = iso_est,
      reason = "no_valid_dose",
      cohorts_completed = n_cohort
    ))
  }

  # Step 5: Identify candidate dose(s) closest to target
  mtd_candidates <- which(diffs == min(diffs, na.rm = TRUE))

  # Step 6: Tiebreaker for multiple candidates
  if (length(mtd_candidates) > 1) {
    candidate_estimates <- iso_est[mtd_candidates]

    above_target <- candidate_estimates > target
    below_target <- candidate_estimates < target

    if (all(above_target)) {
      mtd <- min(mtd_candidates)
    } else if (all(below_target)) {
      mtd <- max(mtd_candidates)
    } else {
      mtd <- max(mtd_candidates)
    }
  } else {
    mtd <- mtd_candidates[1]
  }

  return(list(
    n_pts = n_pts,
    n_tox = n_tox,
    mtd = mtd,
    iso_est = iso_est,
    reason = "trial_completed",
    cohorts_completed = n_cohort
  ))
}

#' PAVA Core Implementation
#'
#' @description
#'   Pool Adjacent Violators Algorithm for isotonic regression.
#'   Enforces non-decreasing monotonicity with weighted averaging.
#'
#' @param y Numeric vector. Adjusted toxicity rates.
#' @param w Numeric vector. Inverse variance weights.
#'
#' @return Numeric vector of isotonic-adjusted estimates.
#'
#' @keywords internal
.pava_core <- function(y, w) {

  n <- length(y)
  if (n == 0) return(numeric(0))
  if (n == 1) return(y)

  # Initialize: each element is its own level
  level_estimates <- y
  level_weights <- w
  level_size <- rep(1, n)

  # Iterate until monotonicity is achieved
  i <- 1
  while (i < length(level_estimates)) {
    # Check if adjacent levels violate monotonicity
    if (level_estimates[i] > level_estimates[i + 1]) {
      # Pool: merge levels i and i+1 with weighted average
      w_total <- level_weights[i] + level_weights[i + 1]
      level_estimates[i] <- (level_weights[i] * level_estimates[i] +
                               level_weights[i + 1] * level_estimates[i + 1]) / w_total
      level_weights[i] <- w_total
      level_size[i] <- level_size[i] + level_size[i + 1]

      # Remove the merged level i+1
      level_estimates <- level_estimates[-(i + 1)]
      level_weights <- level_weights[-(i + 1)]
      level_size <- level_size[-(i + 1)]

      # Backtrack: check previous level
      if (i > 1) i <- i - 1
    } else {
      i <- i + 1
    }
  }

  # Map pooled estimates back to original indices
  result <- numeric(n)
  idx_pointer <- 1
  for (k in 1:length(level_estimates)) {
    n_in_level <- level_size[k]
    result[idx_pointer:(idx_pointer + n_in_level - 1)] <- level_estimates[k]
    idx_pointer <- idx_pointer + n_in_level
  }

  return(result)
}

#' Isotonic Regression for Toxicity Rate Estimation
#'
#' @description
#'   Estimate toxicity rates at each dose level under monotonicity constraint.
#'   Uses PAVA (Pool Adjacent Violators Algorithm) without external packages.
#'   Incorporates pseudocount (Beta-Binomial prior).
#'
#' @param n_pts Numeric vector. Number of patients treated at each dose.
#' @param n_tox Numeric vector. Number of patients with DLTs at each dose.
#' @param min_sample Numeric. Minimum sample size for dose inclusion. Default is 6.
#'
#' @return Numeric vector of isotonic-adjusted toxicity rate estimates.
#'
#' @details
#'   Pseudocounts are added: (y + 0.05) / (n + 0.1).
#'   Patient counts weighted by inverse variance are used as weights in PAVA.
#'   Doses with fewer than min_sample patients return NA.
#'
#' @keywords internal
isotonic_regression <- function(n_pts, n_tox, min_sample = 6) {

  n_doses <- length(n_pts)
  iso_est <- rep(NA_real_, n_doses)

  # Consider only doses with >= min_sample patients
  valid_doses <- n_pts >= min_sample

  # Return all NAs if no dose has sufficient sample size
  if (sum(valid_doses) == 0) {
    return(iso_est)
  }

  # Initialize vectors for pseudocount-adjusted values
  tox_rate_adj <- rep(NA_real_, n_doses)
  variance_inv_weight <- rep(NA_real_, n_doses)

  # Compute pseudocount-adjusted toxicity rates and variance inverse weights
  for (i in which(valid_doses)) {
    # Adjusted toxicity rate: (y + 0.05) / (n + 0.1)
    tox_rate_adj[i] <- (n_tox[i] + 0.05) / (n_pts[i] + 0.1)

    # Inverse variance weight
    variance <- ((n_tox[i] + 0.05) * (n_pts[i] - n_tox[i] + 0.05)) /
      (((n_pts[i] + 0.1)^2) * (n_pts[i] + 0.1 + 1))

    variance_inv_weight[i] <- 1 / variance
  }

  # Apply PAVA to valid doses
  iso_est[valid_doses] <- .pava_core(
    tox_rate_adj[valid_doses],
    variance_inv_weight[valid_doses]
  )

  return(iso_est)
}

#' Summarize BOIN Simulation Results
#'
#' @keywords internal
.summarize_simulation <- function(simulation_results, n_doses) {

  n_trials <- length(simulation_results)

  # Initialize
  mtd_selected <- matrix(0, nrow = n_trials, ncol = n_doses)
  n_pts_all <- matrix(0, nrow = n_trials, ncol = n_doses)
  n_tox_all <- matrix(0, nrow = n_trials, ncol = n_doses)
  mtd_selected_flag <- rep(0, n_trials)

  # Aggregate individual trial results
  for (i in seq_len(n_trials)) {
    result <- simulation_results[[i]]

    n_pts_all[i, ] <- result$n_pts
    n_tox_all[i, ] <- result$n_tox

    # Flag trials where MTD was selected
    if (!is.na(result$mtd)) {
      mtd_selected[i, result$mtd] <- 1
      mtd_selected_flag[i] <- 1
    }
  }

  # Compute summary statistics
  summary <- list(
    mtd_selection_percent = colMeans(mtd_selected) * 100,
    avg_n_pts = colMeans(n_pts_all),
    avg_n_tox = colMeans(n_tox_all),
    percent_no_mtd = (1 - mean(mtd_selected_flag)) * 100,
    avg_total_n_pts = mean(rowSums(n_pts_all)),
    avg_total_n_tox = mean(rowSums(n_tox_all))
  )

  return(summary)
}
