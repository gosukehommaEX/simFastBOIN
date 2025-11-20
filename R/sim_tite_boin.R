#' Run TITE-BOIN Simulation
#'
#' @description
#' Execute multiple TITE-BOIN (Time-to-Event Bayesian Optimal Interval) trial
#' simulations and return aggregated operating characteristics. This function
#' implements the approximated likelihood method from Lin and Yuan (2020),
#' which uses Effective Sample Size (ESS) for decision making.
#'
#' @details
#' This function orchestrates the entire TITE-BOIN simulation workflow:
#' \enumerate{
#'   \item Sets random seed for reproducibility
#'   \item Runs n_trials independent trial simulations
#'   \item Aggregates results into comprehensive summary statistics
#'   \item Returns both detailed trial-level results and summary statistics
#' }
#'
#' The TITE-BOIN design uses ESS (Effective Sample Size) to account for pending
#' patients with incomplete follow-up. Decision making is based on comparing the
#' ESS-adjusted toxicity rate with BOIN boundaries (lambda_e and lambda_d), without
#' requiring pre-computed decision tables.
#'
#' Patient arrival follows an exponential distribution with rate parameter accrual_rate,
#' allowing realistic modeling of enrollment variability.
#'
#' Progress is reported at 10 intervals during execution to monitor simulation status.
#'
#' @param n_trials Numeric. Number of trials to simulate. Default is 1000.
#'   Typically 1000-10000 for detailed operating characteristics.
#' @param target Numeric. The target toxicity probability (e.g., 0.30 for 30%).
#' @param p_true Numeric vector. True toxicity probabilities for each dose.
#'   Length determines number of doses evaluated.
#' @param n_doses Numeric. Number of doses evaluated.
#' @param n_cohort Numeric. Maximum number of cohorts per trial.
#' @param cohort_size Numeric. Number of patients per cohort (scalar only for TITE-BOIN).
#' @param lambda_e Numeric. Escalation boundary from `get_boin_boundary()`.
#' @param lambda_d Numeric. De-escalation boundary from `get_boin_boundary()`.
#' @param stopping_boundaries Character matrix. Trial stopping rule table from
#'   `get_boin_stopping_boundaries()`.
#' @param accrual_rate Numeric. Patient accrual rate (patients per unit time).
#'   For example, accrual_rate = 2 means on average 2 patients per time unit.
#'   Inter-arrival times follow exponential distribution with this rate.
#'   Default is 1.
#' @param assessment_window Numeric. Length of the DLT assessment window
#'   (e.g., 28 days). Default is 28.
#' @param maxpen Numeric. Maximum allowable ratio of pending patients at current dose.
#'   If (n_pending / n_current) > maxpen, suspend enrollment until more patients complete.
#'   Default is 0.5 (suspend if >50% pending).
#' @param alpha1 Numeric. Late-onset parameter (0 < alpha1 < 1). Represents the
#'   proportion of toxicity that occurs in the last fraction (alpha2) of the
#'   assessment window. Default is 0.5.
#' @param alpha2 Numeric. Weibull shape parameter (0 < alpha2 < 1). Defines the
#'   last fraction of the assessment window where alpha1 proportion of toxicity
#'   occurs. Default is 0.5.
#' @param distribution Character. Distribution for time-to-toxicity.
#'   One of "weibull" (default), "log-logistic", or "uniform".
#' @param n_earlystop Numeric. Sample size triggering early stopping check. Default is 18.
#' @param min_mtd_sample Numeric. Minimum sample size for MTD consideration. Default is 6.
#' @param cutoff_eli Numeric. Cutoff probability for dose elimination. Default is 0.95.
#' @param seed Numeric. Random seed for reproducibility. Default is 123.
#'
#' @return A list containing:
#'   \item{detailed_results}{List of results from each individual trial}
#'   \item{summary}{Aggregated summary statistics including:
#'     \itemize{
#'       \item mtd_selection_percent: Percentage of trials selecting each dose as MTD
#'       \item avg_n_pts: Average number of patients treated at each dose
#'       \item avg_n_tox: Average number of DLTs at each dose
#'       \item percent_no_mtd: Percentage of trials where no MTD was selected
#'       \item avg_total_n_pts: Average total patients across all doses
#'       \item avg_total_n_tox: Average total DLTs across all doses
#'       \item avg_trial_duration: Average trial duration in time units
#'     }
#'   }
#'
#' @references
#' Lin, R., & Yuan, Y. (2020). Time-to-event model-assisted designs for
#' dose-finding trials with delayed toxicity. \emph{Biostatistics}, 21(4), 807-824.
#'
#' Yuan, Y., Lin, R., Li, D., Nie, L., & Warren, K. E. (2018).
#' Time-to-Event Bayesian Optimal Interval Design to Accelerate Phase I Trials.
#' \emph{Clinical Cancer Research}, 24(20), 4921-4930.
#'
#' @importFrom stats rbinom pbeta rexp
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic 3-dose TITE-BOIN study with 30% target DLT rate
#' target <- 0.30
#' p_true <- c(0.10, 0.25, 0.40)
#' boin_bound <- get_boin_boundary(target)
#' stopping_boundaries <- get_boin_stopping_boundaries(
#'   target = target,
#'   max_sample_size = 18,
#'   cutoff_stop = 0.90
#' )
#'
#' result <- sim_tite_boin(
#'   n_trials = 1000,
#'   target = target,
#'   p_true = p_true,
#'   n_doses = 3,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   lambda_e = boin_bound$lambda_e,
#'   lambda_d = boin_bound$lambda_d,
#'   stopping_boundaries = stopping_boundaries,
#'   accrual_rate = 2,
#'   assessment_window = 28,
#'   maxpen = 0.5,
#'   seed = 123
#' )
#'
#' # Display summary statistics
#' print(result$summary)
#' }
#'
#' @export
sim_tite_boin <- function(
    n_trials = 1000,
    target,
    p_true,
    n_doses,
    n_cohort,
    cohort_size,
    lambda_e,
    lambda_d,
    stopping_boundaries,
    accrual_rate = 1,
    assessment_window = 28,
    maxpen = 0.5,
    alpha1 = 0.5,
    alpha2 = 0.5,
    distribution = "weibull",
    n_earlystop = 18,
    min_mtd_sample = 6,
    cutoff_eli = 0.95,
    seed = 123
) {

  set.seed(seed)

  cat("========================================\n")
  cat("Starting TITE-BOIN Simulation\n")
  cat("Number of trials:", n_trials, "\n")
  cat("Target DLT rate:", target * 100, "%\n")
  cat("Number of doses:", n_doses, "\n")
  cat("Accrual rate:", accrual_rate, "patients per time unit\n")
  cat("Max pending ratio:", maxpen, "\n")
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

    simulation_results[[trial]] <- .sim_tite_boin_one_trial(
      target = target,
      p_true = p_true,
      n_doses = n_doses,
      n_cohort = n_cohort,
      cohort_size = cohort_size,
      lambda_e = lambda_e,
      lambda_d = lambda_d,
      stopping_boundaries = stopping_boundaries,
      accrual_rate = accrual_rate,
      assessment_window = assessment_window,
      maxpen = maxpen,
      alpha1 = alpha1,
      alpha2 = alpha2,
      distribution = distribution,
      n_earlystop = n_earlystop,
      min_mtd_sample = min_mtd_sample,
      cutoff_eli = cutoff_eli
    )
  }

  cat("\nSimulation completed!\n\n")

  # Compute summary statistics
  summary_result <- .summarize_simulation_tite_boin(simulation_results, n_doses)

  return(list(
    detailed_results = simulation_results,
    summary = summary_result
  ))
}

#' Simulate One TITE-BOIN Trial
#'
#' @description
#' Conduct a single TITE-BOIN trial simulation with time-to-event data generation
#' and ESS-based decision making following Lin and Yuan (2020).
#'
#' @keywords internal
.sim_tite_boin_one_trial <- function(
    target,
    p_true,
    n_doses,
    n_cohort,
    cohort_size,
    lambda_e,
    lambda_d,
    stopping_boundaries,
    accrual_rate,
    assessment_window,
    maxpen,
    alpha1,
    alpha2,
    distribution,
    n_earlystop,
    min_mtd_sample,
    cutoff_eli
) {

  # Initialize tracking variables
  n_pts <- rep(0, n_doses)  # Number of patients enrolled at each dose
  n_tox <- rep(0, n_doses)  # Number of observed DLTs at each dose
  current_dose <- 1
  eliminated_doses <- rep(FALSE, n_doses)

  # Patient-level tracking
  patient_data <- list(
    dose = integer(0),
    enroll_time = numeric(0),
    t_tox = numeric(0),
    tox = integer(0)
  )

  current_time <- 0

  # Main trial loop: enroll cohorts sequentially
  for (cohort in seq_len(n_cohort)) {

    # Early stopping check: if current dose reaches n_earlystop, stop trial
    if (n_pts[current_dose] >= n_earlystop) {
      break
    }

    # Check if current dose is eliminated
    if (eliminated_doses[current_dose]) {
      break
    }

    # Check suspension rule: if too many pending patients, wait
    pts_at_current <- which(patient_data$dose == current_dose)
    if (length(pts_at_current) > 0) {
      completed_mask <- (patient_data$enroll_time[pts_at_current] + assessment_window) <= current_time
      n_completed <- sum(completed_mask)
      n_pending <- sum(!completed_mask)
      n_current <- length(pts_at_current)

      # Suspension rule: if pending ratio > maxpen, wait
      if (n_pending > n_current * maxpen && n_current > 0) {
        # Wait until enough patients complete assessment
        # Find the earliest time when pending ratio drops below maxpen
        pending_complete_times <- patient_data$enroll_time[pts_at_current][!completed_mask] + assessment_window

        if (length(pending_complete_times) > 0) {
          # Sort completion times
          sorted_times <- sort(pending_complete_times)

          # Find the time when pending ratio drops to maxpen
          for (t_check in sorted_times) {
            completed_at_t <- (patient_data$enroll_time[pts_at_current] + assessment_window) <= t_check
            n_completed_at_t <- sum(completed_at_t)
            n_pending_at_t <- sum(!completed_at_t)

            if (n_pending_at_t <= n_current * maxpen) {
              current_time <- t_check
              break
            }
          }
        }
      }
    }

    # Enroll cohort at current dose with exponential inter-arrival times
    for (pt in seq_len(cohort_size)) {

      # Generate inter-arrival time from exponential distribution
      if (length(patient_data$dose) == 0 && pt == 1) {
        # First patient starts at time 0
        inter_arrival <- 0
      } else if (pt == 1) {
        # First patient of new cohort: decision time
        inter_arrival <- 0
      } else {
        # Within cohort: exponential inter-arrival times
        inter_arrival <- rexp(1, rate = accrual_rate)
      }

      current_time <- current_time + inter_arrival

      # Generate time-to-toxicity data for this patient
      tite_data <- rtite(
        n = 1,
        prob = p_true[current_dose],
        alpha1 = alpha1,
        alpha2 = alpha2,
        distribution = distribution,
        maxt = assessment_window
      )

      # Record patient data
      patient_data$dose <- c(patient_data$dose, current_dose)
      patient_data$enroll_time <- c(patient_data$enroll_time, current_time)
      patient_data$t_tox <- c(patient_data$t_tox, tite_data$t_tox)
      patient_data$tox <- c(patient_data$tox, tite_data$tox)

      # Update counts
      n_pts[current_dose] <- n_pts[current_dose] + 1
    }

    # After enrolling cohort, advance time to decision point
    # Add exponential inter-arrival time for next cohort
    current_time <- current_time + rexp(1, rate = accrual_rate)

    # Make dose decision
    # Calculate ESS for current dose
    pts_at_current <- which(patient_data$dose == current_dose)

    # Determine completed and pending patients
    completed_mask <- (patient_data$enroll_time[pts_at_current] + assessment_window) <= current_time
    n_completed <- sum(completed_mask)
    n_tox_completed <- sum(patient_data$tox[pts_at_current][completed_mask])

    # Pending patients: enrolled but not yet completed assessment
    pending_mask <- !completed_mask
    follow_up_times <- numeric(0)
    if (sum(pending_mask) > 0) {
      follow_up_times <- pmin(
        current_time - patient_data$enroll_time[pts_at_current][pending_mask],
        assessment_window
      )
    }

    # Calculate ESS using the calculate_ess function
    ess <- calculate_ess(n_completed, follow_up_times, assessment_window)

    # Update observed toxicities (only count completed patients)
    n_tox[current_dose] <- n_tox_completed

    # Make dose escalation decision using ESS and BOIN boundaries
    next_dose <- current_dose  # Default: stay

    if (ess >= 1) {
      # Calculate adjusted toxicity rate (p_tilde)
      p_tilde <- n_tox_completed / ess

      # Direct comparison with BOIN boundaries (no decision table needed)
      if (p_tilde <= lambda_e) {
        # Toxicity rate is low enough to escalate
        next_dose <- current_dose + 1
      } else if (p_tilde >= lambda_d) {
        # Toxicity rate is too high, de-escalate
        next_dose <- current_dose - 1
      } else {
        # Toxicity rate is in acceptable range, stay
        next_dose <- current_dose
      }
    }

    # Dose elimination check (only if enough data)
    # CRITICAL: Use total enrolled patients (n_pts), not just completed patients
    # This matches TITEgBOIN implementation
    n_current <- n_pts[current_dose]

    if (n_current >= 3) {
      # Check using stopping boundaries
      # Note: stopping_boundaries uses TOTAL enrolled patients, not completed
      if (n_current <= nrow(stopping_boundaries) &&
          n_tox_completed + 1 <= ncol(stopping_boundaries)) {
        stop_decision <- stopping_boundaries[n_tox_completed + 1, n_current]

        if (!is.na(stop_decision) && stop_decision == "STOP") {
          # Eliminate current and all higher doses
          eliminated_doses[current_dose:n_doses] <- TRUE

          if (current_dose == 1) {
            break  # Stop trial if dose 1 is eliminated
          }
          next_dose <- current_dose - 1
        }
      }

      # Additional elimination check using posterior probability
      # Use total enrolled patients (n_current) for consistency with TITEgBOIN
      post_prob <- 1 - pbeta(target, n_tox_completed + 1, n_current - n_tox_completed + 1)
      if (post_prob > cutoff_eli) {
        eliminated_doses[current_dose:n_doses] <- TRUE
        if (current_dose == 1) {
          break
        }
        next_dose <- current_dose - 1
      }
    }

    # Update current dose for next cohort
    if (next_dose > current_dose) {
      # Escalate
      if (next_dose <= n_doses && !eliminated_doses[next_dose]) {
        current_dose <- next_dose
      }
    } else if (next_dose < current_dose) {
      # De-escalate
      if (next_dose >= 1 && !eliminated_doses[next_dose]) {
        current_dose <- next_dose
      }
    }

    # Check if all doses are eliminated
    if (all(eliminated_doses)) {
      break
    }
  }

  # Final toxicity counts (all patients who completed assessment by end of trial)
  for (dose_idx in 1:n_doses) {
    pts_at_dose <- which(patient_data$dose == dose_idx)
    if (length(pts_at_dose) > 0) {
      completed_at_end <- (patient_data$enroll_time[pts_at_dose] + assessment_window) <= current_time
      n_tox[dose_idx] <- sum(patient_data$tox[pts_at_dose][completed_at_end])
    }
  }

  # MTD Selection using isotonic regression (existing function)
  if (eliminated_doses[1] || all(n_pts == 0)) {
    mtd <- NA_integer_
    iso_est <- rep(NA_real_, n_doses)
  } else {
    # Use existing isotonic_regression function from the package
    iso_est <- isotonic_regression(n_pts, n_tox, min_sample = min_mtd_sample)

    # Set eliminated doses to NA
    iso_est[eliminated_doses] <- NA_real_

    # Find dose closest to target
    diffs <- abs(iso_est - target)

    if (all(is.na(diffs))) {
      mtd <- NA_integer_
    } else {
      mtd_candidates <- which(diffs == min(diffs, na.rm = TRUE))

      # Tiebreaker for multiple candidates
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
    }
  }

  return(list(
    n_pts = n_pts,
    n_tox = n_tox,
    mtd = mtd,
    iso_est = iso_est,
    trial_duration = current_time
  ))
}

#' Summarize TITE-BOIN Simulation Results
#'
#' @description
#' Aggregate results from multiple TITE-BOIN trial simulations into comprehensive
#' summary statistics. Uses vectorized operations for efficient computation.
#'
#' @keywords internal
.summarize_simulation_tite_boin <- function(simulation_results, n_doses) {

  n_trials <- length(simulation_results)

  # Initialize
  mtd_selected <- matrix(0, nrow = n_trials, ncol = n_doses)
  n_pts_all <- matrix(0, nrow = n_trials, ncol = n_doses)
  n_tox_all <- matrix(0, nrow = n_trials, ncol = n_doses)
  duration_vector <- numeric(n_trials)
  mtd_selected_flag <- rep(0, n_trials)

  # Aggregate individual trial results
  for (i in seq_len(n_trials)) {
    result <- simulation_results[[i]]

    n_pts_all[i, ] <- result$n_pts
    n_tox_all[i, ] <- result$n_tox
    duration_vector[i] <- result$trial_duration

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
    avg_total_n_tox = mean(rowSums(n_tox_all)),
    avg_trial_duration = mean(duration_vector),
    sd_trial_duration = sd(duration_vector)
  )

  return(summary)
}
