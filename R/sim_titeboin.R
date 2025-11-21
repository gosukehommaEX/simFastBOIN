#' Simulate One TITE-BOIN Trial
#'
#' @description
#' Internal function to simulate a single Time-to-Event Bayesian Optimal Interval
#' (TITE-BOIN) dose-finding trial. This function is called repeatedly by
#' \code{sim_titeboin()} to generate operating characteristics.
#'
#' @param decision_table Decision table in simulation format from
#'   \code{convert_titeboin_decision_to_sim_format()}
#' @param stopping_boundaries Early stopping rule table from
#'   \code{get_boin_stopping_boundaries()}
#' @param target Target toxicity rate
#' @param p_true Vector of true toxicity probabilities for each dose
#' @param n_doses Number of dose levels
#' @param ncohort Number of cohorts to enroll
#' @param cohort_size Number of patients per cohort
#' @param maxt Maximum observation period for toxicity assessment
#' @param accrual_rate Patient accrual rate (patients per time unit)
#' @param alpha1 Late-onset toxicity parameter
#' @param alpha2 Shape parameter for time-to-event distribution
#' @param tite_distribution Distribution type ("weibull" or "log-logistic")
#' @param neli Minimum sample size for dose elimination
#' @param cutoff_eli Posterior probability cutoff for elimination
#' @param n_earlystop Maximum sample size for early stopping assessment
#' @param min_sample Minimum sample size for MTD selection
#'
#' @return List containing:
#' \describe{
#'   \item{n_patients_by_dose}{Number of patients treated at each dose}
#'   \item{n_tox_by_dose}{Number of toxicities observed at each dose}
#'   \item{mtd}{Selected MTD (NA if no MTD selected)}
#'   \item{trial_duration}{Total trial duration}
#' }
#'
#' @keywords internal
.sim_titeboin_one_trial <- function(
    decision_table, stopping_boundaries, target, p_true, n_doses, ncohort, cohort_size,
    maxt, accrual_rate, alpha1, alpha2, tite_distribution = "weibull",
    neli, cutoff_eli, n_earlystop, min_sample
) {

  t_enter <- numeric()
  t_event <- numeric()
  tox_indicator <- integer()
  dose_assigned <- integer()

  current_dose <- 1
  eliminated <- rep(FALSE, n_doses)
  current_time <- 0
  trial_stop <- FALSE
  decision_time <- 0

  for (cohort in 1:ncohort) {

    if (trial_stop || eliminated[current_dose]) break

    if (sum(dose_assigned == current_dose) >= n_earlystop) {
      break
    }

    cohort_data <- rtite(
      n = cohort_size,
      prob = p_true[current_dose],
      alpha1 = alpha1,
      alpha2 = alpha2,
      distribution = tite_distribution,
      maxt = maxt
    )

    for (i in 1:cohort_size) {
      current_time <- current_time + rexp(1, rate = accrual_rate)
      t_enter <- c(t_enter, current_time)
    }

    t_event <- c(t_event, cohort_data$t_tox)
    tox_indicator <- c(tox_indicator, cohort_data$tox)
    dose_assigned <- c(dose_assigned, rep(current_dose, cohort_size))

    decision_time <- current_time + maxt

    at_current <- dose_assigned == current_dose
    if (sum(at_current) == 0) next

    n_curr <- sum(at_current)
    curr_arrive <- t_enter[at_current]
    curr_event <- t_event[at_current]
    curr_tox <- tox_indicator[at_current]

    completion_times <- curr_arrive + curr_event
    completed_mask <- completion_times <= decision_time
    n_completed <- sum(completed_mask)
    ntox_completed <- sum(curr_tox[completed_mask])
    n_pending <- n_curr - n_completed

    mf <- n_pending / n_curr

    if (n_pending == 0) {
      stft <- 0
    } else {
      pending_mask <- !completed_mask
      pending_followup <- pmin(decision_time - curr_arrive[pending_mask], maxt)
      stft <- sum(pending_followup) / maxt
    }

    if (current_dose == 1 && n_completed >= 3) {
      if (ntox_completed <= nrow(stopping_boundaries) - 1 && n_completed <= ncol(stopping_boundaries)) {
        stop_decision <- stopping_boundaries[ntox_completed + 1, n_completed]
        if (!is.na(stop_decision) && stop_decision == "STOP") {
          trial_stop <- TRUE
          break
        }
      }
    }

    if (n_completed >= neli) {
      prob_exceed <- 1 - pbeta(target, ntox_completed + 1, n_completed - ntox_completed + 1)
      if (prob_exceed > cutoff_eli) {
        eliminated[current_dose:n_doses] <- TRUE
        if (current_dose == 1) {
          trial_stop <- TRUE
          break
        } else {
          current_dose <- current_dose - 1
          next
        }
      }
    }

    decision_match <- which(
      decision_table$n_patients == n_curr &
        decision_table$n_tox == ntox_completed &
        decision_table$n_pending == n_pending
    )

    if (length(decision_match) > 0) {
      idx <- decision_match[1]

      suspend_cond <- decision_table$suspend[idx]
      escalate_cond <- decision_table$escalate[idx]
      de_escalate_cond <- decision_table$de_escalate[idx]
      stay_cond <- decision_table$stay[idx]

      mf_thresh <- decision_table$mf_threshold[idx]
      stft_thresh <- decision_table$stft_threshold[idx]

      if (suspend_cond == "Yes" ||
          (!is.na(mf_thresh) && grepl("MF<", suspend_cond) && mf < mf_thresh) ||
          (!is.na(stft_thresh) && grepl("STFT>=", suspend_cond) && stft >= stft_thresh)) {
        next_dose <- current_dose
      }
      else if (escalate_cond == "Yes" ||
               (!is.na(mf_thresh) && grepl("MF>=", escalate_cond) && mf >= mf_thresh) ||
               (!is.na(stft_thresh) && grepl("STFT>=", escalate_cond) && stft >= stft_thresh)) {
        next_dose <- current_dose + 1
      }
      else if (de_escalate_cond == "Yes" || grepl("Yes & Eliminate", de_escalate_cond) ||
               (!is.na(stft_thresh) && grepl("STFT<=", de_escalate_cond) && stft <= stft_thresh)) {
        next_dose <- current_dose - 1
      }
      else if (stay_cond == "Yes" ||
               (!is.na(stft_thresh) && grepl("STFT<", stay_cond) && stft < stft_thresh) ||
               (!is.na(stft_thresh) && grepl("STFT>", stay_cond) && stft > stft_thresh)) {
        next_dose <- current_dose
      }
      else {
        next_dose <- current_dose
      }

      if (next_dose > current_dose) {
        if (next_dose <= n_doses && !eliminated[next_dose]) {
          current_dose <- next_dose
        }
      } else if (next_dose < current_dose) {
        if (next_dose >= 1) {
          current_dose <- next_dose
        }
      }
    }
  }

  n_patients_by_dose <- rep(0, n_doses)
  n_tox_by_dose <- rep(0, n_doses)

  # Calculate trial duration using decision_time (last decision time point)
  # This matches the RShiny app implementation where trial duration is
  # measured to the last decision time, not just the last patient entry
  trial_duration <- if (length(t_enter) > 0) decision_time else 0

  final_completion <- (t_enter + t_event) <= trial_duration

  for (d in 1:n_doses) {
    dmask <- dose_assigned == d
    if (sum(dmask) > 0) {
      n_patients_by_dose[d] <- sum(dmask)
      n_tox_by_dose[d] <- sum(tox_indicator[dmask & final_completion])
    }
  }

  if (trial_stop || eliminated[1]) {
    mtd <- NA_integer_
  } else {
    iso_tox <- isotonic_regression(n_patients_by_dose, n_tox_by_dose, min_sample = min_sample)
    iso_tox[eliminated] <- NA

    valid_doses <- !is.na(iso_tox)
    if (sum(valid_doses) == 0) {
      mtd <- NA_integer_
    } else {
      diff_from_target <- abs(iso_tox[valid_doses] - target)
      mtd_candidates <- which(diff_from_target == min(diff_from_target))

      if (length(mtd_candidates) > 1) {
        candidate_estimates <- iso_tox[valid_doses][mtd_candidates]
        above_target <- candidate_estimates > target
        below_target <- candidate_estimates < target

        if (all(above_target)) {
          mtd_idx <- min(mtd_candidates)
        } else if (all(below_target)) {
          mtd_idx <- max(mtd_candidates)
        } else {
          mtd_idx <- max(mtd_candidates)
        }
      } else {
        mtd_idx <- mtd_candidates[1]
      }

      mtd <- which(valid_doses)[mtd_idx]
    }
  }

  return(list(
    n_patients_by_dose = n_patients_by_dose,
    n_tox_by_dose = n_tox_by_dose,
    mtd = mtd,
    trial_duration = trial_duration
  ))
}


#' Summarize TITE-BOIN Simulation Results
#'
#' @description
#' Internal function to calculate summary statistics from TITE-BOIN simulation results.
#'
#' @param mtd_count Vector of MTD selection counts for each dose
#' @param n_patients_all Matrix of patient counts (trials x doses)
#' @param n_tox_all Matrix of toxicity counts (trials x doses)
#' @param durations Vector of trial durations
#' @param n_trials Number of trials
#'
#' @return List of summary statistics
#'
#' @keywords internal
.summarize_simulation_titeboin <- function(mtd_count, n_patients_all, n_tox_all,
                                           durations, n_trials) {

  no_mtd_count <- n_trials - sum(mtd_count)

  summary <- list(
    mtd_selection_percent = mtd_count / n_trials * 100,
    avg_n_pts = colMeans(n_patients_all),
    avg_total_n_pts = mean(rowSums(n_patients_all)),
    avg_n_tox = colMeans(n_tox_all),
    avg_total_n_tox = mean(rowSums(n_tox_all)),
    percent_no_mtd = no_mtd_count / n_trials * 100,
    avg_duration = mean(durations),
    sd_duration = sd(durations)
  )

  return(summary)
}


#' Simulate TITE-BOIN Trials
#'
#' @description
#' Conduct Monte Carlo simulations of Time-to-Event Bayesian Optimal Interval
#' (TITE-BOIN) dose-finding trials to evaluate operating characteristics.
#' The function simulates multiple trials under specified true toxicity scenarios
#' and summarizes MTD selection rates, patient allocation, and trial duration.
#'
#' @param n_trials Number of simulation trials
#' @param decision_table Decision table in simulation format from
#'   \code{convert_titeboin_decision_to_sim_format()}
#' @param stopping_boundaries Early stopping rule table from
#'   \code{get_boin_stopping_boundaries()}
#' @param target Target toxicity rate (e.g., 0.30 for 30\%)
#' @param p_true Vector of true toxicity probabilities for each dose level
#' @param n_doses Number of dose levels
#' @param ncohort Maximum number of cohorts to enroll
#' @param cohort_size Number of patients per cohort (typically 3)
#' @param maxt Maximum observation period for toxicity assessment (in time units)
#' @param accrual_rate Patient accrual rate (patients per time unit)
#' @param alpha1 Late-onset toxicity parameter (default: 0.5). This parameter
#'   controls the proportion of toxicities occurring in the first half of the
#'   observation period
#' @param alpha2 Shape parameter for the time-to-event distribution (default: 0.5).
#'   Controls the hazard pattern over time
#' @param tite_distribution Distribution for time-to-toxicity. Options:
#'   \itemize{
#'     \item "weibull": Weibull distribution (default)
#'     \item "log-logistic": Log-logistic distribution
#'   }
#' @param neli Minimum sample size required for dose elimination (default: 3)
#' @param cutoff_eli Posterior probability cutoff for dose elimination (default: 0.95)
#' @param n_earlystop Maximum sample size at a dose before forcing advancement
#'   (default: 100, effectively no limit)
#' @param min_sample Minimum sample size at a dose for MTD selection (default: 1)
#' @param seed Random seed for reproducibility (optional)
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{A list of summary statistics including:
#'     \itemize{
#'       \item \code{mtd_selection_percent}: Percentage of trials selecting each dose as MTD
#'       \item \code{avg_n_pts}: Average number of patients treated at each dose
#'       \item \code{avg_total_n_pts}: Average total number of patients per trial
#'       \item \code{avg_n_tox}: Average number of toxicities at each dose
#'       \item \code{avg_total_n_tox}: Average total toxicities per trial
#'       \item \code{percent_no_mtd}: Percentage of trials with no MTD selected
#'       \item \code{avg_duration}: Average trial duration
#'       \item \code{sd_duration}: Standard deviation of trial duration
#'     }
#'   }
#' }
#'
#' @details
#' The TITE-BOIN design extends the standard BOIN design to handle late-onset
#' toxicities by making real-time dose decisions while some patients' toxicity
#' outcomes are still pending. This accelerates trial conduct compared to
#' traditional designs that wait for all patients to complete follow-up before
#' enrolling the next cohort.
#'
#' The simulation assumes:
#' \itemize{
#'   \item Continuous patient accrual at the specified rate
#'   \item Toxicity assessment window of length \code{maxt}
#'   \item Dose decisions made after each cohort enrollment
#'   \item Isotonic regression for MTD selection at trial end
#' }
#'
#' @references
#' Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K.E. (2018).
#' Time-to-Event Bayesian Optimal Interval Design to Accelerate Phase I Trials.
#' Clinical Cancer Research, 24(20): 4921-4930.
#'
#' Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I
#' Clinical Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @importFrom stats pbeta rexp
#'
#' @examples
#' \dontrun{
#' # Set up trial parameters
#' target <- 0.30
#' p_true <- c(0.25, 0.41, 0.45, 0.49, 0.53)
#' n_doses <- length(p_true)
#'
#' # Generate BOIN boundaries
#' boundary <- get_boin_boundary(target)
#'
#' # Generate TITE-BOIN decision table
#' dt_full <- get_titeboin_decision(
#'   lambda_e = boundary$lambda_e,
#'   lambda_d = boundary$lambda_d,
#'   target = target,
#'   n_max = 18,
#'   cohort_size = 3
#' )
#'
#' # Convert to simulation format
#' decision_table <- convert_titeboin_decision_to_sim_format(dt_full)
#'
#' # Generate stopping boundaries
#' stopping_boundaries <- get_boin_stopping_boundaries(
#'   target = target,
#'   max_sample_size = 18,
#'   cutoff_stop = 0.90
#' )
#'
#' # Run simulation
#' result <- sim_titeboin(
#'   n_trials = 1000,
#'   decision_table = decision_table,
#'   stopping_boundaries = stopping_boundaries,
#'   target = target,
#'   p_true = p_true,
#'   n_doses = n_doses,
#'   ncohort = 10,
#'   cohort_size = 3,
#'   maxt = 3,
#'   accrual_rate = 1,
#'   alpha1 = 0.5,
#'   alpha2 = 0.5,
#'   tite_distribution = "weibull",
#'   neli = 3,
#'   cutoff_eli = 0.95,
#'   n_earlystop = 18,
#'   min_sample = 1,
#'   seed = 123
#' )
#'
#' # Display results
#' print_simulation_summary_titeboin(result$summary, n_doses)
#' }
#'
#' @seealso
#' \code{\link{get_titeboin_decision}} for generating decision tables,
#' \code{\link{convert_titeboin_decision_to_sim_format}} for format conversion,
#' \code{\link{print_simulation_summary_titeboin}} for displaying results
#'
#' @export
sim_titeboin <- function(
    n_trials, decision_table, stopping_boundaries, target, p_true, n_doses,
    ncohort, cohort_size, maxt, accrual_rate,
    alpha1 = 0.5, alpha2 = 0.5, tite_distribution = "weibull",
    neli = 3, cutoff_eli = 0.95, n_earlystop = 100, min_sample = 1, seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  cat("========================================\n")
  cat("Starting TITE-BOIN Simulation\n")
  cat("Number of trials:", n_trials, "\n")
  cat("========================================\n\n")

  mtd_count <- rep(0, n_doses)
  n_patients_all <- matrix(0, nrow = n_trials, ncol = n_doses)
  n_tox_all <- matrix(0, nrow = n_trials, ncol = n_doses)
  durations <- numeric(n_trials)

  for (trial in 1:n_trials) {

    if (trial %% 100 == 0) {
      cat("Progress:", trial, "/", n_trials, "\n")
    }

    result <- .sim_titeboin_one_trial(
      decision_table = decision_table,
      stopping_boundaries = stopping_boundaries,
      target = target,
      p_true = p_true,
      n_doses = n_doses,
      ncohort = ncohort,
      cohort_size = cohort_size,
      maxt = maxt,
      accrual_rate = accrual_rate,
      alpha1 = alpha1,
      alpha2 = alpha2,
      tite_distribution = tite_distribution,
      neli = neli,
      cutoff_eli = cutoff_eli,
      n_earlystop = n_earlystop,
      min_sample = min_sample
    )

    n_patients_all[trial, ] <- result$n_patients_by_dose
    n_tox_all[trial, ] <- result$n_tox_by_dose
    durations[trial] <- result$trial_duration

    if (!is.na(result$mtd)) {
      mtd_count[result$mtd] <- mtd_count[result$mtd] + 1
    }
  }

  cat("\nSimulation completed!\n\n")

  # Summarize results
  summary <- .summarize_simulation_titeboin(
    mtd_count = mtd_count,
    n_patients_all = n_patients_all,
    n_tox_all = n_tox_all,
    durations = durations,
    n_trials = n_trials
  )

  return(list(summary = summary))
}
