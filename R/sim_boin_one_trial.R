#' Simulate One BOIN Trial
#'
#' @description
#'   Conduct a complete simulation of one dose-finding trial using the BOIN design.
#'   The trial includes dose escalation/de-escalation and MTD selection phases.
#'   Optimized for computational efficiency using matrix indexing and vectorized operations.
#'
#' @param target
#'   Numeric. The target toxicity probability (e.g., 0.30 for 30%).
#'
#' @param p_true
#'   Numeric vector. True toxicity probabilities for each dose level.
#'
#' @param n_doses
#'   Numeric. Number of doses evaluated.
#'
#' @param n_cohort
#'   Numeric. Maximum number of cohorts allowed in the trial.
#'
#' @param cohort_size
#'   Numeric vector or scalar. Number of patients per cohort.
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
#'   Numeric. Sample size at current dose triggering trial termination.
#'   Default is 18. Trial stops when this sample size is reached at any dose.
#'
#' @param min_mtd_sample
#'   Numeric. Minimum number of patients required for a dose to be
#'   considered as MTD candidate. Default is 6.
#'
#' @return
#'   A list containing:
#'   \item{n_pts}{Numeric vector. Number of patients treated at each dose}
#'   \item{n_tox}{Numeric vector. Number of DLTs at each dose}
#'   \item{mtd}{Numeric. Selected MTD (dose level) or NA if no dose is suitable}
#'   \item{iso_est}{Numeric vector. Isotonic-adjusted toxicity rate estimates}
#'   \item{reason}{Character. Reason for trial termination}
#'   \item{cohorts_completed}{Numeric. Number of cohorts completed}
#'
#' @details
#'   The trial proceeds in three main phases:
#'   1. Dose Escalation/De-escalation: Cohorts are enrolled at the current dose,
#'      and dose adjustments are made based on the decision table.
#'   2. Safety Monitoring: At the lowest dose, if excessive toxicity is observed,
#'      the entire trial is stopped.
#'   3. MTD Selection: After trial completion or early stopping, isotonic regression
#'      is applied to estimate dose-toxicity relationships, and the MTD is selected
#'      as the dose with estimated toxicity rate closest to the target.
#'
#'   Reasons for trial termination include:
#'   - "trial_completed": Normal completion with cohort enrollment
#'   - "lowest_dose_too_toxic": Safety stopping at lowest dose
#'   - "lowest_dose_eliminated": Lowest dose eliminated due to toxicity
#'   - "no_valid_dose": All doses eliminated or insufficient data
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' target <- 0.30
#' p_true <- c(0.10, 0.25, 0.40)
#' boin_bound <- get_boin_boundary(target)
#' decision_table <- get_boin_decision(target, boin_bound$lambda_e,
#'                                 boin_bound$lambda_d, 18, 0.95)
#' stopping_boundaries <- get_boin_stopping_boundaries(target, 18, 0.90)
#'
#' result <- sim_boin_one_trial(
#'   target = target,
#'   p_true = p_true,
#'   n_doses = 3,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   decision_table = decision_table,
#'   stopping_boundaries = stopping_boundaries
#' )
#' print(result)
#' }
#'
#' @export
sim_boin_one_trial <- function(
    target,
    p_true,
    n_doses,
    n_cohort,
    cohort_size,
    decision_table,
    stopping_boundaries,
    n_earlystop = 18,
    min_mtd_sample = 1
) {

  # ========== Initialization ==========
  n_pts <- rep(0L, n_doses)        # Use integer type for patient counts
  n_tox <- rep(0L, n_doses)        # Use integer type for DLT counts
  current_dose <- 1L
  eliminated_doses <- rep(FALSE, n_doses)

  # Pre-compute cohort size vector to avoid repeated indexing
  if (length(cohort_size) == 1) {
    cohort_size <- rep(cohort_size, n_cohort)
  } else if (length(cohort_size) < n_cohort) {
    cohort_size <- c(cohort_size, rep(cohort_size[length(cohort_size)], n_cohort - length(cohort_size)))
  }
  cohort_size <- cohort_size[1:n_cohort]

  # Pre-compute maximum indices for table access to avoid repeated computation
  max_row_decision_table <- nrow(decision_table)
  max_col_decision_table <- ncol(decision_table)
  max_col_stopping_boundaries <- ncol(stopping_boundaries)

  # ========== Main Trial Loop ==========
  for (cohort in 1:n_cohort) {

    # Check early stopping condition
    if (n_pts[current_dose] >= n_earlystop) {
      break
    }

    # Get current cohort size
    current_cohort_size <- cohort_size[cohort]

    # Generate DLT data
    dlt_count <- rbinom(1, current_cohort_size, p_true[current_dose])

    # Update cumulative data
    n_pts[current_dose] <- n_pts[current_dose] + current_cohort_size
    n_tox[current_dose] <- n_tox[current_dose] + dlt_count

    # ========== Safety Stopping Rule at Lowest Dose ==========
    if (current_dose == 1L && n_pts[current_dose] >= 3L) {
      # Direct table access with bounds checking
      col_idx_stop <- min(n_pts[current_dose], max_col_stopping_boundaries)
      decision_stop <- stopping_boundaries[n_tox[current_dose] + 1L, col_idx_stop]

      if (!is.na(decision_stop) && decision_stop == "STOP") {
        return(list(
          n_pts = n_pts,
          n_tox = n_tox,
          mtd = NA_integer_,
          reason = "lowest_dose_too_toxic",
          cohorts_completed = cohort
        ))
      }
    }

    # ========== Obtain Dose Decision from Table ==========
    # Bounds check and direct table access
    if (n_pts[current_dose] <= max_col_decision_table) {
      decision <- decision_table[n_tox[current_dose] + 1L, n_pts[current_dose]]
    } else {
      decision <- NA_character_
    }

    # Default to "S" (Stay) if NA
    if (is.na(decision)) {
      decision <- "S"
    }

    # ========== Dose Adjustment Logic ==========
    switch(decision,
           "E" = {
             # Escalate if not at maximum dose and next dose not eliminated
             if (current_dose < n_doses && !eliminated_doses[current_dose + 1L]) {
               current_dose <- current_dose + 1L
             }
           },
           "D" = {
             # De-escalate to lower dose
             if (current_dose > 1L) {
               current_dose <- current_dose - 1L
               # Skip eliminated doses
               while (current_dose > 1L && eliminated_doses[current_dose]) {
                 current_dose <- current_dose - 1L
               }
             }
           },
           "DE" = {
             # Eliminate current dose and all higher doses
             eliminated_doses[current_dose:n_doses] <- TRUE

             if (current_dose > 1L) {
               current_dose <- current_dose - 1L
               while (current_dose > 1L && eliminated_doses[current_dose]) {
                 current_dose <- current_dose - 1L
               }
             } else {
               # Lowest dose eliminated
               return(list(
                 n_pts = n_pts,
                 n_tox = n_tox,
                 mtd = NA_integer_,
                 reason = "lowest_dose_eliminated",
                 cohorts_completed = cohort
               ))
             }
           }
    )
  }

  # ========== MTD Selection Phase ==========
  iso_est <- isotonic_regression(n_pts, n_tox, min_sample = min_mtd_sample)

  # Mark eliminated doses
  iso_est[eliminated_doses] <- NA_real_

  # Compute distance from target
  diffs <- abs(iso_est - target)

  # Check if any valid dose remains
  if (all(is.na(diffs))) {
    return(list(
      n_pts = n_pts,
      n_tox = n_tox,
      mtd = NA_integer_,
      reason = "no_valid_dose",
      cohorts_completed = n_cohort
    ))
  }

  # Find candidate dose(s) closest to target
  min_diff <- min(diffs, na.rm = TRUE)
  mtd_candidates <- which(diffs == min_diff, useNames = FALSE)

  # Tiebreaker for multiple candidates
  if (length(mtd_candidates) > 1L) {
    candidate_estimates <- iso_est[mtd_candidates]

    # Check if all above target
    all_above <- all(candidate_estimates > target, na.rm = TRUE)
    all_below <- all(candidate_estimates < target, na.rm = TRUE)

    mtd <- if (all_above) {
      min(mtd_candidates)
    } else if (all_below) {
      max(mtd_candidates)
    } else {
      max(mtd_candidates)
    }
  } else {
    mtd <- mtd_candidates[1L]
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
