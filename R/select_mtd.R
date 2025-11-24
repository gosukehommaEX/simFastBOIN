#' Batch MTD Selection for Multiple Trials
#'
#' @description
#'   Select MTD for multiple trials using optimized batch processing. This function
#'   is used internally by \code{sim_boin()} for efficient MTD selection across
#'   simulation results.
#'
#' @param diffs_mat Matrix (n_trials x n_doses). Absolute differences between
#'   isotonic estimates and target toxicity rate for each trial.
#' @param iso_est_mat Matrix (n_trials x n_doses). Isotonic-adjusted toxicity
#'   rate estimates for each trial.
#' @param eliminated_mat Logical matrix (n_trials x n_doses). Whether each dose
#'   has been eliminated in each trial.
#' @param cohorts_completed Integer vector (n_trials). Number of cohorts completed
#'   in each trial.
#' @param stop_reason Character vector (n_trials). Pre-determined stopping reason
#'   for trials that stopped early (e.g., "lowest_dose_too_toxic").
#' @param target Numeric. Target toxicity probability.
#' @param boundMTD Logical. If TRUE, impose the condition that the isotonic estimate
#'   of toxicity probability for the selected MTD must be less than the de-escalation
#'   boundary. Default is FALSE.
#' @param lambda_d Numeric. De-escalation boundary. Only used if boundMTD = TRUE.
#'
#' @return List with two elements:
#'   \item{mtd}{Integer vector (n_trials). Selected MTD for each trial (or NA)}
#'   \item{reason}{Character vector (n_trials). Reason for trial termination}
#'
#' @details
#'   The function first handles trials with pre-determined stopping reasons that
#'   prevent MTD selection (e.g., safety stopping at lowest dose, lowest dose
#'   eliminated), then applies vectorized operations to determine termination
#'   reasons for remaining trials, and finally loops over valid trials to select
#'   MTD using the standard BOIN tiebreaker rules.
#'
#'   Trials that stopped due to reaching n_earlystop are treated as normally
#'   completed trials and undergo MTD selection.
#'
#'   **BOIN MTD Selection Rules:**
#'   \enumerate{
#'     \item Select the dose with isotonic estimate closest to target
#'     \item Tiebreaker when multiple doses are equally close:
#'       \itemize{
#'         \item If all tied doses are above target: select lowest dose
#'         \item If all tied doses are below target: select highest dose
#'         \item If tied doses are on both sides: select highest dose
#'       }
#'   }
#'
#'   **boundMTD Constraint:**
#'   When \code{boundMTD = TRUE}, an additional constraint is applied: the selected
#'   MTD's isotonic-estimated toxicity rate must be strictly less than the
#'   de-escalation boundary (lambda_d). If the initially selected MTD violates this
#'   constraint, the function searches for a lower dose that satisfies the constraint.
#'
#' @keywords internal
#' @export
select_mtd <- function(
    diffs_mat, iso_est_mat, eliminated_mat,
    cohorts_completed, stop_reason, target,
    boundMTD = FALSE,
    lambda_d = NULL
) {
  n_trials <- nrow(diffs_mat)
  n_doses <- ncol(diffs_mat)

  mtd <- rep(NA_integer_, n_trials)
  reason <- rep(NA_character_, n_trials)

  # Identify stopping reasons that prevent MTD selection
  # "n_earlystop_reached" and "n_earlystop_with_stay" should still allow MTD selection
  cannot_select_mtd_reasons <- c("lowest_dose_too_toxic", "lowest_dose_eliminated")
  has_terminal_stop_reason <- !is.na(stop_reason) & stop_reason %in% cannot_select_mtd_reasons

  if (any(has_terminal_stop_reason)) {
    reason[has_terminal_stop_reason] <- stop_reason[has_terminal_stop_reason]
    mtd[has_terminal_stop_reason] <- NA_integer_
  }

  # For trials without terminal stopping reasons, proceed with MTD selection
  # This includes trials with NA stop_reason, "n_earlystop_reached", or "n_earlystop_with_stay"
  can_select_mtd <- !has_terminal_stop_reason

  # Vectorized checks for remaining trials
  no_cohort <- cohorts_completed == 0
  all_na <- apply(diffs_mat, 1, function(x) all(is.na(x)))
  lowest_eliminated <- eliminated_mat[, 1]

  # Set reasons using vectorized operations (only for trials that can potentially select MTD)
  reason[can_select_mtd & no_cohort] <- "no_cohort_completed"
  reason[can_select_mtd & !no_cohort & all_na & lowest_eliminated] <- "lowest_dose_eliminated"
  reason[can_select_mtd & !no_cohort & all_na & !lowest_eliminated] <- "no_valid_dose"

  # Find valid trials (those that can have MTD selected)
  valid_trials <- can_select_mtd & !no_cohort & !all_na

  if (sum(valid_trials) > 0) {
    # For valid trials, set reason to "trial_completed" by default
    # This will be overridden if no dose satisfies boundMTD constraint
    reason[valid_trials] <- "trial_completed"

    # For each valid trial, find MTD
    valid_idx <- which(valid_trials)

    for (trial in valid_idx) {
      diffs <- diffs_mat[trial, ]
      iso_est <- iso_est_mat[trial, ]

      min_diff <- min(diffs, na.rm = TRUE)
      mtd_candidates <- which(diffs == min_diff)

      if (length(mtd_candidates) == 1) {
        mtd_temp <- mtd_candidates[1]
      } else {
        # Tiebreaker: apply BOIN rules
        candidate_estimates <- iso_est[mtd_candidates]
        all_above <- all(candidate_estimates > target, na.rm = TRUE)
        all_below <- all(candidate_estimates < target, na.rm = TRUE)

        if (all_above) {
          mtd_temp <- min(mtd_candidates)
        } else if (all_below) {
          mtd_temp <- max(mtd_candidates)
        } else {
          mtd_temp <- max(mtd_candidates)
        }
      }

      # Apply boundMTD constraint if enabled
      if (boundMTD && !is.null(lambda_d)) {
        # Check if selected MTD satisfies the constraint
        if (!is.na(iso_est[mtd_temp]) && iso_est[mtd_temp] >= lambda_d) {
          # Constraint violated: search for valid doses below lambda_d
          valid_doses <- which(!is.na(iso_est) &
                                 iso_est < lambda_d &
                                 !eliminated_mat[trial, ])

          if (length(valid_doses) == 0) {
            # No dose satisfies the constraint
            mtd[trial] <- NA_integer_
            reason[trial] <- "no_dose_below_lambda_d"
          } else {
            # Select the dose closest to target among valid doses
            valid_diffs <- abs(iso_est[valid_doses] - target)
            best_valid <- valid_doses[which.min(valid_diffs)]
            mtd[trial] <- best_valid
            # Reason remains "trial_completed"
          }
        } else {
          # Constraint satisfied
          mtd[trial] <- mtd_temp
        }
      } else {
        # boundMTD not enabled or lambda_d not provided
        mtd[trial] <- mtd_temp
      }
    }
  }

  return(list(mtd = mtd, reason = reason))
}
