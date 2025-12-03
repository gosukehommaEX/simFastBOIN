#' Batch MTD Selection for Multiple Trials
#'
#' @description
#'   Select MTD for multiple trials using optimized batch processing. This function
#'   is used internally by \code{sim_boin()} for efficient MTD selection across
#'   simulation results.
#'
#' @param n_pts_mat Matrix (n_trials x n_doses). Number of patients treated at each
#'   dose level for each trial.
#' @param n_tox_mat Matrix (n_trials x n_doses). Number of DLTs at each dose level
#'   for each trial.
#' @param eliminated_mat Logical matrix (n_trials x n_doses). Whether each dose
#'   has been eliminated in each trial (used during dose escalation).
#' @param cohorts_completed Integer vector (n_trials). Number of cohorts completed
#'   in each trial.
#' @param stop_reason Character vector (n_trials). Pre-determined stopping reason
#'   for trials that stopped early (e.g., "lowest_dose_too_toxic").
#' @param target Numeric. Target toxicity probability.
#' @param boundMTD Logical. If TRUE, impose the condition that the isotonic estimate
#'   of toxicity probability for the selected MTD must be less than the de-escalation
#'   boundary. Default is FALSE.
#' @param lambda_d Numeric. De-escalation boundary. Only used if boundMTD = TRUE.
#' @param min_mtd_sample Numeric. Minimum number of patients required for a dose to be
#'   considered for MTD selection. Default is 1.
#' @param cutoff_eli Numeric. Cutoff for dose elimination. Default is 0.95.
#' @param extrasafe Logical. Apply extra safety rule at lowest dose. Default is FALSE.
#' @param offset Numeric. Offset for extrasafe rule. Default is 0.05.
#'
#' @return List with two elements:
#'   \item{mtd}{Integer vector (n_trials). Selected MTD for each trial (or NA)}
#'   \item{reason}{Character vector (n_trials). Reason for trial termination}
#'
#' @details
#'   This function follows the BOIN package's MTD selection algorithm:
#'   \enumerate{
#'     \item For each trial, compute elimination status for each dose using
#'       Pr(p > target | data) > cutoff_eli
#'     \item Create admissible set: doses with patients AND not eliminated
#'     \item Apply isotonic regression ONLY to admissible doses
#'     \item Select the dose with isotonic estimate closest to target
#'     \item Apply tiebreaker rules if multiple doses are equally close
#'     \item Apply boundMTD constraint if enabled
#'   }
#'
#'   **BOIN MTD Selection Rules:**
#'   \enumerate{
#'     \item Select the dose with isotonic estimate closest to target
#'     \item Tiebreaker: Add small increments (1e-10 * dose_index) to break ties
#'     \item Use sort() with index.return to find minimum distance
#'   }
#'
#'   **boundMTD Constraint:**
#'   When \code{boundMTD = TRUE}, an additional constraint is applied: the selected
#'   MTD's isotonic-estimated toxicity rate must be strictly less than or equal to the
#'   de-escalation boundary (lambda_d).
#'
#' @keywords internal
#' @importFrom stats pbeta
#' @export
select_mtd <- function(
    n_pts_mat, n_tox_mat, eliminated_mat,
    cohorts_completed, stop_reason, target,
    boundMTD = FALSE,
    lambda_d = NULL,
    min_mtd_sample = 1,
    cutoff_eli = 0.95,
    extrasafe = FALSE,
    offset = 0.05
) {
  n_trials <- nrow(n_pts_mat)
  n_doses <- ncol(n_pts_mat)

  mtd <- rep(NA_integer_, n_trials)
  reason <- rep(NA_character_, n_trials)

  # Identify stopping reasons that prevent MTD selection
  cannot_select_mtd_reasons <- c("lowest_dose_too_toxic", "lowest_dose_eliminated")
  has_terminal_stop_reason <- !is.na(stop_reason) & stop_reason %in% cannot_select_mtd_reasons

  if (any(has_terminal_stop_reason)) {
    reason[has_terminal_stop_reason] <- stop_reason[has_terminal_stop_reason]
    mtd[has_terminal_stop_reason] <- NA_integer_
  }

  # For trials without terminal stopping reasons, proceed with MTD selection
  can_select_mtd <- !has_terminal_stop_reason

  # Vectorized checks for remaining trials
  no_cohort <- cohorts_completed == 0

  # Set reasons using vectorized operations
  reason[can_select_mtd & no_cohort] <- "no_cohort_completed"

  # Find valid trials (those that can have MTD selected)
  valid_trials <- can_select_mtd & !no_cohort

  if (sum(valid_trials) > 0) {
    valid_idx <- which(valid_trials)

    for (trial in valid_idx) {
      n_pts <- n_pts_mat[trial, ]
      n_tox <- n_tox_mat[trial, ]

      # Compute elimination status using BOIN rule
      # elimi[i] = 1 if Pr(p[i] > target | data) > cutoff.eli
      elimi <- rep(0, n_doses)

      for (i in 1:n_doses) {
        if (n_pts[i] >= 3) {
          if (1 - pbeta(target, n_tox[i] + 1, n_pts[i] - n_tox[i] + 1) > cutoff_eli) {
            elimi[i:n_doses] <- 1
            break
          }
        }
      }

      # Apply extrasafe rule for lowest dose
      if (extrasafe) {
        if (n_pts[1] >= 3) {
          if (1 - pbeta(target, n_tox[1] + 1, n_pts[1] - n_tox[1] + 1) > (cutoff_eli - offset)) {
            elimi[1:n_doses] <- 1
          }
        }
      }

      # Check if lowest dose is eliminated
      if (elimi[1] == 1 || sum(n_pts[elimi == 0]) == 0) {
        mtd[trial] <- NA_integer_
        reason[trial] <- "lowest_dose_eliminated"
        next
      }

      # Create admissible set: not eliminated AND has patients
      adm_set <- (n_pts != 0) & (elimi == 0)
      adm_index <- which(adm_set)

      if (length(adm_index) == 0) {
        mtd[trial] <- NA_integer_
        reason[trial] <- "no_admissible_dose"
        next
      }

      # Extract data for admissible doses only
      y_adm <- n_tox[adm_index]
      n_adm <- n_pts[adm_index]

      # Apply isotonic regression to admissible doses ONLY
      phat <- (y_adm + 0.05) / (n_adm + 0.1)
      phat_var <- (y_adm + 0.05) * (n_adm - y_adm + 0.05) /
        ((n_adm + 0.1)^2 * (n_adm + 0.1 + 1))

      # Apply PAVA with inverse variance weights
      phat <- Iso::pava(phat, w = 1 / phat_var)

      # Add small increments to break ties (BOIN standard)
      phat <- phat + (1:length(phat)) * 1e-10

      # Apply boundMTD constraint if enabled
      if (boundMTD && !is.null(lambda_d)) {
        if (all(phat > lambda_d)) {
          mtd[trial] <- NA_integer_
          reason[trial] <- "no_dose_below_lambda_d"
          next
        }

        # BOIN implementation: filter phat and keep corresponding adm_index
        valid_mask <- phat <= lambda_d
        phat_filtered <- phat[valid_mask]
        adm_index_filtered <- adm_index[valid_mask]

        # Find closest to target in filtered set
        selectd <- sort(abs(phat_filtered - target), index.return = TRUE)$ix[1]
        mtd[trial] <- adm_index_filtered[selectd]
        reason[trial] <- "trial_completed"
      } else {
        # Find dose closest to target
        selectd <- sort(abs(phat - target), index.return = TRUE)$ix[1]
        mtd[trial] <- adm_index[selectd]
        reason[trial] <- "trial_completed"
      }
    }
  }

  return(list(mtd = mtd, reason = reason))
}
