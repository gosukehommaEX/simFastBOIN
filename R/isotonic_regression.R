#' Isotonic Regression for Toxicity Rate Estimation
#'
#' @description
#'   Estimate toxicity rates at each dose level under the monotonicity constraint
#'   (toxicity increases with dose) using the Pool Adjacent Violators Algorithm (PAVA).
#'   Incorporates pseudocount (Beta-Binomial prior) to align with BOIN methodology.
#'   Uses vectorized computation for efficiency.
#'
#' @param n_pts Numeric vector. Number of patients treated at each dose level.
#' @param n_tox Numeric vector. Number of patients with dose-limiting toxicity (DLT)
#'   at each dose level.
#' @param min_sample Numeric. Minimum number of patients required for a dose to be
#'   considered for estimation. Default is 1. Doses with fewer patients return NA.
#'
#' @return Numeric vector of isotonic-adjusted toxicity rate estimates for each dose.
#'   Values are NA for doses with insufficient sample size (< min_sample).
#'
#' @details
#'   PAVA enforces the constraint that estimated toxicity rates are monotonically
#'   increasing across doses. Pseudocounts (0.05 to DLT count, 0.1 to total patients)
#'   are added before estimation, reflecting a Beta-Binomial conjugate prior.
#'   Patient counts are weighted by inverse variance in PAVA to account for
#'   estimation uncertainty.
#'
#'   The algorithm uses a backtracking approach with level pooling. When adjacent
#'   levels violate monotonicity, they are merged using weighted averages, and the
#'   algorithm backtracks to check if the new pooled level violates monotonicity
#'   with previous levels. This ensures global monotonicity while maintaining
#'   computational efficiency.
#'
#'   The MTD selection is based on these isotonic-adjusted estimates, ensuring
#'   that the dose-toxicity relationship respects the natural monotonicity assumption.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#'   Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K.E. (2018). Time-to-event Bayesian
#'   Optimal Interval Design to Accelerate Phase I Trials. Clinical Cancer Research,
#'   24(20): 4921-4930.
#'
#' @examples
#' # Estimate isotonic toxicity rates after trial completion
#' n_pts <- c(3, 6, 9, 12)
#' n_tox <- c(0, 1, 3, 4)
#' iso_est <- isotonic_regression(n_pts, n_tox, min_sample = 3)
#' print(iso_est)
#'
#' # Example with some doses having insufficient sample size
#' n_pts <- c(1, 3, 6, 9, 2)
#' n_tox <- c(0, 0, 2, 4, 1)
#' iso_est <- isotonic_regression(n_pts, n_tox, min_sample = 3)
#' print(iso_est)  # First and last doses return NA
#'
#' # Example demonstrating monotonicity constraint
#' n_pts <- c(6, 6, 6, 6)
#' n_tox <- c(1, 0, 2, 3)  # Raw rates: 0.167, 0, 0.333, 0.5
#' iso_est <- isotonic_regression(n_pts, n_tox)
#' print(iso_est)  # Adjusted to be non-decreasing
#'
#' @export
isotonic_regression <- function(n_pts, n_tox, min_sample = 1) {

  n_doses <- length(n_pts)
  iso_est <- rep(NA_real_, n_doses)

  # Identify doses that are not NA and have sufficient sample size
  valid_doses <- !is.na(n_pts) & !is.na(n_tox) & (n_pts >= min_sample)

  # Early exit if no valid doses
  if (sum(valid_doses) == 0) {
    return(iso_est)
  }

  # Extract valid data
  valid_idx <- which(valid_doses)
  n_pts_valid <- n_pts[valid_idx]
  n_tox_valid <- n_tox[valid_idx]

  # Apply pseudocounts: (y + 0.05) / (n + 0.1)
  tox_rate_adj <- (n_tox_valid + 0.05) / (n_pts_valid + 0.1)

  # Compute inverse variance weights
  numerator <- (n_tox_valid + 0.05) * (n_pts_valid - n_tox_valid + 0.05)
  denominator <- ((n_pts_valid + 0.1) ^ 2) * (n_pts_valid + 0.1 + 1)
  variance <- numerator / denominator
  variance_inv_weight <- 1 / variance

  # -------------------------------------------------------------------------
  # PAVA (Pool Adjacent Violators Algorithm) with backtracking
  # -------------------------------------------------------------------------
  n <- length(tox_rate_adj)

  if (n == 0) {
    iso_estimates_valid <- numeric(0)
  } else if (n == 1) {
    iso_estimates_valid <- tox_rate_adj
  } else {
    # Initialize: each element is its own level with index tracking
    level_estimates <- tox_rate_adj
    level_weights <- variance_inv_weight
    level_indices <- as.list(1:n)  # Track which original indices are in each level

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
        level_indices[[i]] <- c(level_indices[[i]], level_indices[[i + 1]])

        # Remove the merged level i+1
        level_estimates <- level_estimates[-(i + 1)]
        level_weights <- level_weights[-(i + 1)]
        level_indices <- level_indices[-(i + 1)]

        # Backtrack: check previous level
        if (i > 1) i <- i - 1
      } else {
        i <- i + 1
      }
    }

    # Map pooled estimates back to original indices
    iso_estimates_valid <- numeric(n)
    for (k in seq_along(level_estimates)) {
      iso_estimates_valid[level_indices[[k]]] <- level_estimates[k]
    }
  }

  # Map back to original dose indices
  iso_est[valid_idx] <- iso_estimates_valid

  return(iso_est)
}


#' Batch Isotonic Regression for Multiple Trials (Internal)
#'
#' @description
#'   Apply isotonic regression to multiple trials simultaneously. This is an internal
#'   function used by \code{sim_boin()} for efficient batch processing of simulation
#'   results.
#'
#' @param n_pts_mat Matrix (n_trials x n_doses). Number of patients treated at each
#'   dose level for each trial.
#' @param n_tox_mat Matrix (n_trials x n_doses). Number of DLTs at each dose level
#'   for each trial.
#' @param eliminated_mat Logical matrix (n_trials x n_doses). Whether each dose has
#'   been eliminated in each trial.
#' @param min_sample Numeric. Minimum number of patients required for a dose to be
#'   considered for estimation.
#'
#' @return Matrix (n_trials x n_doses) of isotonic-adjusted toxicity rate estimates.
#'   Values are NA for doses with insufficient sample size or eliminated doses.
#'
#' @details
#'   This function applies the same PAVA algorithm as \code{isotonic_regression()}
#'   but processes multiple trials in a single function call. While there is still
#'   a loop over trials internally, this approach is more efficient than calling
#'   \code{isotonic_regression()} separately for each trial due to reduced function
#'   call overhead and better memory locality.
#'
#' @keywords internal
.isotonic_regression_batch <- function(
    n_pts_mat, n_tox_mat, eliminated_mat, min_sample
) {
  n_trials <- nrow(n_pts_mat)
  n_doses <- ncol(n_pts_mat)
  iso_est_mat <- matrix(NA_real_, nrow = n_trials, ncol = n_doses)

  # Process each trial with optimized loop
  for (trial in 1:n_trials) {
    n_pts <- n_pts_mat[trial, ]
    n_tox <- n_tox_mat[trial, ]
    eliminated <- eliminated_mat[trial, ]

    # Identify valid doses
    valid_doses <- !eliminated & (n_pts >= min_sample)

    if (sum(valid_doses) == 0) {
      next
    }

    valid_idx <- which(valid_doses)
    n_pts_valid <- n_pts[valid_idx]
    n_tox_valid <- n_tox[valid_idx]

    # Pseudocount adjustment
    tox_rate_adj <- (n_tox_valid + 0.05) / (n_pts_valid + 0.1)

    # Variance weights
    numerator <- (n_tox_valid + 0.05) * (n_pts_valid - n_tox_valid + 0.05)
    denominator <- ((n_pts_valid + 0.1) ^ 2) * (n_pts_valid + 0.1 + 1)
    variance_inv_weight <- 1 / (numerator / denominator)

    # PAVA algorithm with index tracking
    if (length(tox_rate_adj) == 1) {
      iso_est_mat[trial, valid_idx] <- tox_rate_adj
    } else {
      n <- length(tox_rate_adj)
      level_estimates <- tox_rate_adj
      level_weights <- variance_inv_weight
      level_indices <- as.list(1:n)

      i <- 1
      while (i < length(level_estimates)) {
        if (level_estimates[i] > level_estimates[i + 1]) {
          # Pool adjacent violators
          pooled_weight <- level_weights[i] + level_weights[i + 1]
          pooled_estimate <- (level_weights[i] * level_estimates[i] +
                                level_weights[i + 1] * level_estimates[i + 1]) / pooled_weight

          level_estimates[i] <- pooled_estimate
          level_weights[i] <- pooled_weight
          level_indices[[i]] <- c(level_indices[[i]], level_indices[[i + 1]])

          # Remove i+1
          level_estimates <- level_estimates[-(i + 1)]
          level_weights <- level_weights[-(i + 1)]
          level_indices <- level_indices[-(i + 1)]

          # Backtrack
          if (i > 1) {
            i <- i - 1
          }
        } else {
          i <- i + 1
        }
      }

      # Expand pooled results back to original indices
      iso_result <- numeric(n)
      for (j in seq_along(level_estimates)) {
        iso_result[level_indices[[j]]] <- level_estimates[j]
      }

      iso_est_mat[trial, valid_idx] <- iso_result
    }
  }

  return(iso_est_mat)
}


#' Batch MTD Selection for Multiple Trials (Internal)
#'
#' @description
#'   Select MTD for multiple trials using optimized batch processing. This is an
#'   internal function used by \code{sim_boin()} for efficient MTD selection across
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
#'
#' @return List with two elements:
#'   \item{mtd}{Integer vector (n_trials). Selected MTD for each trial (or NA)}
#'   \item{reason}{Character vector (n_trials). Reason for trial termination}
#'
#' @details
#'   The function first handles trials with pre-determined stopping reasons (e.g.,
#'   safety stopping at lowest dose), then applies vectorized operations to determine
#'   termination reasons for remaining trials, and finally loops over valid trials to
#'   select MTD using the standard BOIN tiebreaker rules.
#'
#' @keywords internal
.select_mtd_batch <- function(
    diffs_mat, iso_est_mat, eliminated_mat,
    cohorts_completed, stop_reason, target
) {
  n_trials <- nrow(diffs_mat)
  n_doses <- ncol(diffs_mat)

  mtd <- rep(NA_integer_, n_trials)
  reason <- rep(NA_character_, n_trials)

  # First, handle trials with pre-determined stopping reasons
  has_stop_reason <- !is.na(stop_reason)
  if (any(has_stop_reason)) {
    reason[has_stop_reason] <- stop_reason[has_stop_reason]
    mtd[has_stop_reason] <- NA_integer_
  }

  # Vectorized checks for remaining trials
  no_cohort <- cohorts_completed == 0
  all_na <- apply(diffs_mat, 1, function(x) all(is.na(x)))
  lowest_eliminated <- eliminated_mat[, 1]

  # Set reasons using vectorized operations
  reason[!has_stop_reason & no_cohort] <- "no_cohort_completed"
  reason[!has_stop_reason & !no_cohort & all_na & lowest_eliminated] <- "lowest_dose_eliminated"
  reason[!has_stop_reason & !no_cohort & all_na & !lowest_eliminated] <- "no_valid_dose"

  # Find valid trials (those that can have MTD selected)
  valid_trials <- !has_stop_reason & !no_cohort & !all_na

  if (sum(valid_trials) > 0) {
    reason[valid_trials] <- "trial_completed"

    # For each valid trial, find MTD
    valid_idx <- which(valid_trials)

    for (trial in valid_idx) {
      diffs <- diffs_mat[trial, ]
      iso_est <- iso_est_mat[trial, ]

      min_diff <- min(diffs, na.rm = TRUE)
      mtd_candidates <- which(diffs == min_diff)

      if (length(mtd_candidates) == 1) {
        mtd[trial] <- mtd_candidates[1]
      } else {
        # Tiebreaker: apply BOIN rules
        candidate_estimates <- iso_est[mtd_candidates]
        all_above <- all(candidate_estimates > target, na.rm = TRUE)
        all_below <- all(candidate_estimates < target, na.rm = TRUE)

        if (all_above) {
          mtd[trial] <- min(mtd_candidates)
        } else if (all_below) {
          mtd[trial] <- max(mtd_candidates)
        } else {
          mtd[trial] <- max(mtd_candidates)
        }
      }
    }
  }

  return(list(mtd = mtd, reason = reason))
}
