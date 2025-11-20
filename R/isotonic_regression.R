#' Isotonic Regression for Toxicity Rate Estimation
#'
#' @description
#'   Estimate toxicity rates at each dose level under the monotonicity constraint
#'   (toxicity increases with dose) using the Pool Adjacent Violators Algorithm (PAVA).
#'   Incorporates pseudocount (Beta-Binomial prior) to align with BOIN methodology.
#'   Uses vectorized computation for efficiency.
#'
#' @param n_pts
#'   Numeric vector. Number of patients treated at each dose level.
#'
#' @param n_tox
#'   Numeric vector. Number of patients with dose-limiting toxicity (DLT) at each dose level.
#'
#' @param min_sample
#'   Numeric. Minimum number of patients required for a dose to be considered
#'   for estimation. Default is 1. Doses with fewer patients return NA.
#'
#' @return
#'   Numeric vector of isotonic-adjusted toxicity rate estimates for each dose.
#'   Values are NA for doses with insufficient sample size (< min_sample).
#'
#' @details
#'   PAVA enforces the constraint that estimated toxicity rates are monotonically
#'   increasing across doses. Pseudocounts (0.05 to DLT count, 0.1 to total patients)
#'   are added before estimation, reflecting a Beta-Binomial conjugate prior.
#'   Patient counts are weighted by inverse variance in PAVA to account for
#'   estimation uncertainty.
#'
#'   The MTD selection is based on these isotonic-adjusted estimates, ensuring
#'   that the dose-toxicity relationship respects the natural monotonicity assumption.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' # Estimate isotonic toxicity rates after trial completion
#' n_pts <- c(3, 6, 9, 12)
#' n_tox <- c(0, 1, 3, 4)
#' iso_est <- isotonic_regression(n_pts, n_tox, min_sample = 3)
#' print(iso_est)
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

  # Apply PAVA
  iso_estimates_valid <- .pava_core(tox_rate_adj, variance_inv_weight)

  # Map back to original indices
  iso_est[valid_idx] <- iso_estimates_valid

  return(iso_est)
}


#' PAVA Core Implementation
#'
#' @description
#'   Pool Adjacent Violators Algorithm for isotonic regression.
#'   Enforces non-decreasing monotonicity constraint with weighted averaging.
#'
#' @param y Numeric vector. Adjusted toxicity rates (with pseudocounts).
#' @param w Numeric vector. Inverse variance weights.
#'
#' @return Numeric vector of isotonic-adjusted estimates with non-decreasing property.
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
