#' Isotonic Regression for Toxicity Rate Estimation
#'
#' @description
#'   Estimate toxicity rates at each dose level under the monotonicity constraint
#'   using the Pool Adjacent Violators Algorithm (PAVA) without external packages.
#'   Supports both single dose set (vector input) and multiple dose sets (matrix input).
#'   Incorporates pseudocount (Beta-Binomial prior).
#'
#' @param n_pts
#'   Numeric vector (single dose set) or matrix (multiple dose sets).
#'   Number of patients treated at each dose level.
#'   - Vector: length n_doses
#'   - Matrix: nrow = n_trials, ncol = n_doses
#'
#' @param n_tox
#'   Numeric vector or matrix with same shape as n_pts.
#'   Number of patients with dose-limiting toxicity (DLT) at each dose level.
#'
#' @param min_sample
#'   Numeric. Minimum number of patients required for a dose to be considered.
#'   Default is 6. Doses with fewer patients return NA.
#'
#' @return
#'   - If n_pts is vector: Numeric vector of isotonic-adjusted toxicity rates
#'   - If n_pts is matrix: Matrix with same dimensions, isotonic estimates by row
#'
#' @details
#'   PAVA enforces non-decreasing monotonicity constraint. Pseudocounts
#'   (0.05 to DLT count, 0.1 to total patients) are added before estimation.
#'   Patient counts weighted by inverse variance are used as weights in PAVA.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' # Single dose set
#' n_pts <- c(3, 6, 9, 12)
#' n_tox <- c(0, 1, 3, 4)
#' iso_est <- isotonic_regression(n_pts, n_tox, min_sample = 3)
#' print(iso_est)
#'
#' # Multiple dose sets (matrix)
#' n_pts_matrix <- matrix(c(3, 6, 9, 12, 3, 6, 9, 12), nrow = 2, byrow = TRUE)
#' n_tox_matrix <- matrix(c(0, 1, 3, 4, 1, 2, 4, 5), nrow = 2, byrow = TRUE)
#' iso_est_matrix <- isotonic_regression(n_pts_matrix, n_tox_matrix, min_sample = 3)
#' print(iso_est_matrix)
#'
#' @export
isotonic_regression <- function(n_pts, n_tox, min_sample = 6) {

  # Dispatch based on input type
  if (is.vector(n_pts)) {
    return(.isotonic_regression_single(n_pts, n_tox, min_sample))
  } else if (is.matrix(n_pts) || is.data.frame(n_pts)) {
    return(.isotonic_regression_matrix(n_pts, n_tox, min_sample))
  } else {
    stop("n_pts must be a numeric vector, matrix, or data.frame")
  }
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

#' Single Dose Set Isotonic Regression
#'
#' @description
#'   Apply isotonic regression to a single dose set (vector).
#'   Handles NA values gracefully (e.g., eliminated doses).
#'
#' @param n_pts Numeric vector of patient counts at each dose.
#' @param n_tox Numeric vector of DLT counts at each dose.
#' @param min_sample Minimum sample size for dose inclusion.
#'
#' @return Numeric vector with isotonic-adjusted toxicity rates.
#'
#' @keywords internal
.isotonic_regression_single <- function(n_pts, n_tox, min_sample) {

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

#' Multiple Dose Sets Isotonic Regression
#'
#' @description
#'   Apply isotonic regression to multiple dose sets (matrix rows) in parallel.
#'   Each row is processed independently.
#'
#' @param n_pts Matrix or data.frame of patient counts (nrow = n_trials, ncol = n_doses).
#' @param n_tox Matrix or data.frame of DLT counts (same shape as n_pts).
#' @param min_sample Minimum sample size for dose inclusion.
#'
#' @return Matrix with same dimensions as input, isotonic estimates by row.
#'
#' @keywords internal
.isotonic_regression_matrix <- function(n_pts, n_tox, min_sample) {

  # Convert data.frame to matrix if needed
  if (is.data.frame(n_pts)) {
    n_pts <- as.matrix(n_pts)
  }
  if (is.data.frame(n_tox)) {
    n_tox <- as.matrix(n_tox)
  }

  n_sets <- nrow(n_pts)
  n_doses <- ncol(n_pts)

  # Initialize output matrix
  iso_matrix <- matrix(NA_real_, nrow = n_sets, ncol = n_doses)

  # Apply isotonic regression to each row independently
  for (i in seq_len(n_sets)) {
    iso_matrix[i, ] <- .isotonic_regression_single(
      n_pts[i, ],
      n_tox[i, ],
      min_sample
    )
  }

  # Preserve input structure
  dimnames(iso_matrix) <- dimnames(n_pts)

  return(iso_matrix)
}
