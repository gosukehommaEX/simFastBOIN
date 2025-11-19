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
#'   for estimation. Default is 6. Doses with fewer patients return NA.
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
#' @importFrom Iso pava
#'
#' @examples
#' # Estimate isotonic toxicity rates after trial completion
#' n_pts <- c(3, 6, 9, 12)
#' n_tox <- c(0, 1, 3, 4)
#' iso_est <- isotonic_regression(n_pts, n_tox, min_sample = 3)
#' print(iso_est)
#'
#' @export
isotonic_regression <- function(n_pts, n_tox, min_sample = 6) {

  n_doses <- length(n_pts)
  iso_est <- rep(NA_real_, n_doses)

  # Identify valid doses in one operation
  valid_doses <- n_pts >= min_sample
  n_valid <- sum(valid_doses)

  # Early exit if no valid doses
  if (n_valid == 0) {
    return(iso_est)
  }

  # Vectorized computation for pseudocount-adjusted values
  # Only compute for valid doses to avoid unnecessary operations
  valid_idx <- which(valid_doses)

  # Pseudocount-adjusted toxicity rates: (y + 0.05) / (n + 0.1)
  tox_rate_adj <- (n_tox[valid_idx] + 0.05) / (n_pts[valid_idx] + 0.1)

  # Vectorized variance computation
  numerator <- (n_tox[valid_idx] + 0.05) * (n_pts[valid_idx] - n_tox[valid_idx] + 0.05)
  denominator <- ((n_pts[valid_idx] + 0.1) ^ 2) * (n_pts[valid_idx] + 0.1 + 1)
  variance <- numerator / denominator

  # Inverse variance weights in one operation
  variance_inv_weight <- 1 / variance

  # Apply PAVA to adjusted rates with inverse variance weights
  iso_est[valid_idx] <- Iso::pava(tox_rate_adj, w = variance_inv_weight)

  return(iso_est)
}
