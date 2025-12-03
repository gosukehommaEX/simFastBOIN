#' Isotonic Regression for Toxicity Rate Estimation
#'
#' @description
#'   Apply isotonic regression to estimate toxicity rates under the monotonicity
#'   constraint (toxicity increases with dose) using the Pool Adjacent Violators
#'   Algorithm (PAVA). This function processes multiple trials simultaneously for
#'   efficient batch computation.
#'
#' @param n_pts_mat Matrix (n_trials x n_doses). Number of patients treated at each
#'   dose level for each trial.
#' @param n_tox_mat Matrix (n_trials x n_doses). Number of DLTs at each dose level
#'   for each trial.
#' @param eliminated_mat Logical matrix (n_trials x n_doses). Whether each dose has
#'   been eliminated in each trial.
#' @param min_sample Numeric. Minimum number of patients required for a dose to be
#'   considered for estimation. Default is 1.
#'
#' @return Matrix (n_trials x n_doses) of isotonic-adjusted toxicity rate estimates.
#'   Values are NA for doses with insufficient sample size or eliminated doses.
#'
#'   **Performance Note:** When no valid doses exist for a trial (i.e., all doses
#'   have insufficient sample size or are eliminated), the function returns a row
#'   of NA values without performing PAVA computations. This early exit optimization
#'   significantly improves performance in scenarios where many trials terminate early
#'   or have extensive dose elimination, avoiding unnecessary computational overhead.
#'
#' @details
#'   PAVA enforces the constraint that estimated toxicity rates are monotonically
#'   increasing across doses. Pseudocounts (0.05 to DLT count, 0.1 to total patients)
#'   are added before estimation, reflecting a Beta-Binomial conjugate prior.
#'   Patient counts are weighted by inverse variance in PAVA to account for
#'   estimation uncertainty.
#'
#'   This implementation leverages \code{Iso::pava()}, which is implemented in C
#'   and highly optimized. While we loop over trials (unavoidable for PAVA), each
#'   trial's computation uses the fast C implementation, resulting in significant
#'   speedup compared to pure R PAVA implementations.
#'
#'   Key optimizations:
#'   \itemize{
#'     \item Use Iso::pava() (C implementation) instead of R loops
#'     \item Pre-allocate all vectors to avoid repeated allocations
#'     \item Vectorized pseudocount and variance calculations
#'     \item Early exit for trials with no valid doses
#'   }
#'
#'   The MTD selection is based on these isotonic-adjusted estimates, ensuring
#'   that the dose-toxicity relationship respects the natural monotonicity assumption.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507–523.
#'
#'   Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K.E. (2018). Time-to-event Bayesian
#'   Optimal Interval Design to Accelerate Phase I Trials. Clinical Cancer Research,
#'   24(20): 4921-4930.
#'
#' @importFrom Iso pava
#'
#' @examples
#' # Estimate isotonic toxicity rates for multiple trials
#' n_pts_mat <- matrix(c(3, 6, 9, 12,
#'                       3, 6, 9, 12,
#'                       3, 6, 9, 12), nrow = 3, byrow = TRUE)
#' n_tox_mat <- matrix(c(0, 1, 3, 4,
#'                       0, 0, 2, 3,
#'                       1, 2, 4, 6), nrow = 3, byrow = TRUE)
#' eliminated_mat <- matrix(FALSE, nrow = 3, ncol = 4)
#'
#' iso_est <- isotonic_regression(n_pts_mat, n_tox_mat, eliminated_mat, min_sample = 3)
#' print(iso_est)
#'
#' # For a single trial, provide 1-row matrices
#' n_pts_single <- matrix(c(3, 6, 9, 12), nrow = 1)
#' n_tox_single <- matrix(c(0, 1, 3, 4), nrow = 1)
#' eliminated_single <- matrix(FALSE, nrow = 1, ncol = 4)
#'
#' iso_est_single <- isotonic_regression(n_pts_single, n_tox_single,
#'                                       eliminated_single, min_sample = 3)
#' print(iso_est_single[1, ])  # Extract as vector
#'
#' @export
isotonic_regression <- function(
    n_pts_mat, n_tox_mat, eliminated_mat, min_sample = 1
) {
  n_trials <- nrow(n_pts_mat)
  n_doses <- ncol(n_pts_mat)

  # ========== Vectorized Pre-processing ==========
  # Identify all valid doses across all trials at once
  valid_doses_mat <- !is.na(n_pts_mat) & !is.na(n_tox_mat) &
    (n_pts_mat >= min_sample)

  # Vectorized pseudocount-adjusted toxicity rates for ALL trials
  tox_rate_adj_mat <- (n_tox_mat + 0.05) / (n_pts_mat + 0.1)

  # Vectorized variance calculation for ALL trials
  # Variance = (y + 0.05) * (n - y + 0.05) / ((n + 0.1)^2 * (n + 0.1 + 1))
  numerator_mat <- (n_tox_mat + 0.05) * (n_pts_mat - n_tox_mat + 0.05)
  denominator_mat <- ((n_pts_mat + 0.1) ^ 2) * (n_pts_mat + 0.1 + 1)
  variance_mat <- numerator_mat / denominator_mat

  # FIXED: Inverse variance weights
  # Weight = 1 / variance (not variance itself)
  variance_inv_weight_mat <- 1 / variance_mat

  # Pre-allocate output matrix
  iso_est_mat <- matrix(NA_real_, nrow = n_trials, ncol = n_doses)

  # ========== Apply Iso::pava() to each trial ==========
  # Only process trials with at least one valid dose
  has_valid <- rowSums(valid_doses_mat) > 0

  if (any(has_valid)) {
    trials_with_valid <- which(has_valid)

    # Apply Iso::pava() to each trial (unavoidable loop)
    for (trial in trials_with_valid) {
      valid_idx <- which(valid_doses_mat[trial, ])

      # Apply Iso::pava() - highly optimized C implementation
      iso_est_mat[trial, valid_idx] <- Iso::pava(
        tox_rate_adj_mat[trial, valid_idx],
        w = variance_inv_weight_mat[trial, valid_idx]
      )
    }
  }

  return(iso_est_mat)
}
