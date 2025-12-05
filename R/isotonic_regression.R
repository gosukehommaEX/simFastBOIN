#' Isotonic Regression for Toxicity Rate Estimation (Optimized)
#'
#' @description
#'   Apply isotonic regression to estimate toxicity rates under the monotonicity
#'   constraint (toxicity increases with dose) using the Pool Adjacent Violators
#'   Algorithm (PAVA). This optimized version uses the C-based Iso::pava() for
#'   maximum performance.
#'
#' @usage
#'   isotonic_regression(n_pts_mat, n_tox_mat)
#'
#' @param n_pts_mat
#'   Numeric matrix of size n_trials x n_doses. Number of patients treated at
#'   each dose level for each trial.
#'
#' @param n_tox_mat
#'   Numeric matrix of size n_trials x n_doses. Number of dose-limiting toxicities
#'   (DLTs) observed at each dose level for each trial.
#'
#' @return
#'   Numeric matrix of size n_trials x n_doses. Isotonic regression estimates of
#'   toxicity rates at each dose for each trial. Values are NA for untreated
#'   doses (n_pts = 0).
#'
#' @details
#'   **PERFORMANCE OPTIMIZATION:**
#'   This version replaces the pure R `.pava_core()` implementation with `Iso::pava()`,
#'   which is implemented in C and is significantly faster (typically 5-10x speedup).
#'
#'   Pseudocounts (0.05 for toxicity, 0.1 for total patients) are added before
#'   estimation, reflecting a Beta-Binomial prior. PAVA is applied independently
#'   per trial using weighted pooling to enforce monotonicity.
#'
#' @references
#'   Liu, S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I
#'   Clinical Trials. \emph{Journal of the Royal Statistical Society: Series C},
#'   64, 507–523.
#'
#' @examples
#' # Example 1: Single trial
#' n_pts_single <- matrix(c(3, 6, 9, 12), nrow = 1)
#' n_tox_single <- matrix(c(0, 1, 3, 4), nrow = 1)
#'
#' iso_est <- isotonic_regression(n_pts_single, n_tox_single)
#' print(iso_est)
#'
#' # Example 2: Multiple trials
#' n_pts_mat <- matrix(c(3, 6, 9, 12,
#'                       3, 6, 9, 12,
#'                       3, 6, 9, 12), nrow = 3, byrow = TRUE)
#' n_tox_mat <- matrix(c(0, 1, 3, 4,
#'                       0, 0, 2, 3,
#'                       1, 2, 4, 6), nrow = 3, byrow = TRUE)
#'
#' iso_est <- isotonic_regression(n_pts_mat, n_tox_mat)
#' print(iso_est)
#'
#' @importFrom Iso pava
#'
#' @export
isotonic_regression <- function(n_pts_mat, n_tox_mat) {
  n_trials <- nrow(n_pts_mat)
  n_doses <- ncol(n_pts_mat)

  # Apply pseudocounts
  phat_mat <- (n_tox_mat + 0.05) / (n_pts_mat + 0.1)

  # Calculate inverse variance weights
  wt_mat <- ((n_pts_mat + 0.1)^2 * (n_pts_mat + 0.1 + 1)) /
    ((n_tox_mat + 0.05) * (n_pts_mat - n_tox_mat + 0.05))

  # Set zero weights for untreated doses
  wt_mat[n_pts_mat == 0] <- 0

  # Apply PAVA per trial using C implementation
  iso_est_list <- lapply(seq_len(n_trials), function(i) {
    trt_idx <- which(n_pts_mat[i, ] > 0)

    if (length(trt_idx) == 0) {
      return(rep(NA_real_, n_doses))
    }

    # Extract treated doses and apply PAVA using C implementation
    phat_trt <- phat_mat[i, trt_idx]
    wt_trt <- wt_mat[i, trt_idx]

    # USE Iso::pava() - C IMPLEMENTATION (FAST!)
    iso_trt <- Iso::pava(phat_trt, w = wt_trt)

    # Add perturbation for tiebreaking
    iso_trt <- iso_trt + seq_along(iso_trt) * 1e-10

    # Map back to full dose matrix
    result <- rep(NA_real_, n_doses)
    result[trt_idx] <- iso_trt

    return(result)
  })

  # Combine into matrix
  do.call(rbind, iso_est_list)
}
