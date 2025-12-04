#' Isotonic Regression for Toxicity Rate Estimation
#'
#' @description
#'   Apply isotonic regression to estimate toxicity rates under the monotonicity
#'   constraint (toxicity increases with dose) using the Pool Adjacent Violators
#'   Algorithm (PAVA).
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

  # Apply PAVA per trial
  iso_est_list <- lapply(seq_len(n_trials), function(i) {

    trt_idx <- which(n_pts_mat[i, ] > 0)

    if (length(trt_idx) == 0) {
      return(rep(NA_real_, n_doses))
    }

    # Extract treated doses and apply PAVA
    phat_trt <- phat_mat[i, trt_idx]
    wt_trt <- wt_mat[i, trt_idx]
    iso_trt <- .pava_core(phat_trt, wt_trt)

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

# Internal PAVA implementation with backtracking
.pava_core <- function(y, w) {

  n <- length(y)
  if (n <= 1) return(y)

  # Initialize level sets
  est <- y
  wt <- w
  sz <- rep(1, n)

  # Iteratively enforce monotonicity
  i <- 1
  while (i < length(est)) {
    if (est[i] > est[i + 1]) {
      # Merge adjacent levels
      wt_sum <- wt[i] + wt[i + 1]
      est[i] <- (wt[i] * est[i] + wt[i + 1] * est[i + 1]) / wt_sum
      wt[i] <- wt_sum
      sz[i] <- sz[i] + sz[i + 1]

      # Remove merged level and backtrack
      est <- est[-(i + 1)]
      wt <- wt[-(i + 1)]
      sz <- sz[-(i + 1)]
      if (i > 1) i <- i - 1
    } else {
      i <- i + 1
    }
  }

  # Map pooled estimates back to original indices
  result <- numeric(n)
  idx <- 1
  for (k in seq_along(est)) {
    result[idx:(idx + sz[k] - 1)] <- est[k]
    idx <- idx + sz[k]
  }

  return(result)
}
