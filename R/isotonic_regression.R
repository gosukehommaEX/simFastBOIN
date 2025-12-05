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
#'   per trial using weighted pooling to enforce monotonicity: dose-toxicity
#'   relationships must be non-decreasing across dose levels.
#'
#' @references
#'   Liu, S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I
#'   Clinical Trials. \emph{Journal of the Royal Statistical Society: Series C},
#'   64, 507-523.
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

  # Get dimensions of the input matrices
  n_trials <- nrow(n_pts_mat)
  n_doses <- ncol(n_pts_mat)

  # Apply pseudocounts to observed proportions
  # Adds prior mass: 0.05 for toxicity counts, 0.1 for total patient counts
  # This reflects a Beta-Binomial conjugate prior assumption
  phat_mat <- (n_tox_mat + 0.05) / (n_pts_mat + 0.1)

  # Calculate inverse variance weights for PAVA
  # Weights reflect precision: higher weights for trials with more data
  wt_mat <- '/'(
    (n_pts_mat + 0.1)^2 * (n_pts_mat + 0.1 + 1),
    (n_tox_mat + 0.05) * (n_pts_mat - n_tox_mat + 0.05)
  )

  # Set zero weights for untreated doses (n_pts = 0)
  # These doses are excluded from isotonic regression
  wt_mat[n_pts_mat == 0] <- 0

  # Apply PAVA independently for each trial
  iso_est_list <- lapply(seq_len(n_trials), function(i) {
    # Identify doses that have been treated in trial i
    trt_idx <- which(n_pts_mat[i, ] > 0)

    if (length(trt_idx) == 0) {
      # No doses treated: return NA vector
      return(rep(NA_real_, n_doses))
    }

    # Extract treated doses with pseudocounts and weights
    phat_trt <- phat_mat[i, trt_idx]
    wt_trt <- wt_mat[i, trt_idx]

    # Apply Pool Adjacent Violators Algorithm (PAVA)
    # Enforces monotonicity: toxicity rate must be non-decreasing with dose
    iso_trt <- Iso::pava(phat_trt, w = wt_trt)

    # Add small perturbation for tiebreaking
    # Ensures strict ordering when multiple doses have identical estimates
    iso_trt <- iso_trt + seq_along(iso_trt) * 1e-10

    # Map isotonic estimates back to full dose matrix
    # Untreated doses remain NA
    result <- rep(NA_real_, n_doses)
    result[trt_idx] <- iso_trt

    return(result)
  })

  # Combine results from all trials into matrix
  do.call(rbind, iso_est_list)
}
