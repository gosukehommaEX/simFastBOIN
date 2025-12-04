#' MTD Selection for Multiple Trials
#'
#' @description
#'   Select maximum tolerated dose (MTD) for multiple trials based on isotonic
#'   regression estimates and dose elimination status. This function follows the
#'   BOIN package's MTD selection algorithm.
#'
#' @usage
#'   select_mtd(iso_est_mat, n_pts_mat, eliminated_mat, target,
#'              boundMTD = FALSE, lambda_d = NULL, min_mtd_sample = 1)
#'
#' @param iso_est_mat
#'   Numeric matrix of size n_trials x n_doses. Isotonic regression estimates
#'   of toxicity rates from \code{\link{isotonic_regression}}.
#'
#' @param n_pts_mat
#'   Numeric matrix of size n_trials x n_doses. Number of patients treated at
#'   each dose level for each trial.
#'
#' @param eliminated_mat
#'   Logical matrix of size n_trials x n_doses. Whether each dose has been
#'   eliminated in each trial. Pre-computed by \code{\link{get_pts_and_tox}}.
#'
#' @param target
#'   Numeric. Target toxicity probability.
#'
#' @param boundMTD
#'   Logical. If TRUE, impose constraint that selected MTD's isotonic estimate
#'   must be <= lambda_d. Default is FALSE.
#'
#' @param lambda_d
#'   Numeric. De-escalation boundary. Required if boundMTD = TRUE.
#'
#' @param min_mtd_sample
#'   Numeric. Minimum number of patients required for a dose to be considered
#'   for MTD selection. Default is 1.
#'
#' @return
#'   Data frame with n_trials rows and three columns:
#'   \item{trial}{Integer. Trial ID (1 to n_trials)}
#'   \item{mtd}{Integer. Selected MTD dose number, or NA if no valid dose}
#'   \item{reason}{Character. Reason for trial completion or termination}
#'
#' @details
#'   For each trial, the function performs the following steps:
#'   \enumerate{
#'     \item Identify admissible set: doses with patients AND not eliminated
#'     \item Extract isotonic estimates for admissible doses
#'     \item Select dose with estimate closest to target toxicity rate
#'     \item Apply tiebreaker by adding small dose-index increments
#'     \item If boundMTD = TRUE, ensure selected dose satisfies constraint
#'   }
#'
#'   **BOIN MTD Selection Rules:**
#'   The dose with isotonic estimate closest to the target is selected as MTD.
#'   Ties are broken by small perturbation (1e-10 * dose_index) preferring
#'   lower dose indices when estimates are equally close to target.
#'
#' @examples
#' # Example: Select MTD after isotonic regression
#' n_trials <- 100
#' n_doses <- 5
#' target <- 0.25
#'
#' # Simulated isotonic regression results
#' iso_est_mat <- matrix(c(0.05, 0.15, 0.30, 0.35, 0.40), nrow = 100, ncol = 5, byrow = TRUE)
#'
#' # Number of patients
#' n_pts_mat <- matrix(c(3, 6, 9, 12, 9), nrow = 100, ncol = 5, byrow = TRUE)
#'
#' # Dose elimination status (none eliminated in this example)
#' eliminated_mat <- matrix(FALSE, nrow = 100, ncol = 5)
#'
#' # Select MTD without boundMTD constraint
#' mtd_results <- select_mtd(iso_est_mat, n_pts_mat, eliminated_mat, target = 0.25)
#' print(head(mtd_results))
#'
#' # Select MTD with boundMTD constraint
#' mtd_results_bounded <- select_mtd(iso_est_mat, n_pts_mat, eliminated_mat,
#'                                   target = 0.25, boundMTD = TRUE, lambda_d = 0.30)
#' print(head(mtd_results_bounded))
#'
#' @export
select_mtd <- function(iso_est_mat, n_pts_mat, eliminated_mat, target,
                       boundMTD = FALSE, lambda_d = NULL, min_mtd_sample = 1) {

  n_trials <- nrow(iso_est_mat)
  n_doses <- ncol(iso_est_mat)

  # Vectorized: identify admissible doses
  # Admissible = treated (n_pts > 0) AND not eliminated
  admissible <- (n_pts_mat > 0) & !eliminated_mat

  # Apply minimum sample requirement if specified
  if (min_mtd_sample > 1) {
    admissible <- admissible & (n_pts_mat >= min_mtd_sample)
  }

  # Pre-allocate output vectors
  mtd_vector <- rep(NA_integer_, n_trials)
  reason_vector <- rep(NA_character_, n_trials)

  # Process each trial (vectorized loop-free where possible)
  for (i in seq_len(n_trials)) {

    # Get admissible dose indices for trial i
    adm_idx <- which(admissible[i, ])

    # No admissible doses
    if (length(adm_idx) == 0) {
      reason_vector[i] <- "no_admissible_dose"
      next
    }

    # Extract isotonic estimates for admissible doses only
    iso_adm <- iso_est_mat[i, adm_idx]

    # Apply tiebreaker: add small increments by dose index
    iso_adm_perturb <- iso_adm + seq_along(iso_adm) * 1e-10

    # Calculate distance to target
    dist <- abs(iso_adm_perturb - target)

    # Find dose closest to target (using which.min for O(n) performance)
    best_pos <- which.min(dist)
    mtd_candidate <- adm_idx[best_pos]

    # Apply boundMTD constraint if requested
    if (boundMTD && !is.null(lambda_d)) {

      if (iso_est_mat[i, mtd_candidate] > lambda_d) {

        # Filter admissible doses that satisfy lambda_d constraint
        iso_adm_filtered <- iso_adm[iso_est_mat[i, adm_idx] <= lambda_d]
        adm_idx_filtered <- adm_idx[iso_est_mat[i, adm_idx] <= lambda_d]

        # No valid doses below lambda_d
        if (length(adm_idx_filtered) == 0) {
          reason_vector[i] <- "no_dose_below_lambda_d"
          next
        }

        # Apply perturbation to filtered admissible doses
        iso_filtered_perturb <- iso_adm_filtered + seq_along(iso_adm_filtered) * 1e-10
        dist_filtered <- abs(iso_filtered_perturb - target)

        # Select using which.min (O(n) performance)
        best_filtered <- which.min(dist_filtered)
        mtd_candidate <- adm_idx_filtered[best_filtered]
      }
    }

    mtd_vector[i] <- mtd_candidate
    reason_vector[i] <- "trial_completed"
  }

  # Combine into single data frame (one operation only)
  data.frame(
    trial = seq_len(n_trials),
    mtd = mtd_vector,
    reason = reason_vector,
    stringsAsFactors = FALSE
  )
}
