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
#'   eliminated in each trial.
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
#'     \item Check if lowest dose (dose 1) is eliminated first
#'     \item Identify admissible set: doses with patients AND not eliminated
#'     \item Extract isotonic estimates for admissible doses
#'     \item Select dose with estimate closest to target toxicity rate
#'     \item Apply tiebreaker by adding small dose-index increments
#'     \item If boundMTD = TRUE, ensure selected dose satisfies constraint
#'   }
#'
#'   The dose with isotonic estimate closest to the target is selected as MTD.
#'   Ties are broken by small perturbation (1e-10 * dose_index) preferring
#'   lower dose indices when estimates are equally close to target.
#'
#'   If the lowest dose (dose 1) is eliminated, NO MTD is selected regardless
#'   of other doses' status. This follows BOIN standard: "stop the trial if
#'   the lowest dose is eliminated due to toxicity, and no dose should be
#'   selected as the MTD."
#'
#' @examples
#' target <- 0.30
#' n_pts_mat <- matrix(c(3, 6, 9, 3, 6, 9), nrow = 2, byrow = TRUE)
#' n_tox_mat <- matrix(c(0, 1, 3, 0, 1, 2), nrow = 2, byrow = TRUE)
#' eliminated_mat <- matrix(FALSE, nrow = 2, ncol = 3)
#'
#' iso_est_mat <- isotonic_regression(n_pts_mat, n_tox_mat)
#' mtd_results <- select_mtd(iso_est_mat, n_pts_mat, eliminated_mat, target)
#' print(mtd_results)
#'
#' @importFrom stats pbeta
#'
#' @export
select_mtd <- function(iso_est_mat, n_pts_mat, eliminated_mat, target,
                       boundMTD = FALSE, lambda_d = NULL, min_mtd_sample = 1) {

  n_trials <- nrow(iso_est_mat)
  n_doses <- ncol(iso_est_mat)

  # Pre-allocate result vectors (faster than data.frame in loop)
  mtd_vec <- integer(n_trials)
  reason_vec <- character(n_trials)

  # Vectorized admissibility check
  admissible <- (n_pts_mat >= min_mtd_sample) & !eliminated_mat

  # Vectorized safety check for lowest dose
  lowest_dose_eliminated <- eliminated_mat[, 1]
  if (any(lowest_dose_eliminated)) {
    mtd_vec[lowest_dose_eliminated] <- NA_integer_
    reason_vec[lowest_dose_eliminated] <- "lowest_dose_eliminated"
  }

  # Trials still needing MTD selection
  active_trials <- which(!lowest_dose_eliminated)

  if (length(active_trials) == 0) {
    return(data.frame(
      trial = seq_len(n_trials),
      mtd = mtd_vec,
      reason = reason_vec,
      stringsAsFactors = FALSE
    ))
  }

  # Add perturbation for tiebreaking
  perturb_mat <- matrix(rep(seq_len(n_doses) * 1e-10, n_trials),
                        nrow = n_trials, byrow = TRUE)
  iso_perturbed <- iso_est_mat + perturb_mat

  # Calculate distance to target
  dist_mat <- abs(iso_perturbed - target)

  # Set inadmissible doses to Inf
  dist_mat[!admissible] <- Inf

  # Process each active trial
  for (trial in active_trials) {
    # Check if any admissible dose exists
    if (all(is.infinite(dist_mat[trial, ]))) {
      mtd_vec[trial] <- NA_integer_
      reason_vec[trial] <- "no_admissible_dose"
      next
    }

    # Find dose with minimum distance
    mtd_candidate <- which.min(dist_mat[trial, ])

    # Apply boundMTD constraint if needed
    if (boundMTD && !is.null(lambda_d)) {
      if (iso_est_mat[trial, mtd_candidate] > lambda_d) {
        # Find admissible doses satisfying lambda_d constraint
        valid_mask <- admissible[trial, ] & (iso_est_mat[trial, ] <= lambda_d)

        if (!any(valid_mask)) {
          mtd_vec[trial] <- NA_integer_
          reason_vec[trial] <- "no_dose_below_lambda_d"
          next
        }

        # Set invalid doses to Inf
        dist_constrained <- dist_mat[trial, ]
        dist_constrained[!valid_mask] <- Inf

        mtd_candidate <- which.min(dist_constrained)
      }
    }

    mtd_vec[trial] <- mtd_candidate
    reason_vec[trial] <- "trial_completed"
  }

  # Return as data.frame (single allocation at end)
  data.frame(
    trial = seq_len(n_trials),
    mtd = mtd_vec,
    reason = reason_vec,
    stringsAsFactors = FALSE
  )
}
