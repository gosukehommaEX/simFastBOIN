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
#' # Complete workflow: get_pts_and_tox -> isotonic_regression -> select_mtd
#' target <- 0.30
#' p_true <- c(0.05, 0.10, 0.20, 0.30, 0.45)
#'
#' # Step 1: Generate patient and toxicity data
#' result <- get_pts_and_tox(
#'   n_trials = 10,
#'   target = target,
#'   p_true = p_true,
#'   n_cohort = 5,
#'   cohort_size = 3,
#'   n_earlystop = 18,
#'   seed = 123
#' )
#'
#' n_pts_mat <- result$n_pts_all
#' n_tox_mat <- result$n_tox_all
#' eliminated_mat <- result$eliminated_mat
#'
#' # Step 2: Apply isotonic regression
#' iso_est_mat <- isotonic_regression(n_pts_mat, n_tox_mat)
#'
#' # Step 3: Select MTD
#' mtd_results <- select_mtd(iso_est_mat, n_pts_mat, eliminated_mat, target = target)
#' print(mtd_results)
#'
#' # Step 3b: Select MTD with boundMTD constraint
#' lambda_d <- 0.35
#' mtd_results_bounded <- select_mtd(
#'   iso_est_mat, n_pts_mat, eliminated_mat,
#'   target = target, boundMTD = TRUE, lambda_d = lambda_d
#' )
#' print(mtd_results_bounded)
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

  # Select MTD for each trial
  results <- lapply(seq_len(n_trials), function(i) {

    # Get admissible dose indices for trial i
    adm_idx <- which(admissible[i, ])

    # No admissible doses
    if (length(adm_idx) == 0) {
      return(data.frame(
        trial = i,
        mtd = NA_integer_,
        reason = "no_admissible_dose",
        stringsAsFactors = FALSE
      ))
    }

    # Extract isotonic estimates for admissible doses only
    iso_adm <- iso_est_mat[i, adm_idx]

    # Apply tiebreaker: add small increments by dose index
    iso_adm_perturb <- iso_adm + seq_along(iso_adm) * 1e-10

    # Calculate distance to target
    dist <- abs(iso_adm_perturb - target)

    # Find dose closest to target using sort (matches BOIN)
    best_pos <- sort(dist, index.return = TRUE)$ix[1]
    mtd_candidate <- adm_idx[best_pos]

    # Apply boundMTD constraint if requested
    if (boundMTD && !is.null(lambda_d)) {

      if (iso_est_mat[i, mtd_candidate] > lambda_d) {

        # Filter admissible doses that satisfy lambda_d constraint
        iso_adm_filtered <- iso_adm[iso_est_mat[i, adm_idx] <= lambda_d]
        adm_idx_filtered <- adm_idx[iso_est_mat[i, adm_idx] <= lambda_d]

        # No valid doses below lambda_d
        if (length(adm_idx_filtered) == 0) {
          return(data.frame(
            trial = i,
            mtd = NA_integer_,
            reason = "no_dose_below_lambda_d",
            stringsAsFactors = FALSE
          ))
        }

        # Apply perturbation to filtered admissible doses
        iso_filtered_perturb <- iso_adm_filtered + seq_along(iso_adm_filtered) * 1e-10
        dist_filtered <- abs(iso_filtered_perturb - target)

        # Select using sort (matches BOIN)
        best_filtered <- sort(dist_filtered, index.return = TRUE)$ix[1]
        mtd_candidate <- adm_idx_filtered[best_filtered]
      }
    }

    return(data.frame(
      trial = i,
      mtd = mtd_candidate,
      reason = "trial_completed",
      stringsAsFactors = FALSE
    ))
  })

  # Combine results from all trials
  do.call(rbind, results)
}
