#' MTD Selection for Multiple Trials
#'
#' @description
#'   Select MTD based on isotonic regression estimates.
#'   Replicates BOIN::select.mtd() for multiple trials.
#'
#' @param iso_est_mat
#'   Numeric matrix (n_trials x n_doses). Isotonic regression estimates.
#'
#' @param n_pts_mat
#'   Numeric matrix (n_trials x n_doses). Number of patients at each dose.
#'
#' @param eliminated_mat
#'   Logical matrix (n_trials x n_doses). Whether each dose is eliminated.
#'
#' @param target
#'   Numeric. Target toxicity rate.
#'
#' @param boundMTD
#'   Logical. If TRUE, MTD must satisfy lambda_d constraint.
#'
#' @param lambda_d
#'   Numeric. De-escalation boundary.
#'
#' @param min_mtd_sample
#'   Numeric. Minimum sample size for MTD consideration.
#'
#' @return
#'   Data frame with columns: trial, mtd, reason.
#'
#' @export
select_mtd <- function(iso_est_mat, n_pts_mat, eliminated_mat, target,
                       boundMTD = FALSE, lambda_d = NULL, min_mtd_sample = 1) {

  n_trials <- nrow(iso_est_mat)
  n_doses <- ncol(iso_est_mat)

  mtd_vector <- rep(NA_integer_, n_trials)
  reason_vector <- rep(NA_character_, n_trials)

  for (trial in seq_len(n_trials)) {

    # Admissible: not eliminated AND has patients
    adm_mask <- !eliminated_mat[trial, ] & (n_pts_mat[trial, ] > 0)
    adm_idx <- which(adm_mask)

    # No admissible doses
    if (length(adm_idx) == 0) {
      reason_vector[trial] <- "no_admissible"
      next
    }

    # Get phat for admissible doses
    phat <- iso_est_mat[trial, adm_idx]

    # Add perturbation
    phat_pert <- phat + (1:length(phat)) * 1e-10

    # Select closest to target
    selectd <- sort(abs(phat_pert - target), index.return = TRUE)$ix[1]
    mtd_cand <- adm_idx[selectd]

    # Apply boundMTD if needed
    if (boundMTD && !is.null(lambda_d)) {
      if (phat[selectd] > lambda_d) {
        below_mask <- phat <= lambda_d
        if (!any(below_mask)) {
          reason_vector[trial] <- "no_below_lambda_d"
          next
        }
        phat_filt <- phat[below_mask]
        phat_filt_pert <- phat_filt + (1:length(phat_filt)) * 1e-10
        selectd_filt <- sort(abs(phat_filt_pert - target), index.return = TRUE)$ix[1]
        mtd_cand <- adm_idx[below_mask][selectd_filt]
      }
    }

    mtd_vector[trial] <- mtd_cand
    reason_vector[trial] <- "completed"
  }

  data.frame(
    trial = seq_len(n_trials),
    mtd = mtd_vector,
    reason = reason_vector,
    stringsAsFactors = FALSE
  )
}
