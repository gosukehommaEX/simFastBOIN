#' Summarize Simulated BOIN12 Trials and Their OBDs
#'
#' @description
#'   Build the \code{boin12_oc} object of \code{\link{sim_boin12}} from the
#'   simulated trials and the OBD selection.
#'
#' @param trials A \code{boin12_trials} object.
#' @param selection The result of \code{\link{boin12_select_obd}}.
#' @param est_method The estimate of the utility used.
#' @param keep_trials Logical scalar. Keep the trials and the OBDs.
#' @param call_expr The matched call.
#'
#' @return An object of class \code{boin12_oc}.
#'
#' @noRd
boin12_summarize <- function(trials, selection, est_method, keep_trials,
                             call_expr) {

  set <- trials$settings
  n_doses <- length(set$p_true)
  n_trials <- set$n_trials
  dose_names <- paste0("DL", seq_len(n_doses))

  selected <- factor(selection$obd, levels = seq_len(n_doses))
  sel_percent <- as.numeric(table(selected)) / n_trials * 100
  names(sel_percent) <- dose_names

  n_pts_dose <- colMeans(trials$n_pts)
  n_tox_dose <- colMeans(trials$n_tox)
  n_resp_dose <- colMeans(trials$n_resp)
  names(n_pts_dose) <- dose_names
  names(n_tox_dose) <- dose_names
  names(n_resp_dose) <- dose_names

  reasons <- boin12_stop_reasons()
  stop_reason_percent <- as.numeric(table(factor(trials$stop_reason,
                                                 levels = reasons))) /
    n_trials * 100
  names(stop_reason_percent) <- reasons

  out <- list(
    design = trials$design,
    target = set$target,
    target_resp = set$target_resp,
    p_true = set$p_true,
    p_resp = set$p_resp,
    n_doses = n_doses,
    n_trials = n_trials,
    lambda_e = trials$boundary$lambda_e,
    lambda_d = trials$boundary$lambda_d,
    u_benchmark = trials$boundary$u_benchmark,
    sel_percent = sel_percent,
    percent_no_obd = mean(is.na(selection$obd)) * 100,
    percent_early_stop = mean(trials$early_stop) * 100,
    n_pts_dose = n_pts_dose,
    n_tox_dose = n_tox_dose,
    n_resp_dose = n_resp_dose,
    total_n_pts = sum(n_pts_dose),
    total_n_tox = sum(n_tox_dose),
    total_n_resp = sum(n_resp_dose),
    stop_reason_percent = stop_reason_percent,
    est_method = est_method,
    settings = c(set, list(est_method = est_method)),
    call = call_expr
  )
  if (!is.null(trials$duration)) out$duration <- mean(trials$duration)
  if (keep_trials) {
    out$trials <- trials
    out$obd <- selection$obd
  }

  structure(out, class = "boin12_oc")
}
