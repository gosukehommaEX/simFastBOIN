#' Summarize Simulated U-BOIN Trials
#'
#' @description
#'   Build the \code{uboin_oc} object of \code{\link{sim_uboin}} from the
#'   trials of \code{\link{uboin_simulate}} and their OBDs. Internal.
#'
#' @param trials
#'   An object of class \code{uboin_trials}.
#'
#' @param obd
#'   Integer vector. The OBD of every trial, \code{NA} for none.
#'
#' @param keep_trials
#'   Logical scalar. Keep the trials in the result.
#'
#' @param call_expr
#'   The matched call.
#'
#' @return
#'   An object of class \code{uboin_oc}.
#'
#' @keywords internal
#' @noRd
uboin_summarize <- function(trials, obd, keep_trials, call_expr) {

  set <- trials$settings
  n_doses <- length(set$p_true)
  n_trials <- set$n_trials
  dose_names <- paste0("DL", seq_len(n_doses))

  sel_percent <- as.numeric(table(factor(obd, levels = seq_len(n_doses)))) /
    n_trials * 100
  names(sel_percent) <- dose_names
  n_pts_dose <- colMeans(trials$n_pts)
  n_tox_dose <- colMeans(trials$n_tox)
  n_resp_dose <- colMeans(trials$n_resp)
  names(n_pts_dose) <- dose_names
  names(n_tox_dose) <- dose_names
  names(n_resp_dose) <- dose_names

  true_utility <- as.numeric(set$p_outcome %*% set$utility)
  names(true_utility) <- dose_names

  reasons <- uboin_stop_reasons()
  stop_reason_percent <- as.numeric(table(factor(trials$stop_reason,
                                                 levels = reasons))) /
    n_trials * 100
  names(stop_reason_percent) <- reasons

  out <- list(
    design = trials$design,
    tox_upper = set$tox_upper,
    resp_lower = set$resp_lower,
    target = set$target,
    p_true = set$p_true,
    p_resp = set$p_resp,
    p_outcome = set$p_outcome,
    utility = set$utility,
    true_utility = true_utility,
    n_doses = n_doses,
    n_trials = n_trials,
    lambda_e = trials$boundary$lambda_e,
    lambda_d = trials$boundary$lambda_d,
    sel_percent = sel_percent,
    percent_no_obd = mean(is.na(obd)) * 100,
    percent_stage2 = mean(trials$stage2) * 100,
    n_pts_dose = n_pts_dose,
    n_tox_dose = n_tox_dose,
    n_resp_dose = n_resp_dose,
    total_n_pts = sum(n_pts_dose),
    total_n_tox = sum(n_tox_dose),
    total_n_resp = sum(n_resp_dose),
    stop_reason_percent = stop_reason_percent,
    settings = set,
    call = call_expr
  )
  if (keep_trials) {
    out$trials <- trials
    out$obd <- obd
  }
  structure(out, class = "uboin_oc")
}
