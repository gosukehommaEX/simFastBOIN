#' Summarize Simulated BOIN-ET Trials and Their OBDs
#'
#' @description
#'   Build the \code{boin_et_oc} object of \code{\link{sim_boin_et}} and
#'   \code{\link{sim_tite_boin_et}} from the simulated trials and the OBD
#'   selection.
#'
#' @param trials A \code{boin_et_trials} object.
#' @param selection The result of \code{\link{boin_et_select_obd}}.
#' @param est_method,obd_method The OBD selection used.
#' @param rule Named list of the parameters of the OBD selection, added to the
#'   settings.
#' @param keep_trials Logical scalar. Keep the trials and the OBDs.
#' @param call_expr The matched call.
#'
#' @return An object of class \code{boin_et_oc}.
#'
#' @noRd
boin_et_summarize <- function(trials, selection, est_method, obd_method, rule,
                              keep_trials, call_expr) {

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

  reasons <- boin_et_stop_reasons()
  stop_reasons <- as.numeric(table(factor(trials$stop_reason,
                                          levels = reasons))) / n_trials * 100
  names(stop_reasons) <- reasons

  settings <- c(set, rule)

  out <- list(
    design = trials$design,
    target = set$target,
    target_resp = set$target_resp,
    p_true = set$p_true,
    p_resp = set$p_resp,
    n_doses = n_doses,
    n_trials = n_trials,
    lambda1 = trials$boundary$lambda1,
    lambda2 = trials$boundary$lambda2,
    eta1 = trials$boundary$eta1,
    sel_percent = sel_percent,
    percent_no_obd = mean(is.na(selection$obd)) * 100,
    percent_early_stop = mean(trials$early_stop) * 100,
    n_pts_dose = n_pts_dose,
    n_tox_dose = n_tox_dose,
    n_resp_dose = n_resp_dose,
    total_n_pts = sum(n_pts_dose),
    total_n_tox = sum(n_tox_dose),
    total_n_resp = sum(n_resp_dose),
    stop_reason_percent = stop_reasons,
    est_method = est_method,
    obd_method = obd_method,
    settings = settings,
    call = call_expr
  )
  if (!is.null(trials$duration)) out$duration <- mean(trials$duration)
  if (keep_trials) {
    out$trials <- trials
    out$obd <- selection$obd
  }

  structure(out, class = "boin_et_oc")
}
