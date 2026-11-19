#' Summarize Simulated BOIN-ETC Trials and Their Optimal Dose Combinations
#'
#' @description
#'   Build the \code{boin_etc_oc} object of \code{\link{sim_boin_etc}} from
#'   the simulated trials and the selection.
#'
#' @param trials A \code{boin_etc_trials} object.
#' @param selection The result of \code{\link{boin_etc_select_obd}}.
#' @param psi00,psi11,utility_prior,odc_tie The selection settings used.
#' @param keep_trials Logical scalar. Keep the trials and the selections.
#' @param call_expr The matched call.
#'
#' @return An object of class \code{boin_etc_oc}.
#'
#' @noRd
boin_etc_summarize <- function(trials, selection, psi00, psi11,
                               utility_prior, odc_tie, keep_trials,
                               call_expr) {

  set <- trials$settings
  p_true <- set$p_true
  p_resp <- set$p_resp
  n_rows <- nrow(p_true)
  n_cols <- ncol(p_true)
  n_trials <- set$n_trials
  labels <- list(paste0("A", seq_len(n_rows)), paste0("B", seq_len(n_cols)))

  selected <- !is.na(selection$obd_a)
  if (odc_tie == "all") {
    sel_percent <- matrix(apply(array(as.numeric(selection$tied),
                                      dim(selection$tied)), c(2L, 3L), mean) *
                            100, n_rows, n_cols, dimnames = labels)
  } else {
    cell <- selection$obd_a[selected] +
      n_rows * (selection$obd_b[selected] - 1L)
    sel_percent <- matrix(tabulate(cell, nbins = n_rows * n_cols) /
                            n_trials * 100, n_rows, n_cols, dimnames = labels)
  }

  avg <- function(x) {
    matrix(apply(array(as.numeric(x), dim(x)), c(2L, 3L), mean), n_rows,
           n_cols, dimnames = labels)
  }
  n_pts_dose <- avg(trials$n_pts)
  n_tox_dose <- avg(trials$n_tox)
  n_resp_dose <- avg(trials$n_resp)

  reasons <- boin_etc_stop_reasons()
  stop_reason_percent <- as.numeric(table(factor(trials$stop_reason,
                                                 levels = reasons))) /
    n_trials * 100
  names(stop_reason_percent) <- reasons

  obd <- comb_boin12_true_obd(p_true, p_resp, set$target, psi00, psi11,
                              set$resp_cor, set$resp_cor_type)
  true_obd <- which(obd, arr.ind = TRUE)
  dimnames(true_obd) <- list(NULL, c("dose_a", "dose_b"))
  above <- p_true > set$target + 1e-12
  above_cells <- which(above, arr.ind = TRUE)
  dimnames(above_cells) <- list(NULL, c("dose_a", "dose_b"))
  total_n_pts <- sum(n_pts_dose)

  out <- list(
    design = trials$design,
    target = set$target,
    target_resp = set$target_resp,
    p_true = p_true,
    p_resp = p_resp,
    n_rows = n_rows,
    n_cols = n_cols,
    n_trials = n_trials,
    lambda1 = trials$boundary$lambda1,
    lambda2 = trials$boundary$lambda2,
    eta1 = trials$boundary$eta1,
    sel_percent = sel_percent,
    percent_no_obd = mean(!selected) * 100,
    percent_early_stop = mean(trials$early_stop) * 100,
    n_pts_dose = n_pts_dose,
    n_tox_dose = n_tox_dose,
    n_resp_dose = n_resp_dose,
    total_n_pts = total_n_pts,
    total_n_tox = sum(n_tox_dose),
    total_n_resp = sum(n_resp_dose),
    mean_subtrials = mean(trials$n_subtrials),
    true_obd = true_obd,
    pct_true_obd_selected = if (any(obd)) sum(sel_percent[obd]) else NA_real_,
    pct_patients_true_obd = if (any(obd) && total_n_pts > 0) {
      sum(n_pts_dose[obd]) / total_n_pts * 100
    } else {
      NA_real_
    },
    overdose = list(
      cutoff = set$target,
      combinations = above_cells,
      pct_selected = sum(sel_percent[above]),
      pct_patients = if (total_n_pts > 0) {
        sum(n_pts_dose[above]) / total_n_pts * 100
      } else {
        NA_real_
      }
    ),
    stop_reason_percent = stop_reason_percent,
    odc_tie = odc_tie,
    settings = c(set, list(psi00 = psi00, psi11 = psi11,
                           utility_prior = utility_prior, odc_tie = odc_tie)),
    call = call_expr
  )
  if (keep_trials) {
    out$trials <- trials
    out$obd <- cbind(dose_a = selection$obd_a, dose_b = selection$obd_b)
    out$tied <- selection$tied
  }

  structure(out, class = "boin_etc_oc")
}
