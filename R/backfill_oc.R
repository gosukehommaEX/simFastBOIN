#' Operating Characteristics from Simulated Backfill Trials
#'
#' @description
#'   Internal helper shared by \code{\link{sim_bf_boin}} and
#'   \code{\link{sim_be_boin}}. Select the MTD in every trial and summarize the
#'   trials in the layout of \code{\link{sim_boin}}, with the components that
#'   describe backfilling, timing and responses added.
#'
#' @param trials
#'   An object of class \code{backfill_trials}.
#'
#' @param bound_mtd,mtd_max_estimate,min_mtd_sample,overdose_cutoff
#'   As in \code{\link{sim_boin}}.
#'
#' @param keep_trials
#'   Logical scalar. Keep the trial by trial data in the result.
#'
#' @param call_expr
#'   The matched call of the calling function.
#'
#' @return
#'   An object of class \code{c("bf_boin_oc", "backfill_oc", "boin_oc")} or
#'   \code{c("be_boin_oc", "backfill_oc", "boin_oc")}, as described in
#'   \code{\link{sim_bf_boin}}.
#'
#' @noRd
backfill_oc <- function(trials, bound_mtd, mtd_max_estimate, min_mtd_sample,
                        overdose_cutoff, keep_trials, call_expr) {

  set <- trials$settings
  target <- set$target
  p_true <- set$p_true
  n_trials <- set$n_trials

  selection <- boin_select_mtd(
    n_pts = trials$n_pts, n_tox = trials$n_tox, target = target,
    cutoff_eli = set$cutoff_eli, extrasafe = set$extrasafe,
    offset = set$offset, bound_mtd = bound_mtd, p_tox = set$p_tox,
    mtd_max_estimate = mtd_max_estimate, min_mtd_sample = min_mtd_sample
  )

  # A trial stopped for safety selects no dose, whatever the final data show.
  safety_stop <- trials$stop_reason %in%
    c("lowest_dose_eliminated", "lowest_dose_too_toxic")
  selection$mtd[safety_stop] <- NA_integer_
  selection$reason[safety_stop] <- trials$stop_reason[safety_stop]

  n_doses <- length(p_true)
  dose_names <- paste0("DL", seq_len(n_doses))

  selected <- factor(selection$mtd, levels = seq_len(n_doses))
  sel_percent <- as.numeric(table(selected)) / n_trials * 100
  names(sel_percent) <- dose_names
  percent_no_mtd <- mean(is.na(selection$mtd)) * 100

  dose_means <- function(m) {
    out <- colMeans(m)
    names(out) <- dose_names
    out
  }

  if (is.null(overdose_cutoff)) overdose_cutoff <- target
  check_scalar_prob(overdose_cutoff, "overdose_cutoff")

  above_cutoff <- p_true > overdose_cutoff
  n_per_trial <- rowSums(trials$n_pts)
  max_pts <- set$max_total_pts
  n_above <- if (any(above_cutoff)) {
    rowSums(trials$n_pts[, above_cutoff, drop = FALSE])
  } else {
    rep(0L, nrow(trials$n_pts))
  }

  # Whether the dose the trial ended up recommending is itself above the cutoff.
  selected <- !is.na(selection$mtd)
  mtd_above <- rep(FALSE, length(selection$mtd))
  mtd_above[selected] <- above_cutoff[selection$mtd[selected]]

  overdose <- list(
    cutoff = overdose_cutoff,
    doses = which(above_cutoff),
    pct_patients = sum(n_above) / sum(n_per_trial) * 100,
    pct_patients_by_trial = mean(n_above / n_per_trial) * 100,
    avg_n_patients = mean(n_above),
    pct_trials_any = mean(n_above > 0) * 100,
    pct_trials_over_60 = mean(n_above > 0.6 * max_pts) * 100,
    pct_trials_over_80 = mean(n_above > 0.8 * max_pts) * 100,
    pct_trials_mtd_above = mean(mtd_above) * 100,
    pct_mtd_above_when_selected = if (any(selected)) {
      mean(mtd_above[selected]) * 100
    } else {
      NA_real_
    }
  )

  structure(
    list(
      sel_percent = sel_percent,
      percent_no_mtd = percent_no_mtd,
      n_pts_dose = dose_means(trials$n_pts),
      n_tox_dose = dose_means(trials$n_tox),
      total_n_pts = mean(rowSums(trials$n_pts)),
      total_n_tox = mean(rowSums(trials$n_tox)),
      overdose = overdose,
      stop_reason_percent = 100 * table(trials$stop_reason) / n_trials,
      n_bf_dose = dose_means(trials$n_bf),
      total_n_bf = mean(rowSums(trials$n_bf)),
      pct_trials_backfill = mean(rowSums(trials$n_bf) > 0L) * 100,
      n_resp_dose = dose_means(trials$n_resp),
      total_n_resp = mean(rowSums(trials$n_resp)),
      duration_mean = mean(trials$duration),
      duration_sd = if (n_trials > 1L) stats::sd(trials$duration) else NA_real_,
      pct_trials_suspended = mean(trials$n_suspensions > 0L) * 100,
      avg_n_suspensions = mean(trials$n_suspensions),
      avg_time_suspended = mean(trials$time_suspended),
      avg_n_turned_away = mean(trials$n_turned_away),
      target = target,
      p_true = p_true,
      p_resp = set$p_resp,
      p_saf = set$p_saf,
      p_tox = set$p_tox,
      lambda_e = trials$boundary$lambda_e,
      lambda_d = trials$boundary$lambda_d,
      n_doses = n_doses,
      n_trials = n_trials,
      settings = set,
      trials = if (keep_trials) {
        list(
          n_pts = trials$n_pts,
          n_tox = trials$n_tox,
          n_bf = trials$n_bf,
          n_resp = trials$n_resp,
          eliminated = trials$eliminated,
          cohorts_used = trials$cohorts_used,
          stop_reason = trials$stop_reason,
          duration = trials$duration,
          n_suspensions = trials$n_suspensions,
          time_suspended = trials$time_suspended,
          n_turned_away = trials$n_turned_away,
          mtd = selection$mtd,
          selection_reason = selection$reason
        )
      } else {
        NULL
      },
      call = call_expr
    ),
    class = c(paste0(set$design, "_oc"), "backfill_oc", "boin_oc")
  )
}
