#' Summarize Simulated Trials That Select an MTD
#'
#' @description
#'   Internal helper of \code{sim_iboin()} and \code{sim_gboin()}. Turn the
#'   patients, the toxicity counts or scores and the selected dose of every
#'   trial into the operating characteristics, in the layout of the result
#'   of \code{sim_boin()}.
#'
#' @param trials Simulated trials, with \code{n_pts}, \code{n_tox},
#'   \code{stop_reason} and \code{settings}.
#' @param selection Data frame with the columns \code{mtd} and \code{reason},
#'   one row per trial.
#' @param p_true True DLT probability, or true mean of the toxicity score,
#'   at each dose.
#' @param target Target of the design.
#' @param overdose_cutoff Doses with a true value above this cutoff count as
#'   overdoses; \code{NULL} for \code{target}.
#' @param keep_trials Whether to keep the trial level data.
#'
#' @return A list of the operating characteristics, without class.
#'
#' @noRd
mtd_oc_summary <- function(trials, selection, p_true, target,
                           overdose_cutoff = NULL, keep_trials = FALSE) {

  n_trials <- nrow(trials$n_pts)

  # A trial stopped for safety selects no dose, whatever the final data show.
  safety_stop <- trials$stop_reason %in%
    c("lowest_dose_eliminated", "lowest_dose_too_toxic")
  selection$mtd[safety_stop] <- NA_integer_
  selection$reason[safety_stop] <- trials$stop_reason[safety_stop]

  n_doses <- length(p_true)
  dose_names <- paste0("DL", seq_len(n_doses))

  sel_factor <- factor(selection$mtd, levels = seq_len(n_doses))
  sel_percent <- as.numeric(table(sel_factor)) / n_trials * 100
  names(sel_percent) <- dose_names
  percent_no_mtd <- mean(is.na(selection$mtd)) * 100

  n_pts_dose <- colMeans(trials$n_pts)
  n_tox_dose <- colMeans(trials$n_tox)
  names(n_pts_dose) <- dose_names
  names(n_tox_dose) <- dose_names

  if (is.null(overdose_cutoff)) overdose_cutoff <- target
  check_scalar_prob(overdose_cutoff, "overdose_cutoff")

  above_cutoff <- p_true > overdose_cutoff
  n_per_trial <- rowSums(trials$n_pts)
  max_pts <- trials$settings$max_total_pts
  n_above <- if (any(above_cutoff)) {
    rowSums(trials$n_pts[, above_cutoff, drop = FALSE])
  } else {
    rep(0L, n_trials)
  }

  # Whether the dose the trial ended up recommending is itself above the
  # cutoff.
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

  list(
    sel_percent = sel_percent,
    percent_no_mtd = percent_no_mtd,
    n_pts_dose = n_pts_dose,
    n_tox_dose = n_tox_dose,
    total_n_pts = mean(rowSums(trials$n_pts)),
    total_n_tox = mean(rowSums(trials$n_tox)),
    overdose = overdose,
    stop_reason_percent = 100 * table(trials$stop_reason) / n_trials,
    target = target,
    p_true = p_true,
    p_saf = trials$settings$p_saf,
    p_tox = trials$settings$p_tox,
    n_doses = n_doses,
    n_trials = n_trials,
    settings = trials$settings,
    trials = if (keep_trials) {
      list(
        n_pts = trials$n_pts,
        n_tox = trials$n_tox,
        eliminated = trials$eliminated,
        cohorts_used = trials$cohorts_used,
        stop_reason = trials$stop_reason,
        mtd = selection$mtd,
        selection_reason = selection$reason
      )
    } else {
      NULL
    }
  )
}
