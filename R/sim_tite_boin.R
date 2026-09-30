#' Operating Characteristics of a TITE-BOIN Design
#'
#' @description
#'   Simulate a time-to-event BOIN (TITE-BOIN) trial many times under one
#'   dose-toxicity scenario and summarize how often each dose is selected as the
#'   MTD, how many patients are treated and how many DLTs are observed at each
#'   dose, how long the trials last and how often accrual is suspended.
#'
#' @inheritParams tite_boin_simulate
#'
#' @param bound_mtd
#'   Logical scalar. Require the isotonic estimate at the selected dose to be at
#'   or below the de-escalation boundary. Defaults to \code{FALSE}.
#'
#' @param mtd_max_estimate
#'   Numeric scalar or \code{NULL}. Largest isotonic estimate a dose may have and
#'   still be selected as the MTD. See \code{\link{sim_boin}}.
#'
#' @param overdose_cutoff
#'   Numeric scalar or \code{NULL}. Doses whose true DLT probability exceeds this
#'   value count as overdoses in the \code{overdose} component of the result.
#'   Defaults to \code{NULL}, which uses \code{target}.
#'
#' @param min_mtd_sample
#'   Integer scalar. Smallest number of patients a dose must have received to be
#'   eligible as the MTD. Defaults to 1.
#'
#' @param keep_trials
#'   Logical scalar. Keep the full trial by trial data in the result.
#'   Defaults to \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Report progress while the simulation runs.
#'   Defaults to \code{FALSE}.
#'
#' @return
#'   An object of class \code{c("tite_boin_oc", "boin_oc")}, which has every
#'   component of the result of \code{\link{sim_boin}}, so that code written
#'   for \code{sim_boin()} also reads it, and in addition
#'   \item{duration_mean}{Average trial duration, from the first arrival until
#'     every enrolled patient has completed the assessment, in the time unit of
#'     \code{window} and \code{accrual_rate}.}
#'   \item{duration_sd}{Standard deviation of the trial duration.}
#'   \item{pct_trials_suspended}{Percentage of trials in which accrual was
#'     suspended at least once.}
#'   \item{avg_n_suspensions}{Average number of suspensions per trial.}
#'   \item{avg_time_suspended}{Average total time with accrual suspended per
#'     trial.}
#'   The \code{trials} component, when kept, also holds \code{duration},
#'   \code{n_suspensions} and \code{time_suspended} for every trial.
#'
#' @details
#'   The trials are simulated by \code{\link{tite_boin_simulate}}, whose details
#'   describe the conduct of a trial. Once every patient has completed the
#'   assessment the MTD is selected from the final data exactly as in
#'   \code{\link{sim_boin}}, with \code{\link{boin_select_mtd}}, and a trial
#'   stopped for safety selects no dose.
#'
#'   When no patient is ever pending at a decision, for example with
#'   \code{accrual = "fixed"} and \code{window} shorter than
#'   \code{1 / accrual_rate}, every result agrees with that of
#'   \code{\link{sim_boin}} under the same seed.
#'
#' @references
#'   Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K. E. (2018). Time-to-Event
#'   Bayesian Optimal Interval Design to Accelerate Phase I Trials. Clinical
#'   Cancer Research, 24(20), 4921-4930.
#'
#'   Lin, R. and Yuan, Y. (2020). Time-to-Event Model-Assisted Designs for
#'   Dose-Finding Trials with Delayed Toxicity. Biostatistics, 21(4), 807-824.
#'
#' @examples
#' # A three month window with two patients a month
#' oc <- sim_tite_boin(
#'   target = 0.30,
#'   p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   window = 3,
#'   accrual_rate = 2,
#'   n_trials = 500,
#'   seed = 123
#' )
#' oc
#'
#' \donttest{
#' # The effective sample size of Lin and Yuan (2020)
#' oc_ess <- sim_tite_boin(
#'   target = 0.30,
#'   p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
#'   n_cohort = 12,
#'   cohort_size = 3,
#'   window = 3,
#'   accrual_rate = 2,
#'   method = "ess",
#'   n_trials = 10000,
#'   seed = 123
#' )
#' oc_ess
#' }
#'
#' @seealso \code{\link{tite_boin_simulate}}, \code{\link{sim_tite_boin_multi}},
#'   \code{\link{sim_boin}}
#'
#' @export
sim_tite_boin <- function(target, p_true, n_cohort, cohort_size,
                          window, accrual_rate,
                          method = c("imputation", "ess"),
                          accrual = c("exponential", "uniform", "fixed"),
                          dlt_time = c("weibull", "uniform"),
                          late_fraction = 0.5,
                          max_pending_ratio = NULL, min_completed = NULL,
                          n_trials = 10000, start_dose = 1, n_earlystop = 18,
                          p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                          extrasafe = FALSE, offset = 0.05,
                          bound_mtd = FALSE, mtd_max_estimate = NULL,
                          min_mtd_sample = 1, overdose_cutoff = NULL,
                          n_earlystop_rule = c("with_stay", "simple"),
                          keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  method <- match.arg(method)
  accrual <- match.arg(accrual)
  dlt_time <- match.arg(dlt_time)
  n_earlystop_rule <- match.arg(n_earlystop_rule)

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- tite_boin_simulate(
    target = target, p_true = p_true, n_cohort = n_cohort,
    cohort_size = cohort_size, window = window, accrual_rate = accrual_rate,
    method = method, accrual = accrual, dlt_time = dlt_time,
    late_fraction = late_fraction, max_pending_ratio = max_pending_ratio,
    min_completed = min_completed, n_trials = n_trials,
    start_dose = start_dose, n_earlystop = n_earlystop, p_saf = p_saf,
    p_tox = p_tox, cutoff_eli = cutoff_eli, extrasafe = extrasafe,
    offset = offset, n_earlystop_rule = n_earlystop_rule, seed = seed
  )

  if (verbose) message("Selecting the MTD ...")

  selection <- boin_select_mtd(
    n_pts = trials$n_pts, n_tox = trials$n_tox, target = target,
    cutoff_eli = cutoff_eli, extrasafe = extrasafe, offset = offset,
    bound_mtd = bound_mtd, p_tox = trials$settings$p_tox,
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

  stop_reason_percent <- 100 * table(trials$stop_reason) / n_trials

  if (verbose) message("Done.")

  structure(
    list(
      sel_percent = sel_percent,
      percent_no_mtd = percent_no_mtd,
      n_pts_dose = n_pts_dose,
      n_tox_dose = n_tox_dose,
      total_n_pts = mean(rowSums(trials$n_pts)),
      total_n_tox = mean(rowSums(trials$n_tox)),
      overdose = overdose,
      stop_reason_percent = stop_reason_percent,
      duration_mean = mean(trials$duration),
      duration_sd = if (n_trials > 1L) stats::sd(trials$duration) else NA_real_,
      pct_trials_suspended = mean(trials$n_suspensions > 0L) * 100,
      avg_n_suspensions = mean(trials$n_suspensions),
      avg_time_suspended = mean(trials$time_suspended),
      target = target,
      p_true = p_true,
      p_saf = trials$settings$p_saf,
      p_tox = trials$settings$p_tox,
      lambda_e = trials$boundary$lambda_e,
      lambda_d = trials$boundary$lambda_d,
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
          duration = trials$duration,
          n_suspensions = trials$n_suspensions,
          time_suspended = trials$time_suspended,
          mtd = selection$mtd,
          selection_reason = selection$reason
        )
      } else {
        NULL
      },
      call = call_expr
    ),
    class = c("tite_boin_oc", "boin_oc")
  )
}
