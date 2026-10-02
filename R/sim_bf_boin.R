#' Operating Characteristics of a BF-BOIN Design
#'
#' @description
#'   Simulate the BOIN design with backfilling (BF-BOIN) of Zhao et al. (2024)
#'   many times under one scenario of DLT and response probabilities, and
#'   summarize how often each dose is selected as the MTD, how many patients are
#'   treated and backfilled at each dose, how many DLTs and responses occur, and
#'   how long the trials last.
#'
#' @inheritParams bf_boin_simulate
#' @inheritParams sim_tite_boin
#'
#' @return
#'   An object of class \code{c("bf_boin_oc", "backfill_oc", "boin_oc")}, which
#'   has every component of the result of \code{\link{sim_boin}}, so that code
#'   written for \code{sim_boin()} also reads it, counting backfilled patients
#'   among the patients treated, and in addition
#'   \item{n_bf_dose}{Average number of backfilled patients at each dose.}
#'   \item{total_n_bf}{Average number of backfilled patients per trial.}
#'   \item{pct_trials_backfill}{Percentage of trials with at least one
#'     backfilled patient.}
#'   \item{n_resp_dose}{Average number of patients with a response at each
#'     dose.}
#'   \item{total_n_resp}{Average number of patients with a response per trial.}
#'   \item{duration_mean}{Average trial duration, from the first arrival until
#'     every enrolled patient has completed the DLT assessment, in the time unit
#'     of \code{window} and \code{accrual_rate}.}
#'   \item{duration_sd}{Standard deviation of the trial duration.}
#'   \item{pct_trials_suspended}{Percentage of trials in which the dose
#'     escalation had to wait at least once.}
#'   \item{avg_n_suspensions}{Average number of decisions per trial at which the
#'     dose escalation had to wait.}
#'   \item{avg_time_suspended}{Average total time per trial during which the
#'     dose escalation waited.}
#'   \item{avg_n_turned_away}{Average number of patients turned away per trial,
#'     which is zero unless \code{no_slot = "leave"}.}
#'   \item{p_resp}{The true response probabilities.}
#'   The \code{trials} component, when kept, holds the trial by trial data of
#'   \code{\link{bf_boin_simulate}} together with the selected MTD.
#'
#' @details
#'   The trials are simulated by \code{\link{bf_boin_simulate}}, whose details
#'   describe the conduct of a trial. Once every patient has completed the
#'   assessment the MTD is selected from the data of all patients, backfilled
#'   ones included, exactly as in \code{\link{sim_boin}}, and a trial stopped for
#'   safety selects no dose. The overdose cutoffs of 60 and 80 percent refer to
#'   the size of the dose escalation, \code{n_cohort} cohorts.
#'
#'   When no dose is ever opened for backfilling, for example with every
#'   \code{p_resp} equal to 0, the selection and the patient counts agree with
#'   those of \code{\link{sim_boin}} under the same seed.
#'
#'   To reproduce Table 4 of Zhao et al. (2024), use \code{target = 0.25},
#'   \code{n_cohort = 10}, \code{cohort_size = 3}, \code{window = 1},
#'   \code{accrual_rate = 3}, \code{n_earlystop = 9},
#'   \code{stay_on_1_of_3 = TRUE}, \code{accrual = "uniform"} and
#'   \code{no_slot = "leave"}.
#'
#' @references
#'   Zhao, Y., Yuan, Y., Korn, E. L. and Freidlin, B. (2024). Backfilling
#'   Patients in Phase I Dose-Escalation Trials Using Bayesian Optimal Interval
#'   Design (BOIN). Clinical Cancer Research, 30(4), 673-679.
#'
#' @examples
#' # Scenario 3 of Zhao et al. (2024)
#' oc <- sim_bf_boin(
#'   target = 0.25,
#'   p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
#'   p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   window = 1,
#'   accrual_rate = 3,
#'   n_earlystop = 9,
#'   stay_on_1_of_3 = TRUE,
#'   accrual = "uniform",
#'   no_slot = "leave",
#'   n_trials = 500,
#'   seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{bf_boin_simulate}}, \code{\link{sim_bf_boin_multi}},
#'   \code{\link{sim_be_boin}}, \code{\link{sim_boin}}
#'
#' @export
sim_bf_boin <- function(target, p_true, p_resp, n_cohort, cohort_size,
                        window, accrual_rate, n_cap = 12,
                        backfill_dose = c("highest", "lowest"),
                        conflict_dose = c("highest", "lowest"),
                        no_slot = c("wait", "leave"),
                        accrual = c("exponential", "uniform", "fixed"),
                        dlt_time = c("weibull", "uniform"),
                        late_fraction = 0.5, resp_window = window,
                        resp_late_fraction = 0.5, resp_cor = 0,
                        n_trials = 10000, start_dose = 1, n_earlystop = 18,
                        p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                        extrasafe = FALSE, offset = 0.05,
                        stay_on_1_of_3 = FALSE, bound_mtd = FALSE,
                        mtd_max_estimate = NULL, min_mtd_sample = 1,
                        overdose_cutoff = NULL,
                        n_earlystop_rule = c("with_stay", "simple"),
                        keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- bf_boin_simulate(
    target = target, p_true = p_true, p_resp = p_resp, n_cohort = n_cohort,
    cohort_size = cohort_size, window = window, accrual_rate = accrual_rate,
    n_cap = n_cap, backfill_dose = match.arg(backfill_dose),
    conflict_dose = match.arg(conflict_dose), no_slot = match.arg(no_slot),
    accrual = match.arg(accrual), dlt_time = match.arg(dlt_time),
    late_fraction = late_fraction, resp_window = resp_window,
    resp_late_fraction = resp_late_fraction, resp_cor = resp_cor,
    n_trials = n_trials, start_dose = start_dose, n_earlystop = n_earlystop,
    p_saf = p_saf, p_tox = p_tox, cutoff_eli = cutoff_eli,
    extrasafe = extrasafe, offset = offset, stay_on_1_of_3 = stay_on_1_of_3,
    n_earlystop_rule = match.arg(n_earlystop_rule), seed = seed
  )

  if (verbose) message("Selecting the MTD ...")

  out <- backfill_oc(trials, bound_mtd = bound_mtd,
                     mtd_max_estimate = mtd_max_estimate,
                     min_mtd_sample = min_mtd_sample,
                     overdose_cutoff = overdose_cutoff,
                     keep_trials = keep_trials, call_expr = call_expr)

  if (verbose) message("Done.")
  out
}
