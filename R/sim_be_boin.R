#' Operating Characteristics of a BE-BOIN Design
#'
#' @description
#'   Simulate the BOIN design with backfilling and late-onset toxicity
#'   (BE-BOIN) of Chen et al. (2026) many times under one scenario of DLT and
#'   response probabilities, and summarize how often each dose is selected as
#'   the MTD, how many patients are treated and backfilled at each dose, how
#'   many DLTs and responses occur, and how long the trials last.
#'
#' @inheritParams be_boin_simulate
#' @inheritParams sim_tite_boin
#'
#' @return
#'   An object of class \code{c("be_boin_oc", "backfill_oc", "boin_oc")}, with
#'   the components described in \code{\link{sim_bf_boin}}.
#'
#' @details
#'   The trials are simulated by \code{\link{be_boin_simulate}}, whose details
#'   describe the conduct of a trial. The MTD is selected as in
#'   \code{\link{sim_bf_boin}}.
#'
#'   When no dose is ever opened for backfilling, for example with every
#'   \code{p_resp} equal to 0, and with \code{no_slot = "wait"}, every result
#'   that the two functions share agrees with that of
#'   \code{\link{sim_tite_boin}} with \code{method = "imputation"}, the same
#'   suspension rules and the same seed.
#'
#'   Only the MTD is selected. The second stage of Chen et al. (2026), which
#'   randomizes patients between doses to choose the optimal biological dose, is
#'   not part of the simulation.
#'
#' @references
#'   Chen, K., Zhao, Y., Takeda, K. and Yuan, Y. (2026). BE-BOIN: A Dose
#'   Optimization Design Accommodating Backfill and Late-Onset Toxicity.
#'   Therapeutic Innovation and Regulatory Science.
#'   \doi{10.1007/s43441-026-00994-0}
#'
#' @examples
#' # A three month window with two patients a month
#' oc <- sim_be_boin(
#'   target = 0.25,
#'   p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
#'   p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   window = 3,
#'   accrual_rate = 2,
#'   n_earlystop = 9,
#'   n_trials = 500,
#'   seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{be_boin_simulate}}, \code{\link{sim_be_boin_multi}},
#'   \code{\link{sim_bf_boin}}, \code{\link{sim_tite_boin}}
#'
#' @export
sim_be_boin <- function(target, p_true, p_resp, n_cohort, cohort_size,
                        window, accrual_rate, n_cap = 12,
                        backfill_dose = c("highest", "lowest"),
                        conflict_dose = c("lowest", "highest"),
                        no_slot = c("wait", "leave"),
                        accrual = c("exponential", "uniform", "fixed"),
                        dlt_time = c("weibull", "uniform"),
                        late_fraction = 0.5, resp_window = window,
                        resp_late_fraction = 0.5, resp_cor = 0,
                        max_pending_ratio = 0.49, min_completed = 0,
                        min_follow_up = 0.25,
                        n_trials = 10000, start_dose = 1, n_earlystop = 18,
                        p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                        extrasafe = FALSE, offset = 0.05, bound_mtd = FALSE,
                        mtd_max_estimate = NULL, min_mtd_sample = 1,
                        overdose_cutoff = NULL,
                        n_earlystop_rule = c("with_stay", "simple"),
                        keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- be_boin_simulate(
    target = target, p_true = p_true, p_resp = p_resp, n_cohort = n_cohort,
    cohort_size = cohort_size, window = window, accrual_rate = accrual_rate,
    n_cap = n_cap, backfill_dose = match.arg(backfill_dose),
    conflict_dose = match.arg(conflict_dose), no_slot = match.arg(no_slot),
    accrual = match.arg(accrual), dlt_time = match.arg(dlt_time),
    late_fraction = late_fraction, resp_window = resp_window,
    resp_late_fraction = resp_late_fraction, resp_cor = resp_cor,
    max_pending_ratio = max_pending_ratio, min_completed = min_completed,
    min_follow_up = min_follow_up, n_trials = n_trials,
    start_dose = start_dose, n_earlystop = n_earlystop, p_saf = p_saf,
    p_tox = p_tox, cutoff_eli = cutoff_eli, extrasafe = extrasafe,
    offset = offset, n_earlystop_rule = match.arg(n_earlystop_rule),
    seed = seed
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
