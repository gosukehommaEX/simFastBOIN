#' Simulate Trials of the TITE-Comb-BOIN12 Design
#'
#' @description
#'   Simulate trials of the TITE-Comb-BOIN12 design of Lu et al. (2025) for
#'   two drugs given together, in which the DLT and the response may be
#'   observed late: the rules of \code{\link{comb_boin12_simulate}} applied at
#'   the decision times and with the interim estimates of
#'   \code{\link{tite_boin12_simulate}}. Only the trial conduct is simulated;
#'   the OBD combination of each trial is selected by
#'   \code{\link{comb_boin12_select_obd}} (both are done by
#'   \code{\link{sim_tite_comb_boin12}}).
#'
#' @inheritParams comb_boin12_simulate
#' @inheritParams tite_boin12_simulate
#'
#' @param pending_use
#'   Character scalar. Where the pending outcomes count: \code{"all"} (the
#'   default), in the estimated DLT rate of the decision, in the
#'   admissibility and in the utility, as in
#'   \code{\link{tite_boin12_simulate}}; \code{"utility"}, only in the
#'   quasi-events of the utility, as Sections 2.2 and 3.1 of Lu et al.
#'   (2025) state the rules: the decision uses the observed DLT rate (DLTs
#'   over patients treated) and the admissibility the observed numbers of
#'   DLTs and responses, a pending outcome counting as no event.
#'
#' @return
#'   An object of class \code{comb_boin12_trials}, as returned by
#'   \code{\link{comb_boin12_simulate}}, with \code{design} set to
#'   \code{"TITE-Comb-BOIN12"}, the counts of the final outcomes of all
#'   patients, and one more component
#'   \item{duration}{Numeric vector of the trial durations, in the unit of
#'     \code{window}: the time of the last decision of a trial stopped
#'     because no combination was admissible, and otherwise the time at
#'     which every patient has finished both assessments.}
#'   The timing arguments are added to \code{settings}.
#'
#' @details
#'   Patients arrive, are followed and are counted as in
#'   \code{\link{tite_boin12_simulate}}: the first patient of a cohort
#'   arrives at the time of the previous decision; the decision waits while
#'   fewer than a share \code{min_completed_share} of the patients at the
#'   current combination have finished both assessments; a pending outcome
#'   counts in the estimates of the DLT and response rates, in the
#'   admissibility and in the quasi-events of the desirability by the rules
#'   of Zhou et al. (2022). The decision follows the rules of
#'   \code{\link{comb_boin12_simulate}}, with the estimated DLT rate of the
#'   current combination, the admissibility computed from the interim
#'   posteriors (the DLT rule of a combination also excluding every
#'   combination with both levels at least its own, and \code{extrasafe}
#'   lowering the cutoff at the lowest combination) and the patients treated.
#'
#'   With \code{pending_use = "utility"} the decision uses the observed DLT
#'   rate of the current combination and the admissibility the posteriors
#'   Beta(1 + DLTs, 1 + patients treated - DLTs) and the same for the
#'   responses; only the desirability uses the pending outcomes.
#'
#'   With every outcome observed at every decision these are the quantities
#'   of Comb-BOIN12: with \code{accrual = "fixed"} and
#'   \code{min_completed_share = 1} no variate is drawn for the times and the
#'   trials are exactly those of \code{\link{comb_boin12_simulate}} with the
#'   same seed. With one row or one column of combinations the trials are
#'   those of \code{\link{tite_boin12_simulate}} without exploration.
#'
#'   Lu et al. (2025) state the time-to-event rules of TITE-BOIN12 for
#'   TITE-Comb-BOIN12 without further detail; the accrual suspension and the
#'   estimates are those of \code{\link{tite_boin12_simulate}}, as in
#'   \code{inst/reproduce/reproduce_Lu2025.R}.
#'
#' @references
#'   Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
#'   utility-based Bayesian optimal interval design for dose optimization in
#'   cancer drug-combination trials. Statistics in Biopharmaceutical
#'   Research, 17(2), 266-276. \doi{10.1080/19466315.2024.2370403}
#'
#'   Zhou, Y., Lin, R., Lee, J. J., Li, D., Wang, L., Li, R. and Yuan, Y.
#'   (2022). TITE-BOIN12: A Bayesian phase I/II trial design to find the
#'   optimal biological dose with late-onset toxicity and efficacy. Statistics
#'   in Medicine, 41(11), 1918-1931. \doi{10.1002/sim.9337}
#'
#' @examples
#' # Scenario 16 of Lu et al. (2025), with their utility (wT = 0.11), windows
#' # of 1 and 2 months and two patients a month
#' p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
#' p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
#' trials <- tite_comb_boin12_simulate(
#'   target = 0.33, target_resp = 0.2, p_true = p_true, p_resp = p_resp,
#'   n_cohort = 14, cohort_size = 3, window = 1, resp_window = 2,
#'   accrual_rate = 2, event_time = "uniform", n_trials = 200,
#'   psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11, benchmark_weight = 0.6,
#'   extrasafe = TRUE, offset = 0.10, seed = 123
#' )
#' trials
#'
#' @seealso \code{\link{sim_tite_comb_boin12}},
#'   \code{\link{comb_boin12_simulate}}, \code{\link{tite_boin12_simulate}}
#'
#' @export
tite_comb_boin12_simulate <- function(target, target_resp, p_true, p_resp,
                                      n_cohort, cohort_size, window,
                                      accrual_rate, resp_window = window,
                                      accrual = c("uniform", "exponential",
                                                  "fixed"),
                                      event_time = c("weibull", "uniform"),
                                      late_fraction = 0.5,
                                      resp_late_fraction = 0.5,
                                      min_completed_share = 0.5,
                                      pending_estimate = c("conditional",
                                                           "marginal"),
                                      pending_use = c("all", "utility"),
                                      n_trials = 10000, start_dose = c(1, 1),
                                      n_earlystop = NULL,
                                      n_earlystop_rule = c("with_stay",
                                                           "simple"),
                                      p_saf = NULL, p_tox = NULL,
                                      cutoff_eli = 0.95,
                                      cutoff_eli_resp = 0.90,
                                      extrasafe = FALSE, offset = 0.05,
                                      psi00 = 40, psi11 = 60,
                                      benchmark_weight = 0.5,
                                      utility_prior = c(1, 1), n_star = 6,
                                      tie_break = c("higher", "random"),
                                      no_candidate = c("nearest", "stop"),
                                      resp_cor = 0,
                                      resp_cor_type = c("copula", "binary"),
                                      seed = 123) {

  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  pending_estimate <- match.arg(pending_estimate)
  pending_use <- match.arg(pending_use)
  n_earlystop_rule <- match.arg(n_earlystop_rule)
  tie_break <- match.arg(tie_break)
  no_candidate <- match.arg(no_candidate)
  resp_cor_type <- match.arg(resp_cor_type)

  check_p_true_matrix(p_true)
  check_p_resp_matrix(p_resp, p_true)
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_start_combination(start_dose, p_true)
  check_count(n_star, "n_star", 1L)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  check_flag(extrasafe, "extrasafe")
  check_comb_offset(offset, cutoff_eli)
  check_scalar_prob(resp_cor, "resp_cor", lower = -1, upper = 1)
  if (!is.numeric(min_completed_share) || length(min_completed_share) != 1L ||
      !is.finite(min_completed_share) || min_completed_share < 0 ||
      min_completed_share > 1) {
    stop("'min_completed_share' must be a single number from 0 to 1",
         call. = FALSE)
  }
  boin_et_check_timing(window, resp_window, accrual_rate, late_fraction,
                       resp_late_fraction)
  par_t <- boin_et_event_params(as.vector(p_true), window, late_fraction,
                                event_time, "p_true")
  par_e <- boin_et_event_params(as.vector(p_resp), resp_window,
                                resp_late_fraction, event_time, "p_resp")
  des <- boin12_check_design(target, target_resp, p_saf = p_saf,
                             p_tox = p_tox, psi00 = psi00, psi11 = psi11,
                             benchmark_weight = benchmark_weight,
                             utility_prior = utility_prior)

  n_rows <- nrow(p_true)
  n_cols <- ncol(p_true)
  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)
  if (!is.null(n_earlystop)) check_count(n_earlystop, "n_earlystop", 1L)

  if (resp_cor_type == "binary") {
    link <- boin_et_resp_link(as.vector(p_true), as.vector(p_resp), resp_cor)
  } else {
    link <- list(q0 = rep(0, n_rows * n_cols), q1 = rep(0, n_rows * n_cols))
  }

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- tite_comb_boin12_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    p_resp = as.numeric(p_resp),
    q0 = as.numeric(link$q0),
    q1 = as.numeric(link$q1),
    cor_type = if (resp_cor_type == "binary") 1L else 0L,
    resp_cor = as.numeric(resp_cor),
    n_rows = as.integer(n_rows),
    n_cols = as.integer(n_cols),
    cohort_size = cohort_size_vec,
    start_row = as.integer(start_dose[1L]),
    start_col = as.integer(start_dose[2L]),
    n_earlystop = if (is.null(n_earlystop)) 0L else as.integer(n_earlystop),
    early_stop_simple = identical(n_earlystop_rule, "simple"),
    n_star = as.integer(n_star),
    tie_break = if (tie_break == "random") 1L else 0L,
    no_candidate = if (no_candidate == "stop") 1L else 0L,
    lambda_e = des$lambda_e,
    lambda_d = des$lambda_d,
    psi00 = as.numeric(psi00),
    psi11 = as.numeric(psi11),
    alpha = des$utility_prior[1L],
    beta = des$utility_prior[2L],
    u_benchmark = des$u_benchmark / 100,
    target = as.numeric(target),
    target_resp = as.numeric(target_resp),
    cutoff_eli = as.numeric(cutoff_eli),
    cutoff_lowest = as.numeric(if (extrasafe) cutoff_eli - offset else
      cutoff_eli),
    cutoff_eli_resp = as.numeric(cutoff_eli_resp),
    pending = if (pending_estimate == "marginal") 1L else 0L,
    pending_use = if (pending_use == "utility") 1L else 0L,
    window = as.numeric(window),
    resp_window = as.numeric(resp_window),
    accrual_rate = as.numeric(accrual_rate),
    arrival = match(accrual, c("uniform", "exponential", "fixed")) - 1L,
    event_time = match(event_time, c("weibull", "uniform")) - 1L,
    shape_t = par_t$shape,
    scale_t = par_t$scale,
    shape_e = par_e$shape,
    scale_e = par_e$scale,
    min_completed_share = as.numeric(min_completed_share)
  )

  trial_names <- list(NULL, paste0("A", seq_len(n_rows)),
                      paste0("B", seq_len(n_cols)))
  for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
    dimnames(res[[nm]]) <- trial_names
  }

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$n_tox,
      n_resp = res$n_resp,
      n_tox_resp = res$n_tox_resp,
      cohorts_used = res$cohorts_used,
      stop_reason = boin12_stop_reasons()[res$stop_code + 1L],
      early_stop = res$stop_code %in% c(0L, 3L),
      duration = res$duration,
      design = "TITE-Comb-BOIN12",
      boundary = list(lambda_e = des$lambda_e, lambda_d = des$lambda_d,
                      u_low = des$u_low, u_benchmark = des$u_benchmark),
      settings = list(
        target = target,
        target_resp = target_resp,
        p_true = p_true,
        p_resp = p_resp,
        n_cohort = n_cohort,
        cohort_size = cohort_size_vec,
        window = window,
        resp_window = resp_window,
        accrual_rate = accrual_rate,
        accrual = accrual,
        event_time = event_time,
        late_fraction = late_fraction,
        resp_late_fraction = resp_late_fraction,
        min_completed_share = min_completed_share,
        pending_estimate = pending_estimate,
        pending_use = pending_use,
        n_trials = n_trials,
        start_dose = as.integer(start_dose),
        n_earlystop = n_earlystop,
        n_earlystop_rule = n_earlystop_rule,
        p_saf = des$p_saf,
        p_tox = des$p_tox,
        cutoff_eli = cutoff_eli,
        cutoff_eli_resp = cutoff_eli_resp,
        extrasafe = extrasafe,
        offset = offset,
        psi00 = psi00,
        psi11 = psi11,
        benchmark_weight = benchmark_weight,
        utility_prior = des$utility_prior,
        n_star = n_star,
        tie_break = tie_break,
        no_candidate = no_candidate,
        resp_cor = resp_cor,
        resp_cor_type = resp_cor_type,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "comb_boin12_trials"
  )
}
