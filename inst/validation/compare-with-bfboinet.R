# Comparison of sim_bf_boin_et() with the bfboinet package
#
# Takeda, K., Zhu, J. and Hirakawa, A. (2025). BF-BOIN-ET: A backfill Bayesian
# optimal interval design using efficacy and toxicity outcomes for dose
# optimization. Pharmaceutical Statistics, 24, e2470. doi:10.1002/pst.2470
#
# bfboinet::get.oc.backboinet() draws its random numbers differently from
# simFastBOIN (event times from a Gaussian copula, one seed per trial), so the
# two agree in distribution, not trial by trial. For every dose, the
# percentage of trials selecting it as the OBD, the average number of patients
# treated there and the average number backfilled there are compared with a z
# statistic, the standard error being that of the difference between two
# independent simulations; the values of bfboinet are computed from the trial
# data it returns. The average trial duration, which bfboinet rounds to one
# decimal, is compared as well. A cell is judged PASS when |z| <= 4 (or, for
# the duration, the difference is within the rounding of 0.05), and FAIL
# otherwise. The trials start at the lowest dose, where bfboinet reports the
# OBD as a dose (see ?sim_bf_boin_et).
#
# Usage, from the package root after devtools::load_all(), with bfboinet
# installed:
#
#   source("inst/validation/compare-with-bfboinet.R")
#   res <- run_compare_with_bfboinet()
#
# The judgements are written to
# inst/validation/results/compare-with-bfboinet.csv.

bfboinet_scenarios <- function() {
  list(
    list(name = "Example of bfboinet", target = 0.30, target_resp = 0.25,
         p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
         p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70)),
    list(name = "Plateau", target = 0.30, target_resp = 0.25,
         p_true = c(0.05, 0.10, 0.15, 0.30, 0.45),
         p_resp = c(0.10, 0.30, 0.45, 0.50, 0.52)),
    list(name = "Toxic", target = 0.25, target_resp = 0.30,
         p_true = c(0.30, 0.45, 0.60, 0.70),
         p_resp = c(0.20, 0.30, 0.40, 0.50)),
    list(name = "Futile high doses", target = 0.30, target_resp = 0.40,
         p_true = c(0.02, 0.05, 0.08, 0.10, 0.12, 0.15),
         p_resp = c(0.45, 0.50, 0.05, 0.05, 0.05, 0.05))
  )
}

compare_one_with_bfboinet <- function(scenario, n_trials_bfboinet = 1000,
                                      n_trials = 100000, seed = 1,
                                      tau = 1, resp_tau = 1, accrual_rate = 3,
                                      event_time = "weibull",
                                      enroll = "uniform", suspend = 0,
                                      label = "") {
  n_doses <- length(scenario$p_true)
  ref <- bfboinet::get.oc.backboinet(
    target_T = scenario$target, toxprob = scenario$p_true,
    target_E = scenario$target_resp, effprob = scenario$p_resp,
    n.dose = n_doses, startdose = 1, ncohort = 10, cohortsize = 3,
    pT.saf = 0.6 * scenario$target, pT.tox = 1.4 * scenario$target,
    pE.saf = 0.6 * scenario$target_resp, alpha.T1 = 0.5, alpha.E1 = 0.5,
    tau.T = tau, tau.E = resp_tau, te.corr = 0.2,
    gen.event.time = event_time, accrual = accrual_rate,
    gen.enroll.time = enroll, n.elimination = 6, stopping.npts = 12,
    suspend = suspend, stopping.prob.T = 0.95, stopping.prob.E = 0.90,
    ppsi01 = 0, ppsi00 = 40, ppsi11 = 60, ppsi10 = 100,
    n.sim = n_trials_bfboinet, seed.sim = seed
  )
  oc <- sim_bf_boin_et(
    target = scenario$target, target_resp = scenario$target_resp,
    p_true = scenario$p_true, p_resp = scenario$p_resp, n_cohort = 10,
    cohort_size = 3, window = tau, resp_window = resp_tau,
    accrual_rate = accrual_rate, accrual = enroll, event_time = event_time,
    suspend = suspend, n_trials = n_trials, keep_trials = TRUE, seed = seed
  )

  ref_obd <- as.vector(ref$obd)
  ref_sel <- c(vapply(seq_len(n_doses), function(j) mean(ref_obd == j), 1),
               mean(ref_obd == 0)) * 100
  ours_sel <- c(unname(oc$sel_percent), oc$percent_no_obd)
  p <- (ours_sel * n_trials + ref_sel * n_trials_bfboinet) /
    (100 * (n_trials + n_trials_bfboinet))
  se_sel <- 100 * sqrt(p * (1 - p) * (1 / n_trials + 1 / n_trials_bfboinet))

  ref_pts <- colMeans(ref$data.obs.n)
  ref_bf <- colMeans(ref$backfillcount)
  sd2 <- function(a, b) sqrt(stats::var(a) / n_trials +
                               stats::var(b) / n_trials_bfboinet)
  se_pts <- vapply(seq_len(n_doses), function(j) {
    sd2(oc$trials$n_pts[, j], ref$data.obs.n[, j])
  }, 1)
  se_bf <- vapply(seq_len(n_doses), function(j) {
    sd2(oc$trials$n_bf[, j], ref$backfillcount[, j])
  }, 1)
  se_dur <- stats::sd(oc$trials$duration) *
    sqrt(1 / n_trials + 1 / n_trials_bfboinet)

  out <- data.frame(
    scenario = paste0(scenario$name, label),
    quantity = c(rep("OBD selected (%)", n_doses), "No OBD (%)",
                 rep("patients treated", n_doses),
                 rep("patients backfilled", n_doses), "trial duration"),
    dose = c(seq_len(n_doses), NA, seq_len(n_doses), seq_len(n_doses), NA),
    bfboinet = c(ref_sel, unname(ref_pts), unname(ref_bf),
                 unname(ref$duration)),
    simFastBOIN = c(ours_sel, unname(oc$n_pts_dose), unname(oc$n_bf_dose),
                    oc$duration),
    se = c(se_sel, se_pts, se_bf, se_dur),
    rounding = c(rep(0, 2 * n_doses + 2), 0.05),
    stringsAsFactors = FALSE
  )
  out$difference <- out$simFastBOIN - out$bfboinet
  out$z <- ifelse(out$se > 0, out$difference / out$se, NA_real_)
  out$judgement <- ifelse(abs(out$difference) <= out$rounding + 1e-12 |
                            (!is.na(out$z) & abs(out$z) <= 4), "PASS", "FAIL")
  out$rounding <- NULL
  out
}

run_compare_with_bfboinet <- function(n_trials_bfboinet = 1000,
                                      n_trials = 100000,
                                      out_dir = "inst/validation/results") {
  scenarios <- bfboinet_scenarios()
  runs <- list()
  for (s in seq_along(scenarios)) {
    runs[[length(runs) + 1L]] <- compare_one_with_bfboinet(
      scenarios[[s]], n_trials_bfboinet = n_trials_bfboinet,
      n_trials = n_trials, seed = s
    )
  }
  # Uniform event times, exponential arrivals, a longer response window and
  # the suspension rule.
  runs[[length(runs) + 1L]] <- compare_one_with_bfboinet(
    scenarios[[2]], n_trials_bfboinet = n_trials_bfboinet,
    n_trials = n_trials, seed = 11, resp_tau = 2, accrual_rate = 2,
    event_time = "uniform", enroll = "exponential", suspend = 0.33,
    label = " (uniform, exponential, suspend 0.33)"
  )
  res <- do.call(rbind, runs)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "compare-with-bfboinet.csv"),
                   row.names = FALSE)
  res
}
