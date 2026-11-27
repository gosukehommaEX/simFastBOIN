# Comparison of the BOIN-ET and TITE-BOIN-ET functions of simFastBOIN with the
# boinet package
#
# Yamaguchi Y, Takeda K, Maruo K. boinet: An R Package for Bayesian Optimal
# Interval Design for Dose Finding Based on Both Efficacy and Toxicity
# Outcomes. The R Journal 2025; 17(2): 15-34. doi:10.32614/RJ-2025-012
#
# boinet::boinet() and boinet::tite.boinet() draw their random numbers
# differently from simFastBOIN (event times from a Gaussian copula and arrival
# times), so the two agree in distribution, not trial by trial. simFastBOIN is
# run with boinet_compatible = TRUE, which applies the rules of boinet 1.6.0.
# For every dose, the percentage of trials selecting it as the OBD and the
# average number of patients treated there are compared with a z statistic,
# the standard error being that of the difference between two independent
# simulations; for TITE-BOIN-ET the average trial duration is compared as
# well. boinet rounds its results to one decimal, so a cell is judged PASS
# when the difference is within four standard errors plus the rounding (0.05),
# as in inst/reproduce/reproduce_Takeda2018.R, and FAIL otherwise.
# rejudge_compare_with_boinet() applies the rule again to a saved CSV.
# The decision rules, the admissibility rules and the OBD selection themselves
# are compared exactly in tests/testthat/test-boin_et_boundary.R and
# test-boin_et_select_obd.R.
#
# Usage, from the package root after devtools::load_all(), with boinet
# installed:
#
#   source("inst/validation/compare-with-boinet.R")
#   res <- run_compare_with_boinet()
#
# The judgements are written to
# inst/validation/results/compare-with-boinet.csv.

boinet_scenarios <- function() {
  source_file <- "inst/reproduce/reproduce_Takeda2018.R"
  if (!exists("takeda2018_scenarios", mode = "function")) source(source_file)
  # The scenarios as printed in Table 3 of Takeda et al. (2018), with 0.55 at
  # the lowest dose of scenario 11; both packages simulate the same ones.
  takeda2018_scenarios("printed")
}

compare_one_with_boinet <- function(scenario, est_method = "observed",
                                    obd_method = "max_resp",
                                    n_trials_boinet = 5000,
                                    n_trials = 100000, seed = 1,
                                    design = c("BOIN-ET", "TITE-BOIN-ET"),
                                    tau_t = 30, tau_e = 45, accrual = 10,
                                    event_time = "weibull",
                                    enroll = "uniform") {
  design <- match.arg(design)
  boinet_est <- c(observed = "obs.prob", multi_iso = "multi.iso",
                  fp_logistic = "fp.logistic")[[est_method]]
  boinet_obd <- c(max_resp = "max.effprob",
                  utility_weighted = "utility.weighted",
                  utility_truncated_linear = "utility.truncated.linear",
                  utility_scoring = "utility.scoring")[[obd_method]]
  n_doses <- length(scenario$p_true)

  if (design == "BOIN-ET") {
    ref <- boinet::boinet(
      n.dose = n_doses, start.dose = 1, size.cohort = 3, n.cohort = 12,
      toxprob = scenario$p_true, effprob = scenario$p_resp, phi = 0.30,
      delta = 0.60, tau.T = 30, tau.E = 30, te.corr = 0.2, accrual = 10,
      estpt.method = boinet_est, obd.method = boinet_obd,
      n.sim = n_trials_boinet, seed.sim = seed
    )
    oc <- sim_boin_et(
      target = 0.30, target_resp = 0.60, p_true = scenario$p_true,
      p_resp = scenario$p_resp, n_cohort = 12, cohort_size = 3,
      n_trials = n_trials, resp_cor = 0.2, est_method = est_method,
      obd_method = obd_method, boinet_compatible = TRUE, keep_trials = TRUE,
      seed = seed
    )
  } else {
    ref <- boinet::tite.boinet(
      n.dose = n_doses, start.dose = 1, size.cohort = 3, n.cohort = 12,
      toxprob = scenario$p_true, effprob = scenario$p_resp, phi = 0.30,
      delta = 0.60, tau.T = tau_t, tau.E = tau_e, te.corr = 0.2,
      gen.event.time = event_time, accrual = accrual,
      gen.enroll.time = enroll, estpt.method = boinet_est,
      obd.method = boinet_obd, n.sim = n_trials_boinet, seed.sim = seed
    )
    oc <- sim_tite_boin_et(
      target = 0.30, target_resp = 0.60, p_true = scenario$p_true,
      p_resp = scenario$p_resp, n_cohort = 12, cohort_size = 3,
      window = tau_t, resp_window = tau_e, accrual_rate = 1 / accrual,
      accrual = enroll, event_time = event_time, n_trials = n_trials,
      resp_cor = 0.2, est_method = est_method, obd_method = obd_method,
      boinet_compatible = TRUE, keep_trials = TRUE, seed = seed
    )
  }

  ours_sel <- c(oc$sel_percent, oc$percent_no_obd)
  ref_sel <- c(unname(ref$prop.select), unname(ref$prop.stop))
  p <- (ours_sel * n_trials + ref_sel * n_trials_boinet) /
    (100 * (n_trials + n_trials_boinet))
  se_sel <- 100 * sqrt(p * (1 - p) * (1 / n_trials + 1 / n_trials_boinet))
  se_pts <- apply(oc$trials$n_pts, 2, stats::sd) *
    sqrt(1 / n_trials + 1 / n_trials_boinet)

  out <- data.frame(
    design = design,
    scenario = scenario$name,
    est_method = est_method,
    obd_method = obd_method,
    quantity = c(rep("OBD selected (%)", n_doses), "No OBD (%)",
                 rep("patients treated", n_doses)),
    dose = c(seq_len(n_doses), NA, seq_len(n_doses)),
    boinet = c(ref_sel, unname(ref$n.patient)),
    simFastBOIN = c(unname(ours_sel), unname(oc$n_pts_dose)),
    se = c(se_sel, se_pts),
    stringsAsFactors = FALSE
  )
  if (design == "TITE-BOIN-ET") {
    out <- rbind(out, data.frame(
      design = design, scenario = scenario$name, est_method = est_method,
      obd_method = obd_method, quantity = "trial duration", dose = NA,
      boinet = unname(ref$duration), simFastBOIN = oc$duration,
      se = stats::sd(oc$trials$duration) *
        sqrt(1 / n_trials + 1 / n_trials_boinet),
      stringsAsFactors = FALSE
    ))
  }
  out$difference <- out$simFastBOIN - out$boinet
  out$z <- ifelse(out$se > 0, out$difference / out$se, NA_real_)
  out$judgement <- compare_judgement(out$difference, out$se, 0.05)
  out
}

# PASS when a difference lies within four standard errors plus the rounding
# of the reference.
compare_judgement <- function(difference, se, rounding) {
  ifelse(abs(difference) <= 4 * se + rounding + 1e-12, "PASS", "FAIL")
}

rejudge_compare_with_boinet <- function(
    path = "inst/validation/results/compare-with-boinet.csv") {
  res <- utils::read.csv(path, stringsAsFactors = FALSE)
  res$judgement <- compare_judgement(res$difference, res$se, 0.05)
  utils::write.csv(res, path, row.names = FALSE)
  res
}

run_compare_with_boinet <- function(n_trials_boinet = 5000,
                                    n_trials = 100000,
                                    n_trials_tite_boinet = 2000,
                                    out_dir = "inst/validation/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_tite_boin_et", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  if (!requireNamespace("boinet", quietly = TRUE)) {
    stop("the 'boinet' package is needed", call. = FALSE)
  }
  scenarios <- boinet_scenarios()
  runs <- list()
  # Default selection in every scenario of Takeda et al. (2018).
  for (s in seq_along(scenarios)) {
    runs[[length(runs) + 1L]] <- compare_one_with_boinet(
      scenarios[[s]], n_trials_boinet = n_trials_boinet, n_trials = n_trials,
      seed = s
    )
  }
  # The other estimators and selection rules in three scenarios.
  for (s in c(1, 7, 9)) {
    for (m in c("utility_weighted", "utility_truncated_linear",
                "utility_scoring")) {
      runs[[length(runs) + 1L]] <- compare_one_with_boinet(
        scenarios[[s]], obd_method = m, n_trials_boinet = n_trials_boinet,
        n_trials = n_trials, seed = 100 + s
      )
    }
    runs[[length(runs) + 1L]] <- compare_one_with_boinet(
      scenarios[[s]], est_method = "multi_iso",
      n_trials_boinet = n_trials_boinet, n_trials = n_trials, seed = 200 + s
    )
  }
  # TITE-BOIN-ET: a 30-day DLT window, a 45-day response window and one
  # patient every 10 days, in every scenario; then uniform event times,
  # exponential arrivals and the utility selection in two scenarios.
  for (s in seq_along(scenarios)) {
    runs[[length(runs) + 1L]] <- compare_one_with_boinet(
      scenarios[[s]], design = "TITE-BOIN-ET",
      n_trials_boinet = n_trials_tite_boinet, n_trials = n_trials,
      seed = 300 + s
    )
  }
  for (s in c(1, 7)) {
    runs[[length(runs) + 1L]] <- compare_one_with_boinet(
      scenarios[[s]], design = "TITE-BOIN-ET", event_time = "uniform",
      enroll = "exponential", obd_method = "utility_scoring", tau_e = 60,
      accrual = 5, n_trials_boinet = n_trials_tite_boinet,
      n_trials = n_trials, seed = 400 + s
    )
  }
  res <- do.call(rbind, runs)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "compare-with-boinet.csv"),
                   row.names = FALSE)
  res
}
