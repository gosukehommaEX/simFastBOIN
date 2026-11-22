# Comparison of the gBOIN-ET functions of simFastBOIN with the boinet package
#
# Yamaguchi Y, Takeda K, Maruo K. boinet: An R Package for Bayesian Optimal
# Interval Design for Dose Finding Based on Both Efficacy and Toxicity
# Outcomes. The R Journal 2025; 17(2): 15-34. doi:10.32614/RJ-2025-012
#
# boinet::gboinet() draws its random numbers differently from simFastBOIN
# (event times from a Gaussian copula, arrival times, and the categories
# with rmultinom()), so the two agree in distribution, not trial by trial.
# simFastBOIN is run with boinet_compatible = TRUE, which applies the rules
# of boinet 1.6.0; the rows with the default rules of the package are INFO.
# The scenarios are those of Table S1 of Takeda et al. (2022) in which every
# dose has a probability of a toxicity event and of an efficacy event
# strictly between 0 and 1 (1, 3, 6, 7, 8 and 10): gboinet() derives the
# parameters of its event times from these probabilities and cannot run at
# 0 or 1. The scores are the equivalent toxicity and efficacy scores of the
# article, the targets 0.47 / 1.5 and 0.95 / 1.5, and the boundaries those
# of the grid search, which is the same in the two packages.
#
# For every dose, the percentage of trials selecting it as the OBD and the
# average number of patients treated there are compared; the standard error
# is that of the difference between two independent simulations. boinet
# rounds its results to one decimal, so a cell is judged PASS when the
# difference is within four standard errors plus the rounding (0.05), as in
# inst/validation/compare-with-boinet.R, and FAIL otherwise.
#
# Usage, from the package root after devtools::load_all(), with boinet
# installed:
#
#   source("inst/validation/compare-with-gboinet.R")
#   res <- run_compare_with_gboinet()
#
# The judgements are written to
# inst/validation/results/compare-with-gboinet.csv.

gboinet_scenarios <- function() {
  source_file <- "inst/reproduce/reproduce_Takeda2022.R"
  if (!exists("takeda2022_scenarios", mode = "function")) source(source_file)
  sc <- takeda2022_scenarios("table")
  keep <- c(1, 3, 6, 7, 8, 10)
  out <- sc[keep]
  names(out) <- paste("Scenario", keep)
  out
}

compare_one_with_gboinet <- function(scenario, name,
                                     obd_method = "max_resp",
                                     n_trials_boinet = 5000,
                                     n_trials = 100000, seed = 1) {
  boinet_obd <- c(max_resp = "max.effprob",
                  utility_weighted = "utility.weighted",
                  utility_truncated_linear = "utility.truncated.linear",
                  utility_scoring = "utility.scoring")[[obd_method]]
  n_doses <- ncol(scenario$p_grade)
  score <- c(0, 0.5, 1, 1.5)
  score_resp <- c(0, 1, 1.5)

  ref <- boinet::gboinet(
    n.dose = n_doses, start.dose = 1, size.cohort = 3, n.cohort = 10,
    toxprob = scenario$p_grade, effprob = scenario$p_resp_grade,
    sev.weight = score, res.weight = score_resp, phi = 0.47 / 1.5,
    delta = 0.95 / 1.5, tau.T = 30, tau.E = 30, te.corr = 0.2, accrual = 10,
    estpt.method = "obs.prob", obd.method = boinet_obd,
    n.sim = n_trials_boinet, seed.sim = seed
  )
  ref_sel <- c(unname(ref$prop.select), unname(ref$prop.stop))

  blocks <- lapply(c(TRUE, FALSE), function(compat) {
    oc <- sim_gboin_et(
      target = 0.47 / 1.5, target_resp = 0.95 / 1.5,
      p_grade = scenario$p_grade, p_resp_grade = scenario$p_resp_grade,
      score = score, score_resp = score_resp, n_cohort = 10,
      cohort_size = 3, n_trials = n_trials, resp_cor = 0.2,
      obd_method = obd_method, boinet_compatible = compat,
      keep_trials = TRUE, seed = seed
    )
    ours_sel <- c(oc$sel_percent, oc$percent_no_obd)
    p <- (ours_sel * n_trials + ref_sel * n_trials_boinet) /
      (100 * (n_trials + n_trials_boinet))
    se_sel <- 100 * sqrt(p * (1 - p) * (1 / n_trials + 1 / n_trials_boinet))
    se_pts <- apply(oc$trials$n_pts, 2, stats::sd) *
      sqrt(1 / n_trials + 1 / n_trials_boinet)
    out <- data.frame(
      design = if (compat) "gBOIN-ET" else "gBOIN-ET (rules of the package)",
      scenario = name,
      obd_method = obd_method,
      quantity = c(rep("OBD selected (%)", n_doses), "No OBD (%)",
                   rep("patients treated", n_doses)),
      dose = c(seq_len(n_doses), NA, seq_len(n_doses)),
      boinet = c(ref_sel, unname(ref$n.patient)),
      simFastBOIN = c(unname(ours_sel), unname(oc$n_pts_dose)),
      se = c(se_sel, se_pts),
      stringsAsFactors = FALSE
    )
    out$difference <- out$simFastBOIN - out$boinet
    out$z <- ifelse(out$se > 0, out$difference / out$se, NA_real_)
    out$judgement <- if (compat) {
      gboinet_judgement(out$difference, out$se, 0.05)
    } else {
      "INFO"
    }
    out
  })
  do.call(rbind, blocks)
}

# PASS when a difference lies within four standard errors plus the rounding
# of the reference.
gboinet_judgement <- function(difference, se, rounding) {
  ifelse(abs(difference) <= 4 * se + rounding + 1e-12, "PASS", "FAIL")
}

run_compare_with_gboinet <- function(n_trials_boinet = 5000,
                                     n_trials = 100000,
                                     out_dir = "inst/validation/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_gboin_et", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  if (!requireNamespace("boinet", quietly = TRUE)) {
    stop("the 'boinet' package is needed", call. = FALSE)
  }
  scenarios <- gboinet_scenarios()
  runs <- list()
  # Default selection in the six scenarios.
  for (s in seq_along(scenarios)) {
    runs[[length(runs) + 1L]] <- compare_one_with_gboinet(
      scenarios[[s]], names(scenarios)[s],
      n_trials_boinet = n_trials_boinet, n_trials = n_trials, seed = s
    )
  }
  # The utilities in two scenarios.
  for (s in c(1, 4)) {
    for (m in c("utility_weighted", "utility_truncated_linear",
                "utility_scoring")) {
      runs[[length(runs) + 1L]] <- compare_one_with_gboinet(
        scenarios[[s]], names(scenarios)[s], obd_method = m,
        n_trials_boinet = n_trials_boinet, n_trials = n_trials,
        seed = 100 + s
      )
    }
  }
  res <- do.call(rbind, runs)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "compare-with-gboinet.csv"),
                   row.names = FALSE)
  res
}
