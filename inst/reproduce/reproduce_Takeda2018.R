# Reproduction of Takeda, Taguri and Morita (2018) with simFastBOIN
#
# Takeda K, Taguri M, Morita S. BOIN-ET: Bayesian optimal interval design for
# dose finding based on both efficacy and toxicity outcomes. Pharmaceutical
# Statistics 2018; 17(4): 383-395. doi:10.1002/pst.1864
#
# Three results of the article are checked:
#   1. the boundaries (lambda1, lambda2, eta1) = (0.14, 0.35, 0.48) of
#      Section 2.2.1, an exact comparison;
#   2. the decision tables of Table 2 for 3 and 6 patients, an exact
#      comparison;
#   3. the BOIN-ET rows of Table 3 (12 scenarios, 1000 trials each in the
#      article), a comparison within Monte Carlo error.
#
# For Table 3 the settings of Section 3.1 are used: 12 cohorts of 3 starting
# at the lowest dose, a binary correlation of 0.2 between the DLT and the
# response of a patient (the formulas of Sato et al. 2016), and a dose removed
# when Pr(DLT rate < 0.40) < 0.10 or Pr(response rate > 0.20) < 0.10. The OD is
# the dose with the highest fractional polynomial estimate of the response
# rate up to the MTD, which needs the 'mfp' package. The article does not give
# the stopping sample size; the maximum, 36, is used. The dose rules are the
# defaults of sim_boin_et(), which keep to the statement of the article that
# doses are not skipped during escalation (see ?boin_et_simulate).
#
# The article rounds its percentages and patient numbers to whole numbers and
# simulated 1000 trials, so a cell is judged PASS when it lies within four
# standard errors (of the difference between two independent simulations)
# plus 0.5 for the rounding, and FAIL otherwise. The column "ET" of the
# article (early termination) is compared with the percentage of trials
# without an OD.
#
# The scenarios and the published values were transcribed from the images of
# Table 3, not from extracted text. Scenario 11 is as printed, with a
# response rate of 0.55 at the lowest dose.
#
# Usage, from the package root after devtools::load_all():
#
#   source("inst/reproduce/reproduce_Takeda2018.R")
#   res <- run_reproduce_takeda2018()
#
# The judgements are written to inst/reproduce/results/reproduce_Takeda2018.csv.
# The fractional polynomial fits dominate the run time of Table 3; they are
# computed once for every distinct set of final data.

takeda2018_scenarios <- function() {
  tox <- list(
    c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
    c(0.05, 0.10, 0.15, 0.20, 0.25, 0.40),
    c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30),
    c(0.05, 0.10, 0.15, 0.25, 0.45, 0.65),
    c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30),
    c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
    c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
    c(0.05, 0.08, 0.22, 0.35, 0.55, 0.75),
    c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45),
    c(0.20, 0.40, 0.50, 0.55, 0.62, 0.64),
    c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30),
    c(0.70, 0.75, 0.80, 0.85, 0.90, 0.95)
  )
  eff <- list(
    c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
    c(0.05, 0.10, 0.25, 0.40, 0.65, 0.70),
    c(0.05, 0.20, 0.35, 0.50, 0.45, 0.40),
    c(0.05, 0.15, 0.25, 0.60, 0.65, 0.70),
    c(0.05, 0.25, 0.50, 0.45, 0.40, 0.35),
    c(0.05, 0.30, 0.55, 0.57, 0.59, 0.61),
    c(0.30, 0.70, 0.60, 0.50, 0.40, 0.30),
    c(0.30, 0.68, 0.70, 0.72, 0.74, 0.76),
    c(0.65, 0.60, 0.55, 0.50, 0.45, 0.40),
    c(0.50, 0.55, 0.60, 0.62, 0.64, 0.66),
    c(0.55, 0.25, 0.50, 0.60, 0.70, 0.80),
    c(0.20, 0.30, 0.40, 0.50, 0.60, 0.70)
  )
  lapply(seq_along(tox), function(s) {
    list(name = paste("Scenario", s), p_true = tox[[s]], p_resp = eff[[s]])
  })
}

# BOIN-ET rows of Table 3: early termination (%), recommendation (%) and
# allocation (patients) at d1 to d6.
takeda2018_table3 <- function() {
  et <- c(0, 1, 1, 1, 1, 1, 0, 0, 2, 4, 1, 100)
  rec <- rbind(
    c(0, 12, 56, 25, 7, 0),
    c(1, 2, 6, 26, 58, 8),
    c(0, 3, 15, 47, 22, 12),
    c(1, 2, 10, 75, 12, 0),
    c(0, 4, 45, 27, 16, 7),
    c(0, 10, 63, 21, 4, 1),
    c(1, 82, 13, 3, 0, 0),
    c(1, 76, 19, 4, 0, 0),
    c(72, 17, 5, 2, 1, 0),
    c(64, 28, 4, 1, 0, 0),
    c(0, 2, 33, 42, 19, 4),
    c(0, 0, 0, 0, 0, 0)
  )
  alloc <- rbind(
    c(4, 9, 15, 6, 2, 0),
    c(3, 5, 6, 9, 11, 2),
    c(3, 5, 7, 11, 5, 4),
    c(3, 5, 7, 16, 4, 1),
    c(4, 5, 13, 7, 4, 3),
    c(4, 8, 16, 6, 2, 0),
    c(5, 26, 4, 1, 0, 0),
    c(5, 25, 5, 1, 0, 0),
    c(27, 6, 2, 1, 0, 0),
    c(21, 10, 3, 1, 0, 0),
    c(3, 5, 13, 10, 4, 1),
    c(4, 1, 0, 0, 0, 0)
  )
  list(et = et, rec = rec, alloc = alloc)
}

reproduce_takeda2018_boundary <- function() {
  bd <- boin_et_boundary(target = 0.30, target_resp = 0.60)
  data.frame(
    item = "Boundaries (Section 2.2.1)",
    quantity = c("lambda1", "lambda2", "eta1"),
    dose = NA_integer_,
    published = c(0.14, 0.35, 0.48),
    reproduced = c(bd$lambda1, bd$lambda2, bd$eta1),
    tolerance = 1e-12,
    stringsAsFactors = FALSE
  )
}

reproduce_takeda2018_table2 <- function() {
  table_3 <- rbind(
    c("E", "E", "S", "S"),
    c("E/S/D", "E/S/D", "S", "S"),
    c("D", "D", "D", "D"),
    c("D", "D", "D", "D")
  )
  esd <- c("E/S/D", "E/S/D", "E/S/D", "S", "S", "S", "S")
  table_6 <- rbind(c("E", "E", "E", "S", "S", "S", "S"), esd, esd,
                   rep("D", 7), rep("D", 7), rep("D", 7), rep("D", 7))
  tab <- boin_et_decision_table(target = 0.30, target_resp = 0.60, n = c(3, 6))
  rows <- list()
  for (k in c(3, 6)) {
    published <- if (k == 3) table_3 else table_6
    part <- tab[tab$n == k, ]
    rows[[length(rows) + 1L]] <- data.frame(
      item = paste0("Table 2 (", k, " patients)"),
      quantity = paste0(part$n_tox, " DLT, ", part$n_resp, " responses"),
      dose = NA_integer_,
      published = published[cbind(part$n_tox + 1, part$n_resp + 1)],
      reproduced = part$decision,
      tolerance = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

reproduce_takeda2018_table3 <- function(n_trials = 2000, seed = 2018,
                                        est_method = "fp_logistic") {
  scenarios <- takeda2018_scenarios()
  pub <- takeda2018_table3()
  rows <- list()
  for (s in seq_along(scenarios)) {
    oc <- sim_boin_et(
      target = 0.30, target_resp = 0.60, p_true = scenarios[[s]]$p_true,
      p_resp = scenarios[[s]]$p_resp, n_cohort = 12, cohort_size = 3,
      n_trials = n_trials, cutoff_eli = 0.90, cutoff_eli_resp = 0.90,
      p_eli = 0.40, p_eli_resp = 0.20, resp_cor = 0.2,
      resp_cor_type = "binary", est_method = est_method,
      obd_method = "max_resp", keep_trials = TRUE, seed = seed + s
    )
    # Standard errors of the difference between the article (1000 trials)
    # and this run.
    p_rec <- c(oc$percent_no_obd, oc$sel_percent) / 100
    se_rec <- 100 * sqrt(p_rec * (1 - p_rec) * (1 / 1000 + 1 / n_trials))
    sd_alloc <- apply(oc$trials$n_pts, 2, stats::sd)
    se_alloc <- sd_alloc * sqrt(1 / 1000 + 1 / n_trials)
    rows[[length(rows) + 1L]] <- data.frame(
      item = scenarios[[s]]$name,
      quantity = c("ET (%)", rep("recommendation (%)", 6),
                   rep("allocation", 6)),
      dose = c(NA, 1:6, 1:6),
      published = c(pub$et[s], pub$rec[s, ], pub$alloc[s, ]),
      reproduced = c(oc$percent_no_obd, unname(oc$sel_percent),
                     unname(oc$n_pts_dose)),
      tolerance = c(4 * se_rec + 0.5, 4 * se_alloc + 0.5),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

run_reproduce_takeda2018 <- function(n_trials = 2000, seed = 2018,
                                     out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_boin_et", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  bd <- reproduce_takeda2018_boundary()
  bd$judgement <- ifelse(abs(bd$published - bd$reproduced) <= bd$tolerance,
                         "PASS", "FAIL")
  bd$published <- as.character(bd$published)
  bd$reproduced <- as.character(bd$reproduced)
  t2 <- reproduce_takeda2018_table2()
  t2$judgement <- ifelse(t2$published == t2$reproduced, "PASS", "FAIL")
  exact <- rbind(bd, t2)

  t3 <- reproduce_takeda2018_table3(n_trials = n_trials, seed = seed)
  t3$judgement <- ifelse(abs(t3$published - t3$reproduced) <= t3$tolerance,
                         "PASS", "FAIL")
  t3$reproduced <- round(t3$reproduced, 2)
  t3$tolerance <- round(t3$tolerance, 2)
  t3$published <- as.character(t3$published)
  t3$reproduced <- as.character(t3$reproduced)

  res <- rbind(exact, t3)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Takeda2018.csv"),
                   row.names = FALSE)
  res
}
