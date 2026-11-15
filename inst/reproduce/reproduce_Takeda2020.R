# Reproduction of Takeda, Morita and Taguri (2020) with simFastBOIN
#
# Takeda K, Morita S, Taguri M. TITE-BOIN-ET: Time-to-event Bayesian optimal
# interval design to accelerate dose-finding based on both efficacy and
# toxicity outcomes. Pharmaceutical Statistics 2020; 19: 335-349.
# doi:10.1002/pst.1995
#
# The TITE-BOIN-ET and BOIN-ET rows of Table 2 (12 scenarios, 1000 trials each
# in the article) are compared within Monte Carlo error: the early termination
# (%), the recommendation (%) and the allocation (patients) at d1 to d6, and,
# for TITE-BOIN-ET, the average trial duration in weeks.
#
# The settings of Section 3.1 are used: 12 cohorts of 3 starting at the lowest
# dose, a DLT window of 28 days, a response window of 56 days, one patient per
# 14 days, Weibull times to DLT and to response with half of the events in the
# second half of each window, and a dose removed when Pr(DLT rate < 0.40) <
# 0.10 or Pr(response rate > 0.20) < 0.10. The decisions do not wait for
# finished assessments (min_completed_share = 0): the article used no such
# rule and recommended one in its discussion (Section 5); boinet waits for
# half of the patients at the current dose. The OD is the dose with the highest
# fractional polynomial estimate of the response rate up to the MTD, which
# needs the 'mfp' package. The article does not state the correlation between
# the DLT and the response or the distribution of the times between arrivals;
# the defaults of boinet::tite.boinet() are used for TITE-BOIN-ET (a Gaussian
# copula with correlation 0.2 and uniform times between arrivals), and the
# binary correlation of 0.2 of Takeda et al. (2018) for BOIN-ET, as in
# reproduce_Takeda2018.R. The dose rules are the defaults of the package,
# which do not skip doses during escalation (see ?boin_et_simulate).
#
# The article rounds its percentages, patient numbers and durations to whole
# numbers, so a cell is judged PASS when it lies within four standard errors
# (of the difference between two independent simulations) plus 0.5 for the
# rounding, and FAIL otherwise. The column "ET" of the article (early
# termination) is compared with the percentage of trials without an OD.
#
# The durations of the article include the interval before the arrival of the
# first patient: Table 2 gives 167 or 168 weeks for BOIN-ET in scenarios 1 to
# 8 and 11, that is 12 cohorts of three intervals of 2 weeks and a response
# window of 8 weeks. The durations of sim_tite_boin_et() start at the arrival
# of the first patient, so the mean interval (14 days) is added before the
# comparison, and its variance (uniform times between arrivals) to that of the
# duration.
#
# The package, like boinet, removes a dose using every enrolled patient, the
# pending ones counted as without DLT and without response. The article does
# not say which patients this rule uses; in scenarios 9, 10 and 12 its results
# agree better with a rule that uses only the patients whose outcomes are
# known.
#
# The scenarios and the published values were extracted from the text of the
# article (pages 343-345). Scenario 11 is as printed, with a response rate of
# 0.55 at the lowest dose, as in Table 3 of Takeda et al. (2018).
#
# Usage, from the package root:
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Takeda2020.R")
#   res <- run_reproduce_takeda2020()
#
# The judgements are written to inst/reproduce/results/reproduce_Takeda2020.csv.

takeda2020_scenarios <- function() {
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

# Rows of Table 2 for one design: early termination (%), recommendation (%)
# and allocation (patients) at d1 to d6, and the duration in weeks.
takeda2020_table2 <- function(design = c("TITE-BOIN-ET", "BOIN-ET")) {
  design <- match.arg(design)
  if (design == "TITE-BOIN-ET") {
    et <- c(1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 1, 100)
    rec <- rbind(
      c(2, 19, 50, 23, 6, 1),
      c(3, 4, 13, 31, 41, 7),
      c(2, 6, 23, 44, 18, 7),
      c(2, 4, 13, 69, 11, 0),
      c(2, 9, 47, 25, 13, 4),
      c(2, 16, 61, 18, 3, 0),
      c(1, 76, 18, 4, 1, 0),
      c(1, 63, 29, 7, 0, 0),
      c(67, 23, 4, 2, 1, 0),
      c(67, 27, 2, 0, 0, 0),
      c(2, 3, 34, 37, 19, 6),
      c(0, 0, 0, 0, 0, 0)
    )
    alloc <- rbind(
      c(5, 11, 12, 5, 1, 0),
      c(5, 6, 8, 8, 7, 2),
      c(5, 6, 9, 9, 5, 2),
      c(5, 6, 8, 13, 3, 0),
      c(5, 7, 11, 6, 4, 2),
      c(5, 10, 14, 5, 1, 0),
      c(6, 22, 5, 2, 1, 0),
      c(5, 21, 8, 2, 0, 0),
      c(24, 8, 2, 1, 0, 0),
      c(22, 11, 2, 0, 0, 0),
      c(4, 6, 12, 8, 4, 1),
      c(6, 0, 0, 0, 0, 0)
    )
    dur <- c(79, 79, 79, 79, 79, 79, 79, 80, 77, 76, 79, 14)
  } else {
    et <- c(1, 1, 1, 1, 1, 1, 0, 0, 1, 2, 1, 100)
    rec <- rbind(
      c(0, 13, 54, 27, 5, 0),
      c(1, 2, 6, 25, 59, 7),
      c(1, 4, 15, 46, 24, 9),
      c(1, 2, 9, 77, 10, 0),
      c(1, 5, 47, 27, 15, 6),
      c(0, 10, 63, 23, 3, 0),
      c(1, 82, 14, 3, 1, 0),
      c(0, 76, 20, 3, 0, 0),
      c(74, 18, 5, 2, 0, 0),
      c(67, 28, 3, 0, 0, 0),
      c(0, 2, 34, 42, 18, 3),
      c(0, 0, 0, 0, 0, 0)
    )
    alloc <- rbind(
      c(4, 10, 15, 6, 1, 0),
      c(3, 5, 6, 8, 11, 2),
      c(3, 5, 7, 11, 6, 3),
      c(3, 5, 7, 17, 4, 1),
      c(4, 5, 13, 7, 4, 3),
      c(4, 8, 16, 6, 2, 0),
      c(5, 25, 4, 1, 0, 0),
      c(5, 25, 6, 1, 0, 0),
      c(27, 6, 2, 1, 0, 0),
      c(22, 11, 2, 0, 0, 0),
      c(3, 5, 13, 10, 4, 1),
      c(6, 1, 0, 0, 0, 0)
    )
    dur <- NULL
  }
  list(et = et, rec = rec, alloc = alloc, dur = dur)
}

reproduce_takeda2020_table2 <- function(n_trials = 2000, seed = 2020,
                                        est_method = "fp_logistic") {
  scenarios <- takeda2020_scenarios()
  rows <- list()
  for (design in c("TITE-BOIN-ET", "BOIN-ET")) {
    pub <- takeda2020_table2(design)
    for (s in seq_along(scenarios)) {
      common <- list(
        target = 0.30, target_resp = 0.60, p_true = scenarios[[s]]$p_true,
        p_resp = scenarios[[s]]$p_resp, n_cohort = 12, cohort_size = 3,
        n_trials = n_trials, cutoff_eli = 0.90, cutoff_eli_resp = 0.90,
        p_eli = 0.40, p_eli_resp = 0.20, resp_cor = 0.2,
        est_method = est_method, obd_method = "max_resp", keep_trials = TRUE,
        seed = seed + s
      )
      if (design == "TITE-BOIN-ET") {
        oc <- do.call(sim_tite_boin_et, c(common, list(
          window = 28, resp_window = 56, accrual_rate = 1 / 14,
          min_completed_share = 0
        )))
      } else {
        oc <- do.call(sim_boin_et, c(common, list(resp_cor_type = "binary")))
      }
      # Standard errors of the difference between the article (1000 trials)
      # and this run.
      p_rec <- c(oc$percent_no_obd, oc$sel_percent) / 100
      se_rec <- 100 * sqrt(p_rec * (1 - p_rec) * (1 / 1000 + 1 / n_trials))
      sd_alloc <- apply(oc$trials$n_pts, 2, stats::sd)
      se_alloc <- sd_alloc * sqrt(1 / 1000 + 1 / n_trials)
      out <- data.frame(
        item = scenarios[[s]]$name,
        design = design,
        quantity = c("ET (%)", rep("recommendation (%)", 6),
                     rep("allocation", 6)),
        dose = c(NA, 1:6, 1:6),
        published = c(pub$et[s], pub$rec[s, ], pub$alloc[s, ]),
        reproduced = c(oc$percent_no_obd, unname(oc$sel_percent),
                       unname(oc$n_pts_dose)),
        tolerance = c(4 * se_rec + 0.5, 4 * se_alloc + 0.5),
        stringsAsFactors = FALSE
      )
      if (design == "TITE-BOIN-ET") {
        # The interval before the first patient (see the header): uniform
        # between 0 and 28 days, with mean 14 days and variance 28^2 / 12.
        gap_mean <- 14
        gap_var <- 28^2 / 12
        weeks_mean <- (mean(oc$trials$duration) + gap_mean) / 7
        weeks_sd <- sqrt(stats::var(oc$trials$duration) + gap_var) / 7
        out <- rbind(out, data.frame(
          item = scenarios[[s]]$name, design = design,
          quantity = "duration (weeks)", dose = NA, published = pub$dur[s],
          reproduced = weeks_mean,
          tolerance = 4 * weeks_sd * sqrt(1 / 1000 + 1 / n_trials) + 0.5,
          stringsAsFactors = FALSE
        ))
      }
      rows[[length(rows) + 1L]] <- out
    }
  }
  do.call(rbind, rows)
}

run_reproduce_takeda2020 <- function(n_trials = 2000, seed = 2020,
                                     out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_tite_boin_et", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  if (!requireNamespace("mfp", quietly = TRUE)) {
    stop("the 'mfp' package is needed", call. = FALSE)
  }
  res <- reproduce_takeda2020_table2(n_trials = n_trials, seed = seed)
  res$judgement <- ifelse(abs(res$published - res$reproduced) <= res$tolerance,
                          "PASS", "FAIL")
  res$reproduced <- round(res$reproduced, 2)
  res$tolerance <- round(res$tolerance, 2)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Takeda2020.csv"),
                   row.names = FALSE)
  res
}
