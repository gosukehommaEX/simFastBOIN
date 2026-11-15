# Reproduction of Takeda, Zhu and Hirakawa (2025) with simFastBOIN
#
# Takeda K, Zhu J, Hirakawa A. BF-BOIN-ET: A backfill Bayesian optimal interval
# design using efficacy and toxicity outcomes for dose optimization.
# Pharmaceutical Statistics 2025; 24: e2470. doi:10.1002/pst.2470
#
# The BF-BOIN-ET-C and BF-BOIN-ET-A rows of Table 1 (response window of 30
# days) and Table 2 (60 days) are compared within Monte Carlo error, for the
# 14 fixed scenarios of Table S1 (1000 trials each in the article): the
# percentage of trials selecting the OD (PCS), the percentage of patients
# treated at the OD (PCA), the percentage of trials selecting a toxic dose
# (PTS), the percentage of patients treated at toxic doses (PTA), the average
# number of patients (N) and the average trial duration in months (Dur).
#
# The settings of Section 3.1 are used: a target DLT rate of 0.30, a minimum
# response rate of 0.25, 10 cohorts of 3 in the dose escalation starting at
# the lowest dose, at most 12 patients per dose, a DLT window of 30 days, one
# patient per 10 days, Weibull times with half of the events in the second
# half of each window, a Gaussian copula with correlation 0.2, the exclusion
# rules of Section 2.4 (6 patients, cutoffs 0.95 and 0.90) and the utility
# scores 100, 60, 40 and 0. BF-BOIN-ET-C waits until the DLT outcomes of the
# cohort are known (suspend = 0) and BF-BOIN-ET-A allows a third of them to be
# pending (suspend = 0.33). Time is counted in months of 30 days, as in the
# example of bfboinet::get.oc.backboinet() (windows of 1 or 2 months, three
# patients a month), which gives the durations of the article in months and
# its backfilling period of 100 months after the last cohort. The article does
# not state the distribution of the times between arrivals; the uniform
# distribution of bfboinet is used.
#
# The OD of a scenario is the dose in boldface in Table S1: the dose with a DLT
# probability of at most 0.30 and a response probability of at least 0.25
# that maximizes the utility; Scenarios 13 and 14 have none, and PCS is then
# the percentage of trials selecting no OD. A toxic dose is a dose with a DLT
# probability above 0.30. PCA and PTA are ratios of the total numbers of
# patients over all trials, as in bfboinet.
#
# The article rounds to one decimal, so a cell is judged PASS when it lies
# within four standard errors (of the difference between two independent
# simulations) plus 0.05 for the rounding, and FAIL otherwise.
#
# The scenarios and the published values were extracted from the text of the
# article (Tables 1 and 2) and of its Supporting Information (Table S1); the
# boldface of Table S1 was read from the rendered page.
#
# Usage, from the package root:
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Takeda2025.R")
#   res <- run_reproduce_takeda2025()
#
# The judgements are written to inst/reproduce/results/reproduce_Takeda2025.csv.

takeda2025_scenarios <- function() {
  tox <- list(
    c(0.03, 0.05, 0.20, 0.22, 0.45),
    c(0.02, 0.05, 0.07, 0.10, 0.15),
    c(0.01, 0.05, 0.10, 0.12, 0.27),
    c(0.05, 0.07, 0.10, 0.20, 0.35),
    c(0.03, 0.06, 0.10, 0.30, 0.45),
    c(0.02, 0.05, 0.10, 0.20, 0.30),
    c(0.01, 0.08, 0.10, 0.20, 0.25),
    c(0.15, 0.20, 0.25, 0.35, 0.45),
    c(0.05, 0.12, 0.27, 0.35, 0.50),
    c(0.01, 0.15, 0.20, 0.40, 0.50),
    c(0.15, 0.30, 0.35, 0.40, 0.45),
    c(0.25, 0.35, 0.55, 0.60, 0.65),
    c(0.50, 0.55, 0.60, 0.65, 0.70),
    c(0.25, 0.35, 0.55, 0.60, 0.65)
  )
  eff <- list(
    c(0.05, 0.10, 0.50, 0.68, 0.70),
    c(0.05, 0.08, 0.15, 0.30, 0.45),
    c(0.10, 0.15, 0.25, 0.50, 0.50),
    c(0.05, 0.15, 0.30, 0.55, 0.55),
    c(0.10, 0.20, 0.40, 0.45, 0.50),
    c(0.05, 0.15, 0.40, 0.40, 0.40),
    c(0.07, 0.10, 0.35, 0.20, 0.20),
    c(0.25, 0.55, 0.40, 0.30, 0.20),
    c(0.20, 0.35, 0.36, 0.37, 0.38),
    c(0.45, 0.45, 0.45, 0.65, 0.80),
    c(0.30, 0.30, 0.30, 0.30, 0.30),
    c(0.40, 0.41, 0.42, 0.43, 0.44),
    c(0.05, 0.10, 0.20, 0.30, 0.40),
    c(0.01, 0.02, 0.03, 0.04, 0.05)
  )
  od <- c(4, 5, 4, 4, 3, 3, 3, 2, 2, 1, 1, 1, 0, 0)
  lapply(seq_along(tox), function(s) {
    list(name = paste("Scenario", s), p_true = tox[[s]], p_resp = eff[[s]],
         od = od[s])
  })
}

# Rows of Tables 1 and 2 (PCS, PCA, PTS, PTA, N, Dur), one per scenario.
takeda2025_tables <- function() {
  list(
    "Table 1 (30 days)" = list(
      "BF-BOIN-ET-C" = rbind(
        c(48.9, 26.2, 6.5, 14.5, 36.3, 16.9),
        c(65.6, 33.6, 0.0, 0.0, 34.3, 16.0),
        c(57.7, 26.6, 0.0, 0.0, 39.2, 17.0),
        c(59.9, 28.8, 15.6, 19.6, 38.3, 17.3),
        c(53.8, 29.6, 5.0, 11.6, 38.7, 17.0),
        c(53.3, 26.6, 0.0, 0.0, 38.9, 17.1),
        c(64.4, 26.7, 0.0, 0.0, 36.9, 16.9),
        c(68.8, 30.7, 1.2, 19.3, 34.8, 15.2),
        c(50.7, 30.0, 7.0, 22.1, 37.3, 16.2),
        c(46.4, 24.9, 15.6, 21.8, 41.5, 16.4),
        c(63.6, 37.3, 7.6, 27.2, 29.2, 13.8),
        c(74.7, 53.3, 18.3, 46.7, 22.9, 11.2),
        c(78.8, 0.0, 21.2, 100.0, 9.7, 5.7),
        c(56.3, 0.0, 22.3, 53.4, 17.9, 10.8)
      ),
      "BF-BOIN-ET-A" = rbind(
        c(46.9, 25.1, 5.4, 15.2, 34.4, 13.7),
        c(51.4, 35.7, 0.0, 0.0, 30.8, 12.8),
        c(57.0, 22.8, 0.0, 0.0, 34.9, 13.8),
        c(59.4, 26.4, 13.1, 20.2, 35.0, 13.9),
        c(54.6, 28.0, 3.9, 12.0, 35.7, 13.9),
        c(51.4, 25.3, 0.0, 0.0, 35.0, 13.9),
        c(62.2, 25.0, 0.0, 0.0, 33.6, 13.5),
        c(66.2, 29.5, 1.2, 21.2, 32.7, 12.5),
        c(49.9, 28.7, 6.7, 23.1, 34.2, 13.1),
        c(46.4, 23.1, 14.9, 24.1, 37.0, 13.7),
        c(68.0, 36.7, 8.9, 28.3, 27.1, 10.6),
        c(78.7, 51.0, 17.6, 49.0, 21.7, 8.6),
        c(76.4, 0.0, 23.6, 100.0, 10.4, 4.6),
        c(58.4, 0.0, 22.5, 51.6, 17.6, 7.7)
      )
    ),
    "Table 2 (60 days)" = list(
      "BF-BOIN-ET-C" = rbind(
        c(47.8, 26.1, 6.8, 14.5, 36.6, 17.7),
        c(63.7, 33.3, 0.0, 0.0, 34.5, 17.0),
        c(56.3, 26.5, 0.0, 0.0, 39.5, 17.9),
        c(58.1, 28.7, 15.1, 19.6, 38.5, 18.2),
        c(53.7, 29.5, 4.8, 11.8, 39.1, 18.0),
        c(54.2, 26.7, 0.0, 0.0, 39.1, 18.0),
        c(64.5, 27.0, 0.0, 0.0, 37.0, 17.8),
        c(67.7, 30.8, 0.7, 18.8, 34.4, 16.0),
        c(50.3, 29.7, 7.1, 22.3, 37.5, 17.2),
        c(48.6, 24.7, 13.5, 21.9, 41.7, 17.3),
        c(64.1, 37.2, 7.0, 27.2, 29.1, 14.7),
        c(76.1, 53.4, 16.8, 46.6, 22.8, 12.2),
        c(79.8, 0.0, 20.2, 100.0, 9.7, 6.7),
        c(57.0, 0.0, 22.3, 53.4, 17.9, 11.8)
      ),
      "BF-BOIN-ET-A" = rbind(
        c(45.5, 24.9, 5.3, 15.2, 34.3, 14.6),
        c(43.8, 35.3, 0.0, 0.0, 30.8, 13.7),
        c(59.3, 22.6, 0.0, 0.0, 34.9, 14.6),
        c(59.1, 26.3, 12.0, 20.1, 34.9, 14.8),
        c(52.7, 27.6, 4.7, 12.4, 36.3, 15.1),
        c(52.1, 25.3, 0.0, 0.0, 35.0, 14.8),
        c(63.6, 25.1, 0.0, 0.0, 33.7, 14.5),
        c(65.8, 29.5, 1.3, 20.8, 32.7, 13.4),
        c(50.6, 28.6, 6.1, 23.1, 34.2, 14.1),
        c(46.9, 22.7, 14.7, 24.7, 37.2, 14.6),
        c(68.3, 36.3, 7.7, 29.0, 27.2, 11.6),
        c(78.9, 50.7, 16.8, 49.2, 21.7, 9.4),
        c(78.7, 0.0, 21.3, 100.0, 10.3, 5.6),
        c(59.0, 0.0, 22.9, 51.5, 17.6, 8.7)
      )
    )
  )
}

# PCS, PCA, PTS, PTA, N and Dur of one scenario from the trials kept by
# sim_bf_boin_et(), with the standard errors of the difference between the
# article (1000 trials) and this run.
takeda2025_summary <- function(oc, scenario, target = 0.30) {
  tr <- oc$trials
  n_trials <- nrow(tr$n_pts)
  factor_se <- sqrt(1 / 1000 + 1 / n_trials)
  total <- rowSums(tr$n_pts)
  toxic <- which(scenario$p_true > target)
  # Percentage of trials meeting a condition.
  pct <- function(hit) {
    p <- mean(hit)
    c(100 * p, 100 * sqrt(p * (1 - p)) * factor_se)
  }
  # Ratio of the total numbers of patients, with the standard error of the
  # ratio estimator.
  ratio <- function(part) {
    r <- sum(part) / sum(total)
    c(100 * r, 100 * stats::sd(part - r * total) / mean(total) * factor_se)
  }
  if (scenario$od > 0) {
    pcs <- pct(!is.na(tr$obd) & tr$obd == scenario$od)
    pca <- ratio(tr$n_pts[, scenario$od])
  } else {
    pcs <- pct(is.na(tr$obd))
    pca <- c(0, 0)
  }
  pts <- pct(!is.na(tr$obd) & tr$obd %in% toxic)
  pta <- if (length(toxic) > 0) {
    ratio(rowSums(tr$n_pts[, toxic, drop = FALSE]))
  } else {
    c(0, 0)
  }
  n <- c(mean(total), stats::sd(total) * factor_se)
  dur <- c(mean(tr$duration), stats::sd(tr$duration) * factor_se)
  out <- rbind(pcs, pca, pts, pta, n, dur)
  colnames(out) <- c("value", "se")
  out
}

reproduce_takeda2025_tables <- function(n_trials = 10000, seed = 2025) {
  scenarios <- takeda2025_scenarios()
  pub <- takeda2025_tables()
  windows <- c("Table 1 (30 days)" = 1, "Table 2 (60 days)" = 2)
  suspend <- c("BF-BOIN-ET-C" = 0, "BF-BOIN-ET-A" = 0.33)
  rows <- list()
  for (tb in names(windows)) {
    for (design in names(suspend)) {
      for (s in seq_along(scenarios)) {
        oc <- sim_bf_boin_et(
          target = 0.30, target_resp = 0.25, p_true = scenarios[[s]]$p_true,
          p_resp = scenarios[[s]]$p_resp, n_cohort = 10, cohort_size = 3,
          window = 1, resp_window = windows[[tb]], accrual_rate = 3,
          resp_cor = 0.2, n_trials = n_trials, n_earlystop = 12,
          n_elimination = 6, suspend = suspend[[design]], cutoff_eli = 0.95,
          cutoff_eli_resp = 0.90, keep_trials = TRUE, seed = seed + s
        )
        sm <- takeda2025_summary(oc, scenarios[[s]])
        rows[[length(rows) + 1L]] <- data.frame(
          item = scenarios[[s]]$name,
          table = tb,
          design = design,
          quantity = c("PCS", "PCA", "PTS", "PTA", "N", "Dur"),
          published = pub[[tb]][[design]][s, ],
          reproduced = sm[, "value"],
          tolerance = 4 * sm[, "se"] + 0.05,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

run_reproduce_takeda2025 <- function(n_trials = 10000, seed = 2025,
                                     out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_bf_boin_et", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- reproduce_takeda2025_tables(n_trials = n_trials, seed = seed)
  res$judgement <- ifelse(abs(res$published - res$reproduced) <=
                            res$tolerance + 1e-12, "PASS", "FAIL")
  res$reproduced <- round(res$reproduced, 2)
  res$tolerance <- round(res$tolerance, 2)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Takeda2025.csv"),
                   row.names = FALSE)
  res
}
