# Comparison of simFastBOIN with Zhao and colleagues (2024)
#
# Zhao Y, Yuan Y, Korn EL, Freidlin B. Backfilling patients in phase I
# dose-escalation trials using Bayesian optimal interval design (BOIN).
# Clinical Cancer Research 2024; 30(4): 673-679.
# doi:10.1158/1078-0432.CCR-23-2585
#
# The published operating characteristics come from an independent
# implementation with its own random numbers, so agreement is only expected
# within Monte Carlo error. As a yardstick, the BOIN rows of the main Table 4
# and of Supplementary Table A2 are two independent runs of the same design
# (the accrual rate does not affect BOIN without backfilling) and differ by up
# to 1.3 percentage points in the correct selection.
#
# Usage, from the package root after devtools::load_all():
#
#   source("inst/validation/compare-with-Zhao2024.R")
#   t1 <- compare_zhao2024_table1()
#   oc <- compare_zhao2024_boin()
#   bf <- compare_zhao2024_bf_boin()
#
# compare_zhao2024_boin() settles how n_stop is read (n_earlystop = 9 with the
# rule "with_stay"). compare_zhao2024_bf_boin() then compares sim_bf_boin()
# with the BF-BOIN and BOIN rows of Table 4, durations included, averaged over
# several seeds.
#
# All values were transcribed from the images of the article and its
# supplementary appendix, not from extracted text.

# Main Table 1: BOIN decision table for a target of 0.25 with an elimination
# cutoff of 0.95, where one DLT out of three stays instead of de-escalating.
zhao2024_table1 <- function() {
  data.frame(
    n = 3:30,
    escalate = c(0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3,
                 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5),
    deescalate = c(2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5,
                   6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9),
    eliminate = c(3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7,
                  8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12)
  )
}

# Main Table 3: true DLT and response rates. The MTD is the dose at 0.25.
zhao2024_scenarios <- function() {
  dlt <- list(
    c(0.25, 0.41, 0.45, 0.49, 0.53),
    c(0.12, 0.25, 0.42, 0.49, 0.55),
    c(0.04, 0.12, 0.25, 0.43, 0.63),
    c(0.02, 0.06, 0.10, 0.25, 0.40),
    c(0.02, 0.05, 0.08, 0.11, 0.25)
  )
  response <- list(
    c(0.30, 0.40, 0.45, 0.50, 0.55),
    c(0.20, 0.30, 0.40, 0.50, 0.60),
    c(0.10, 0.20, 0.30, 0.45, 0.58),
    c(0.05, 0.10, 0.15, 0.30, 0.45),
    c(0.05, 0.10, 0.15, 0.20, 0.30),
    c(0.30, 0.35, 0.36, 0.36, 0.36),
    c(0.15, 0.30, 0.35, 0.36, 0.36),
    c(0.10, 0.20, 0.30, 0.35, 0.35),
    c(0.10, 0.15, 0.20, 0.30, 0.35),
    c(0.30, 0.32, 0.35, 0.36, 0.36),
    c(0.10, 0.30, 0.32, 0.35, 0.36),
    c(0.10, 0.15, 0.30, 0.32, 0.35)
  )
  # Scenarios 6 to 12 reuse the DLT rates of scenarios 2 to 5.
  dlt_index <- c(1, 2, 3, 4, 5, 2, 3, 4, 5, 3, 4, 5)
  lapply(seq_along(response), function(s) {
    list(p_true = dlt[[dlt_index[s]]], p_resp = response[[s]])
  })
}

# BOIN rows of the main Table 4 and of Supplementary Tables A2 and A16.
# Table 4 and A2 stay on one DLT out of three, A16 does not. A2 uses an
# accrual rate of 6 per month, the others 3 per month. Patient counts are
# averages per trial; NA marks cells the article reports as not applicable.
zhao2024_boin_published <- function() {
  data.frame(
    table = rep(c("Table 4", "Table A2", "Table A16"), each = 5),
    scenario = rep(1:5, times = 3),
    stay_on_1_of_3 = rep(c(TRUE, TRUE, FALSE), each = 5),
    correct = c(77.0, 55.6, 56.7, 53.3, 63.9,
                76.7, 56.2, 56.2, 52.0, 63.7,
                78.6, 55.2, 54.5, 53.2, 59.5),
    over_selection = c(16.7, 14.5, 12.3, 17.1, NA,
                       16.5, 14.3, 12.5, 17.0, NA,
                       15.0, 12.6, 11.2, 14.9, NA),
    n_below = c(NA, 7.2, 11.4, 14.0, 17.6,
                NA, 7.3, 11.4, 14.1, 17.7,
                NA, 8.1, 13.2, 15.8, 19.7),
    n_at = c(8.7, 8.4, 8.0, 7.7, 7.3,
             8.8, 8.4, 7.9, 7.6, 7.2,
             9.1, 8.4, 7.6, 7.5, 6.4),
    n_over = c(6.2, 5.1, 4.6, 4.0, NA,
               6.2, 5.0, 4.5, 3.9, NA,
               5.6, 4.4, 3.8, 3.2, NA),
    n_total = c(15.1, 20.7, 23.9, 25.7, 24.9,
                15.0, 20.7, 23.9, 25.6, 24.9,
                14.6, 20.8, 24.6, 26.5, 26.1),
    duration = c(8.8, 12.3, 14.3, 15.5, 15.1,
                 6.8, 9.5, 11.0, 11.9, 11.7,
                 8.5, 12.4, 14.8, 16.1, 15.9)
  )
}

# Compare boin_boundary() with the main Table 1.
compare_zhao2024_table1 <- function() {
  published <- zhao2024_table1()
  bound <- boin_boundary(target = 0.25, max_n = 30, cutoff_eli = 0.95,
                         stay_on_1_of_3 = TRUE)
  idx <- published$n
  out <- data.frame(
    n = idx,
    escalate_published = published$escalate,
    escalate_package = bound$b_esc[idx],
    deescalate_published = published$deescalate,
    deescalate_package = bound$b_deesc[idx],
    eliminate_published = published$eliminate,
    eliminate_package = bound$b_elim[idx]
  )
  out$match <- out$escalate_published == out$escalate_package &
    out$deescalate_published == out$deescalate_package &
    out$eliminate_published == out$eliminate_package
  out
}

# Operating characteristics of sim_boin() in the layout of Table 4.
zhao2024_boin_metrics <- function(oc) {
  mtd <- which(abs(oc$p_true - 0.25) < 1e-8)
  n_doses <- length(oc$p_true)
  below <- seq_len(mtd - 1)
  above <- if (mtd < n_doses) (mtd + 1):n_doses else integer(0)
  c(
    correct = unname(oc$sel_percent[mtd]),
    over_selection = if (length(above) > 0) sum(oc$sel_percent[above]) else NA,
    n_below = if (length(below) > 0) sum(oc$n_pts_dose[below]) else NA,
    n_at = unname(oc$n_pts_dose[mtd]),
    n_over = if (length(above) > 0) sum(oc$n_pts_dose[above]) else NA,
    n_total = oc$total_n_pts
  )
}

# Run sim_boin() for every published BOIN row under each candidate reading of
# the n_stop rule and return published and simulated values side by side.
# n_earlystop = 9 stops once 9 patients have been treated at the current dose
# (the main text and the trial example of Supplementary Section B);
# n_earlystop = 12 stops once more than 9 have been treated, since cohorts of 3
# move from 9 to 12 (Supplementary Section C says "exceeded").
compare_zhao2024_boin <- function(n_earlystop = c(9, 12),
                                  n_earlystop_rule = c("with_stay", "simple"),
                                  n_trials = 10000, seed = 123) {
  published <- zhao2024_boin_published()
  scenarios <- zhao2024_scenarios()
  metrics <- c("correct", "over_selection", "n_below", "n_at", "n_over",
               "n_total")
  grid <- expand.grid(n_earlystop = n_earlystop,
                      n_earlystop_rule = n_earlystop_rule,
                      stay_on_1_of_3 = c(TRUE, FALSE),
                      scenario = 1:5,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    oc <- sim_boin(target = 0.25, p_true = scenarios[[g$scenario]]$p_true,
                   n_cohort = 10, cohort_size = 3, n_trials = n_trials,
                   n_earlystop = g$n_earlystop, cutoff_eli = 0.95,
                   stay_on_1_of_3 = g$stay_on_1_of_3,
                   n_earlystop_rule = g$n_earlystop_rule, seed = seed)
    simulated <- zhao2024_boin_metrics(oc)
    pub <- published[published$scenario == g$scenario &
                       published$stay_on_1_of_3 == g$stay_on_1_of_3, ]
    rows[[i]] <- do.call(rbind, lapply(seq_len(nrow(pub)), function(j) {
      data.frame(
        table = pub$table[j],
        scenario = g$scenario,
        stay_on_1_of_3 = g$stay_on_1_of_3,
        n_earlystop = g$n_earlystop,
        n_earlystop_rule = g$n_earlystop_rule,
        metric = metrics,
        published = unlist(pub[j, metrics]),
        simulated = unname(simulated[metrics]),
        row.names = NULL
      )
    }))
  }
  out <- do.call(rbind, rows)
  out$difference <- out$simulated - out$published
  out
}

# BF-BOIN rows of the main Table 4, scenarios 1 to 12, with an accrual rate of
# 3 per month. NA marks cells the article reports as not applicable.
zhao2024_bf_boin_published <- function() {
  data.frame(
    table = "Table 4",
    scenario = 1:12,
    correct = c(79.9, 57.8, 57.6, 56.7, 62.8, 57.5, 57.1, 57.0, 62.9, 57.2,
                57.0, 62.8),
    over_selection = c(13.5, 10.9, 9.7, 13.6, NA, 10.8, 9.6, 13.2, NA, 9.5,
                       13.1, NA),
    n_below = c(NA, 10.4, 16.9, 20.4, 27.2, 11.1, 18.3, 23.5, 29.5, 20.1,
                24.3, 30.0),
    n_at = c(10.8, 10.3, 9.8, 9.6, 7.2, 10.3, 9.9, 9.7, 7.2, 9.8, 9.7, 7.1),
    n_over = c(6.1, 4.8, 4.2, 3.7, NA, 4.8, 4.2, 3.6, NA, 4.2, 3.5, NA),
    n_total = c(17.0, 25.5, 30.9, 33.8, 34.4, 26.2, 32.3, 36.7, 36.7, 34.2,
                37.5, 37.1),
    duration = c(8.4, 12.1, 14.3, 15.6, 15.5, 12.1, 14.3, 15.6, 15.5, 14.3,
                 15.6, 15.5)
  )
}

# Run sim_bf_boin() with the settings of the article for every BF-BOIN row of
# Table 4, and for the BOIN rows with no response so that nothing is
# backfilled, and return published and simulated values side by side. The
# simulated values are averages over the seeds. The article states a Poisson
# process, but its durations agree with uniform times between arrivals and
# patients turned away while the escalation waits, which is what is used here.
compare_zhao2024_bf_boin <- function(seeds = 1:5, n_trials = 10000) {
  scenarios <- zhao2024_scenarios()
  metrics <- c("correct", "over_selection", "n_below", "n_at", "n_over",
               "n_total", "duration")
  bf <- zhao2024_bf_boin_published()
  boin <- zhao2024_boin_published()
  boin <- boin[boin$table == "Table 4", ]
  rows <- rbind(
    data.frame(method = "BF-BOIN", bf[, c("scenario", metrics)]),
    data.frame(method = "BOIN", boin[, c("scenario", metrics)])
  )

  out <- vector("list", nrow(rows))
  for (i in seq_len(nrow(rows))) {
    s <- rows$scenario[i]
    p_resp <- if (rows$method[i] == "BOIN") rep(0, 5) else scenarios[[s]]$p_resp
    runs <- sapply(seeds, function(seed) {
      oc <- sim_bf_boin(
        target = 0.25, p_true = scenarios[[s]]$p_true, p_resp = p_resp,
        n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
        n_earlystop = 9, stay_on_1_of_3 = TRUE, accrual = "uniform",
        no_slot = "leave", n_trials = n_trials, seed = seed
      )
      c(zhao2024_boin_metrics(oc), duration = oc$duration_mean)
    })
    out[[i]] <- data.frame(
      method = rows$method[i],
      scenario = s,
      metric = metrics,
      published = unlist(rows[i, metrics], use.names = FALSE),
      simulated = unname(rowMeans(runs[metrics, , drop = FALSE])),
      row.names = NULL
    )
  }
  out <- do.call(rbind, out)
  out$difference <- out$simulated - out$published
  out
}
