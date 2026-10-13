# Comparison of the drug-combination design with Lin and Yin (2017)
# -----------------------------------------------------------------------------
#
# Lin, R. and Yin, G. (2017). Bayesian optimal interval design for dose finding
# in drug-combination trials. Statistical Methods in Medical Research, 26(5),
# 2155-2167.
#
# The article reports two tables that can be set beside the package:
#
#   Table 3: for scenarios 1 to 4 and the safety cutoffs lambda = 0.95 and 0.9,
#            the percentage of trials selecting an MTD, the percentage of
#            patients allocated to the MTDs, the percentage of patients with a
#            DLT and the percentage of trials terminated early because the
#            lowest combination was overly toxic.
#   Table 2: the accuracy index of the ten scenarios, for lambda = 0.95.
#
# Every row is reported as INFO, not as PASS or FAIL, because the two cannot
# agree exactly:
#
#   * The package follows BOIN::get.oc.comb(), not the article, where the two
#     differ: a candidate's posterior probability is increased by 0.0005 per
#     patient treated there, the MTD is selected from (y + 0.05) / (n + 0.1)
#     with ties broken towards lower combinations, whereas the article uses
#     y / n and its own tie rules.
#   * The article simulated 1,000 trials per scenario, so its percentages carry
#     a Monte Carlo standard error of up to 1.6 points.
#   * Table 2 was simulated without early termination ("we do not impose the
#     early stopping rules"), whereas the package always stops a trial whose
#     lowest combination is eliminated, as Table 3 of the article does. This
#     matters in scenario 4, where the lowest combination is the MTD.
#   * The article states its elimination rule only as
#     Pr(p > 0.3 | data) >= lambda. In the reference implementation, two DLTs
#     among the first three patients give 0.916 under its uniform prior (0.911
#     under the Jeffreys prior of the article) and so eliminate a combination
#     when lambda = 0.9. In scenario 2, where the lowest combination has a DLT
#     rate of 0.15, the first cohort shows two or more DLTs with probability
#     0.061, whereas Table 3 reports early termination in 1.2% of the trials
#     for lambda = 0.9. The article does not say what else its rule requires.
#
# The settings are those of Section 3 of the article: target 0.3, 20 cohorts of
# three (60 patients), no early stopping on the sample size at a combination
# (n_earlystop = 100), the boundaries of p_saf = 0.6 * target and
# p_tox = 1.4 * target, and the first cohort at the lowest combination.
#
# Usage
#
#   library(simFastBOIN)
#   source(system.file("validation", "compare-with-LinYin2017.R",
#                      package = "simFastBOIN"))
#
#   results <- compare_with_lin_yin()             # 10,000 trials per setting
#   summarize_lin_yin(results)
#   write.csv(results, "simFastBOIN-vs-LinYin2017.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Table 1 of the article: rows are the levels 1 to 3 of agent 2, columns the
# levels 1 to 5 of agent 1 (the article prints level 3 of agent 2 at the top).
# -----------------------------------------------------------------------------

lin_yin_scenarios <- function() {
  rows <- function(level1, level2, level3) {
    matrix(c(level1, level2, level3), nrow = 3, byrow = TRUE)
  }
  list(
    "1" = rows(c(0.05, 0.10, 0.15, 0.30, 0.45), c(0.10, 0.15, 0.30, 0.45, 0.55),
               c(0.15, 0.30, 0.45, 0.50, 0.60)),
    "2" = rows(c(0.15, 0.30, 0.45, 0.50, 0.60), c(0.30, 0.45, 0.50, 0.60, 0.75),
               c(0.45, 0.55, 0.60, 0.70, 0.80)),
    "3" = rows(c(0.02, 0.07, 0.10, 0.15, 0.30), c(0.07, 0.10, 0.15, 0.30, 0.45),
               c(0.10, 0.15, 0.30, 0.45, 0.55)),
    "4" = rows(c(0.30, 0.45, 0.60, 0.70, 0.80), c(0.45, 0.55, 0.65, 0.75, 0.85),
               c(0.50, 0.60, 0.70, 0.80, 0.90)),
    "5" = rows(c(0.01, 0.02, 0.08, 0.10, 0.11), c(0.03, 0.05, 0.10, 0.13, 0.15),
               c(0.07, 0.09, 0.12, 0.15, 0.30)),
    "6" = rows(c(0.05, 0.08, 0.10, 0.13, 0.15), c(0.09, 0.12, 0.15, 0.30, 0.45),
               c(0.15, 0.30, 0.45, 0.50, 0.60)),
    "7" = rows(c(0.07, 0.10, 0.12, 0.15, 0.30), c(0.15, 0.30, 0.45, 0.52, 0.60),
               c(0.30, 0.50, 0.60, 0.65, 0.75)),
    "8" = rows(c(0.02, 0.10, 0.15, 0.50, 0.60), c(0.05, 0.12, 0.30, 0.55, 0.70),
               c(0.08, 0.15, 0.45, 0.60, 0.80)),
    "9" = rows(c(0.005, 0.01, 0.02, 0.04, 0.07), c(0.02, 0.05, 0.08, 0.12, 0.15),
               c(0.15, 0.30, 0.45, 0.55, 0.65)),
    "10" = rows(c(0.05, 0.10, 0.15, 0.30, 0.45), c(0.45, 0.50, 0.60, 0.65, 0.70),
                c(0.70, 0.75, 0.80, 0.85, 0.90))
  )
}

# Published values, read from Tables 2 and 3 of the article.
lin_yin_published <- function() {
  table3 <- rbind(
    data.frame(lambda = 0.95, scenario = c("1", "2", "3", "4"),
               mtd_selection = c(68.1, 70.9, 67.5, 62.5),
               patients_at_mtd = c(41.9, 50.7, 39.5, 68.6),
               patients_with_dlt = c(26.5, 30.6, 24.7, 30.2),
               early_termination = c(0.0, 0.9, 0.0, 20.3),
               stringsAsFactors = FALSE),
    data.frame(lambda = 0.90, scenario = c("1", "2", "3", "4"),
               mtd_selection = c(68.3, 67.6, 64.8, 54.3),
               patients_at_mtd = c(42.8, 49.7, 39.4, 65.2),
               patients_with_dlt = c(26.5, 30.3, 24.4, 28.3),
               early_termination = c(0.0, 1.2, 0.0, 28.6),
               stringsAsFactors = FALSE)
  )
  table2 <- data.frame(
    lambda = 0.95, scenario = as.character(1:10),
    accuracy_index = c(66.3, 81.1, 68.1, 93.5, 80.7, 59.9, 77.2, 60.8, 57.8,
                       72.3),
    stringsAsFactors = FALSE
  )
  list(table2 = table2, table3 = table3)
}

# Accuracy index of the article, in percent.
lin_yin_accuracy_index <- function(sel_percent, p_true, target) {
  q <- sel_percent / 100
  distance <- abs(p_true - target)
  (1 - length(p_true) * sum(q * distance) / sum(distance)) * 100
}

# -----------------------------------------------------------------------------
# Simulations
# -----------------------------------------------------------------------------

compare_with_lin_yin <- function(n_trials = 10000, seed = 6, verbose = TRUE) {

  scenarios <- lin_yin_scenarios()
  published <- lin_yin_published()
  target <- 0.30

  run <- function(scenario, lambda) {
    sim_comb_boin(
      target = target, p_true = scenarios[[scenario]], n_cohort = 20,
      cohort_size = 3, n_earlystop = 100, cutoff_eli = lambda,
      n_trials = n_trials, seed = seed
    )
  }

  rows <- list()

  for (i in seq_len(nrow(published$table3))) {
    pub <- published$table3[i, ]
    if (verbose) {
      message("Table 3: scenario ", pub$scenario, ", lambda ", pub$lambda)
    }
    oc <- run(pub$scenario, pub$lambda)
    stops <- oc$stop_reason_percent
    terminated <- if ("lowest_dose_eliminated" %in% names(stops)) {
      stops[["lowest_dose_eliminated"]]
    } else {
      0
    }
    ours <- c(
      mtd_selection = oc$pct_true_mtd_selected,
      patients_at_mtd = oc$pct_patients_true_mtd,
      patients_with_dlt = oc$total_n_tox / oc$total_n_pts * 100,
      early_termination = terminated
    )
    for (quantity in names(ours)) {
      rows[[length(rows) + 1L]] <- data.frame(
        table = "Table 3", scenario = pub$scenario, lambda = pub$lambda,
        quantity = quantity, article = pub[[quantity]],
        simFastBOIN = round(ours[[quantity]], 1),
        stringsAsFactors = FALSE
      )
    }
  }

  for (i in seq_len(nrow(published$table2))) {
    pub <- published$table2[i, ]
    if (verbose) message("Table 2: scenario ", pub$scenario)
    oc <- run(pub$scenario, pub$lambda)
    ai <- lin_yin_accuracy_index(oc$sel_percent, scenarios[[pub$scenario]],
                                 target)
    rows[[length(rows) + 1L]] <- data.frame(
      table = "Table 2", scenario = pub$scenario, lambda = pub$lambda,
      quantity = "accuracy_index", article = pub$accuracy_index,
      simFastBOIN = round(ai, 1), stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  out$difference <- round(out$simFastBOIN - out$article, 1)
  out$status <- "INFO"
  attr(out, "n_trials") <- n_trials
  attr(out, "seed") <- seed
  out
}

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------

summarize_lin_yin <- function(results) {
  cat("simFastBOIN drug-combination design against Lin and Yin (2017)\n")
  cat("  trials per setting : ", attr(results, "n_trials"), "\n", sep = "")
  cat("  seed               : ", attr(results, "seed"), "\n", sep = "")
  cat("  status             : INFO for every row (see the header of the",
      "script for why\n                       the two cannot agree exactly)\n")
  cat("\n")
  print(results, row.names = FALSE)
  cat("\nLargest absolute difference by quantity:\n")
  worst <- tapply(abs(results$difference), results$quantity, max)
  print(worst)
  invisible(results)
}
