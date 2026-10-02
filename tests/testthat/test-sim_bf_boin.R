zhao2024_scenarios <- function() {
  dlt <- list(
    c(0.25, 0.41, 0.45, 0.49, 0.53),
    c(0.12, 0.25, 0.42, 0.49, 0.55),
    c(0.04, 0.12, 0.25, 0.43, 0.63),
    c(0.02, 0.06, 0.10, 0.25, 0.40),
    c(0.02, 0.05, 0.08, 0.11, 0.25)
  )
  response <- list(
    c(0.30, 0.40, 0.45, 0.50, 0.55), c(0.20, 0.30, 0.40, 0.50, 0.60),
    c(0.10, 0.20, 0.30, 0.45, 0.58), c(0.05, 0.10, 0.15, 0.30, 0.45),
    c(0.05, 0.10, 0.15, 0.20, 0.30), c(0.30, 0.35, 0.36, 0.36, 0.36),
    c(0.15, 0.30, 0.35, 0.36, 0.36), c(0.10, 0.20, 0.30, 0.35, 0.35),
    c(0.10, 0.15, 0.20, 0.30, 0.35), c(0.30, 0.32, 0.35, 0.36, 0.36),
    c(0.10, 0.30, 0.32, 0.35, 0.36), c(0.10, 0.15, 0.30, 0.32, 0.35)
  )
  dlt_index <- c(1, 2, 3, 4, 5, 2, 3, 4, 5, 3, 4, 5)
  lapply(seq_along(response), function(s) {
    list(p_true = dlt[[dlt_index[s]]], p_resp = response[[s]])
  })
}

# Compare sim_bf_boin() with the published rows of one method and return the
# number of published cells checked. The published values come from a
# separate implementation with its own random numbers, so the tolerances allow
# for simulation error: over 20 seeds of 10,000 trials each, the largest gaps
# were 2.1 points for the selection, 0.4 patients and 0.17 months.
check_zhao2024_rows <- function(method) {
  published <- read.csv(test_path("fixtures", "zhao2024-table4.csv"),
                        stringsAsFactors = FALSE)
  published <- published[published$method == method, ]
  scenarios <- zhao2024_scenarios()
  metrics <- c("correct", "over_selection", "n_below", "n_at", "n_over",
               "n_total", "duration")
  tolerance <- c(correct = 3.5, over_selection = 3.5, n_below = 0.8,
                 n_at = 0.8, n_over = 0.8, n_total = 0.8, duration = 0.3)

  checked <- 0L
  for (i in seq_len(nrow(published))) {
    row <- published[i, ]
    s <- row$scenario
    p_resp <- if (method == "BOIN") rep(0, 5) else scenarios[[s]]$p_resp
    oc <- sim_bf_boin(
      target = 0.25, p_true = scenarios[[s]]$p_true, p_resp = p_resp,
      n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
      n_earlystop = 9, stay_on_1_of_3 = TRUE, accrual = "uniform",
      no_slot = "leave", n_trials = 10000, seed = 100 + s
    )
    mtd <- which(abs(scenarios[[s]]$p_true - 0.25) < 1e-8)
    below <- seq_len(mtd - 1L)
    above <- if (mtd < 5L) (mtd + 1L):5L else integer(0)
    simulated <- c(
      correct = unname(oc$sel_percent[mtd]),
      over_selection = if (length(above) > 0L) sum(oc$sel_percent[above]) else NA,
      n_below = if (length(below) > 0L) sum(oc$n_pts_dose[below]) else NA,
      n_at = unname(oc$n_pts_dose[mtd]),
      n_over = if (length(above) > 0L) sum(oc$n_pts_dose[above]) else NA,
      n_total = oc$total_n_pts,
      duration = oc$duration_mean
    )
    expected <- unlist(row[metrics])
    # The not applicable cells agree as well.
    expect_identical(is.na(simulated[metrics]), is.na(expected),
                     info = paste(method, "scenario", s))
    ok <- !is.na(expected)
    gap <- abs(simulated[metrics][ok] - expected[ok])
    off <- names(gap)[gap > tolerance[names(gap)]]
    expect_identical(off, character(0), info = paste(method, "scenario", s))
    checked <- checked + sum(ok)
  }
  checked
}

test_that("the BF-BOIN rows of Table 4 of Zhao et al. (2024) are reproduced", {
  # Twelve rows, with the cells that are not applicable left out.
  expect_equal(check_zhao2024_rows("BF-BOIN"), 12L * 7L - 7L)
})

test_that("the BOIN rows of Table 4 of Zhao et al. (2024) are reproduced", {
  # Without responses BF-BOIN is BOIN, and the durations of the article
  # require uniform times between arrivals and patients turned away.
  expect_equal(check_zhao2024_rows("BOIN"), 5L * 7L - 3L)
})

test_that("without any response the results are those of sim_boin", {
  p_true <- c(0.04, 0.12, 0.25, 0.43, 0.63)
  common <- list(target = 0.25, p_true = p_true, n_cohort = 10,
                 cohort_size = 3, n_earlystop = 9, stay_on_1_of_3 = TRUE,
                 extrasafe = TRUE, bound_mtd = TRUE, min_mtd_sample = 3,
                 overdose_cutoff = 0.3, n_trials = 1000, seed = 8)
  bf <- do.call(sim_bf_boin, c(common, list(p_resp = rep(0, 5), window = 1,
                                            accrual_rate = 3)))
  boin <- do.call(sim_boin, common)

  for (component in c("sel_percent", "percent_no_mtd", "n_pts_dose",
                      "n_tox_dose", "total_n_pts", "total_n_tox", "overdose",
                      "stop_reason_percent")) {
    expect_equal(bf[[component]], boin[[component]], info = component)
  }
  expect_equal(bf$total_n_bf, 0)
  expect_equal(bf$pct_trials_backfill, 0)
})

test_that("the summaries agree with the trial by trial data", {
  oc <- sim_bf_boin(
    target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
    p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58), n_cohort = 10, cohort_size = 3,
    window = 1, accrual_rate = 3, n_trials = 500, keep_trials = TRUE,
    seed = 9
  )
  trials <- oc$trials

  expect_s3_class(oc, "bf_boin_oc")
  expect_s3_class(oc, "backfill_oc")
  expect_s3_class(oc, "boin_oc")
  expect_equal(unname(oc$n_bf_dose), unname(colMeans(trials$n_bf)))
  expect_equal(oc$total_n_bf, sum(oc$n_bf_dose))
  expect_equal(oc$total_n_pts, mean(rowSums(trials$n_pts)))
  expect_equal(unname(oc$n_resp_dose), unname(colMeans(trials$n_resp)))
  expect_equal(oc$pct_trials_backfill, 100 * mean(rowSums(trials$n_bf) > 0))
  expect_equal(oc$duration_mean, mean(trials$duration))
  expect_equal(sum(oc$sel_percent) + oc$percent_no_mtd, 100)
  expect_equal(oc$p_resp, c(0.10, 0.20, 0.30, 0.45, 0.58))
  expect_true(is.call(oc$call))
  expect_null(sim_bf_boin(
    target = 0.25, p_true = c(0.10, 0.25, 0.40), p_resp = c(0.2, 0.3, 0.4),
    n_cohort = 5, cohort_size = 3, window = 1, accrual_rate = 3,
    n_trials = 20, seed = 1
  )$trials)
})

test_that("verbose reports progress", {
  expect_message(
    sim_bf_boin(target = 0.25, p_true = c(0.10, 0.25, 0.40),
                p_resp = c(0.2, 0.3, 0.4), n_cohort = 5, cohort_size = 3,
                window = 1, accrual_rate = 3, n_trials = 20, verbose = TRUE,
                seed = 1),
    "Simulating"
  )
})
