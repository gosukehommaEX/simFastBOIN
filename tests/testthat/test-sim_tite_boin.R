test_that("sim_tite_boin returns a coherent summary", {
  oc <- sim_tite_boin(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 500, seed = 123
  )

  expect_s3_class(oc, "tite_boin_oc")
  expect_s3_class(oc, "boin_oc")
  expect_length(oc$sel_percent, 5L)
  expect_equal(sum(oc$sel_percent) + oc$percent_no_mtd, 100, tolerance = 1e-10)
  expect_equal(sum(oc$n_pts_dose), oc$total_n_pts, tolerance = 1e-10)
  expect_equal(sum(oc$n_tox_dose), oc$total_n_tox, tolerance = 1e-10)
  expect_true(oc$total_n_tox <= oc$total_n_pts)
  expect_true(oc$total_n_pts <= 30)
  expect_gt(oc$duration_mean, 0)
  expect_true(oc$pct_trials_suspended >= 0 && oc$pct_trials_suspended <= 100)
  expect_true(oc$avg_time_suspended < oc$duration_mean)
  expect_null(oc$trials)
})

test_that("the summary agrees with the trial level data it is built from", {
  oc <- sim_tite_boin(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 400, keep_trials = TRUE, seed = 4
  )
  tr <- oc$trials

  expect_equal(unname(oc$n_pts_dose), unname(colMeans(tr$n_pts)))
  expect_equal(unname(oc$n_tox_dose), unname(colMeans(tr$n_tox)))
  expect_equal(oc$percent_no_mtd, mean(is.na(tr$mtd)) * 100)
  for (dose in seq_len(5)) {
    expect_equal(unname(oc$sel_percent[dose]),
                 sum(tr$mtd == dose, na.rm = TRUE) / length(tr$mtd) * 100,
                 tolerance = 1e-10)
  }
  expect_equal(oc$duration_mean, mean(tr$duration))
  expect_equal(oc$duration_sd, stats::sd(tr$duration))
  expect_equal(oc$pct_trials_suspended, mean(tr$n_suspensions > 0) * 100)
  expect_equal(oc$avg_n_suspensions, mean(tr$n_suspensions))
  expect_equal(oc$avg_time_suspended, mean(tr$time_suspended))

  # The trials are those of tite_boin_simulate() with the same arguments.
  raw <- tite_boin_simulate(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 400, seed = 4
  )
  expect_identical(tr$n_pts, raw$n_pts)
  expect_identical(tr$duration, raw$duration)
})

test_that("without pending patients the results are those of sim_boin", {
  # Arrivals one unit apart and a window of half a unit, under several options
  # that act on the conduct of the trial and on the selection of the MTD.
  settings <- list(
    list(target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
         n_cohort = 10, cohort_size = 3),
    list(target = 0.30, p_true = c(0.30, 0.40, 0.50, 0.60), n_cohort = 12,
         cohort_size = 3, extrasafe = TRUE, bound_mtd = TRUE,
         overdose_cutoff = 1 / 3),
    list(target = 0.25, p_true = c(0.05, 0.10, 0.20, 0.30, 0.45),
         n_cohort = 15, cohort_size = 2, mtd_max_estimate = 0.25,
         min_mtd_sample = 3, n_earlystop = 12, n_earlystop_rule = "simple")
  )

  for (i in seq_along(settings)) {
    set <- settings[[i]]
    boin <- do.call(sim_boin, c(set, list(n_trials = 500, seed = 200 + i)))
    for (method in c("imputation", "ess")) {
      tite <- do.call(sim_tite_boin, c(set, list(
        window = 0.5, accrual_rate = 1, accrual = "fixed", method = method,
        n_trials = 500, seed = 200 + i
      )))
      info <- paste("setting", i, method)
      expect_identical(tite$sel_percent, boin$sel_percent, info = info)
      expect_identical(tite$percent_no_mtd, boin$percent_no_mtd, info = info)
      expect_identical(tite$n_pts_dose, boin$n_pts_dose, info = info)
      expect_identical(tite$n_tox_dose, boin$n_tox_dose, info = info)
      expect_identical(tite$total_n_pts, boin$total_n_pts, info = info)
      expect_identical(tite$total_n_tox, boin$total_n_tox, info = info)
      expect_identical(tite$overdose, boin$overdose, info = info)
      expect_identical(tite$stop_reason_percent, boin$stop_reason_percent,
                       info = info)
      expect_identical(tite$pct_trials_suspended, 0, info = info)
    }
  }
})

test_that("the design finds the MTD it is meant to find", {
  # The same run in a separate Python implementation selects dose 3 in about
  # 82 percent of trials and dose 2 in about 15 percent.
  oc <- sim_tite_boin(
    target = 0.30, p_true = c(0.02, 0.05, 0.30, 0.55, 0.70),
    n_cohort = 20, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 2000, seed = 21
  )
  expect_equal(which.max(oc$sel_percent), 3L, ignore_attr = TRUE)
})

test_that("faster accrual suspends more often and shortens the trial", {
  # In a separate Python implementation with the same seed, half a patient
  # per unit of time gives 1.4 suspensions per trial and a mean duration of
  # 61, and four patients give 5.0 suspensions and a duration of 19.
  args <- list(target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
               n_cohort = 10, cohort_size = 3, window = 3, n_trials = 500,
               seed = 31)
  slow <- do.call(sim_tite_boin, c(args, list(accrual_rate = 0.5)))
  fast <- do.call(sim_tite_boin, c(args, list(accrual_rate = 4)))

  expect_lt(slow$avg_n_suspensions, fast$avg_n_suspensions)
  expect_gt(slow$duration_mean, fast$duration_mean)
})

test_that("the settings and the call are kept", {
  oc <- sim_tite_boin(
    target = 0.30, p_true = c(0.10, 0.25, 0.40), n_cohort = 5,
    cohort_size = 3, window = 2, accrual_rate = 1.5, method = "ess",
    accrual = "uniform", dlt_time = "uniform", n_trials = 50, seed = 1
  )

  expect_identical(oc$settings$method, "ess")
  expect_identical(oc$settings$accrual, "uniform")
  expect_identical(oc$settings$dlt_time, "uniform")
  expect_equal(oc$settings$window, 2)
  expect_equal(oc$settings$accrual_rate, 1.5)
  expect_identical(oc$settings$min_completed, 2L)
  expect_true(is.call(oc$call))
})

test_that("invalid arguments are rejected", {
  args <- list(target = 0.30, p_true = c(0.10, 0.25, 0.40), n_cohort = 5,
               cohort_size = 3, window = 3, accrual_rate = 2, n_trials = 10)

  expect_error(do.call(sim_tite_boin, c(args, list(overdose_cutoff = 1.2))),
               "overdose_cutoff")
  expect_error(do.call(sim_tite_boin, c(args, list(min_mtd_sample = 0))),
               "min_mtd_sample")
  expect_error(do.call(sim_tite_boin, c(args, list(method = "crm"))),
               "should be one of")
  expect_error(do.call(sim_tite_boin,
                       c(args[names(args) != "window"], list(window = -1))),
               "window")
})
