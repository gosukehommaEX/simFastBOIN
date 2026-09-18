test_that("sim_boin returns a coherent summary", {
  oc <- sim_boin(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, n_trials = 500, seed = 123
  )

  expect_s3_class(oc, "boin_oc")
  expect_length(oc$sel_percent, 5L)
  expect_equal(sum(oc$sel_percent) + oc$percent_no_mtd, 100, tolerance = 1e-10)
  expect_true(all(oc$sel_percent >= 0))
  expect_equal(sum(oc$n_pts_dose), oc$total_n_pts, tolerance = 1e-10)
  expect_equal(sum(oc$n_tox_dose), oc$total_n_tox, tolerance = 1e-10)
  expect_true(oc$total_n_tox <= oc$total_n_pts)
  expect_true(oc$total_n_pts <= 30)
  expect_null(oc$trials)
})

test_that("the summary agrees with the trial level data it is built from", {
  oc <- sim_boin(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, n_trials = 400, keep_trials = TRUE, seed = 4
  )

  expect_false(is.null(oc$trials))
  expect_equal(unname(oc$n_pts_dose), unname(colMeans(oc$trials$n_pts)))
  expect_equal(unname(oc$n_tox_dose), unname(colMeans(oc$trials$n_tox)))
  expect_equal(oc$percent_no_mtd, mean(is.na(oc$trials$mtd)) * 100)
  for (dose in seq_len(5)) {
    expect_equal(unname(oc$sel_percent[dose]),
                 sum(oc$trials$mtd == dose, na.rm = TRUE) /
                   length(oc$trials$mtd) * 100,
                 tolerance = 1e-10)
  }
})

test_that("the design finds the MTD it is meant to find", {
  oc <- sim_boin(
    target = 0.30, p_true = c(0.02, 0.05, 0.30, 0.55, 0.70),
    n_cohort = 20, cohort_size = 3, n_trials = 2000, seed = 21
  )
  expect_equal(which.max(oc$sel_percent), 3L, ignore_attr = TRUE)
})

test_that("an entirely toxic scenario usually selects no dose", {
  oc <- sim_boin(
    target = 0.30, p_true = c(0.70, 0.80, 0.90, 0.95, 0.97),
    n_cohort = 20, cohort_size = 3, n_trials = 1000, seed = 22
  )
  expect_gt(oc$percent_no_mtd, 80)
  expect_true("lowest_dose_eliminated" %in% names(oc$stop_reason_percent))
})

test_that("extrasafe stops more often than the default rule", {
  # Under this scenario the reference design stops without an MTD in about 56
  # percent of trials, rising to about 74 percent with the safety rule, so the
  # ordering below has a wide margin.
  args <- list(target = 0.30, p_true = c(0.45, 0.55, 0.65, 0.75, 0.85),
               n_cohort = 20, cohort_size = 3, n_trials = 1000, seed = 23)

  plain <- do.call(sim_boin, args)
  safe <- do.call(sim_boin, c(args, list(extrasafe = TRUE)))

  expect_gt(safe$percent_no_mtd, plain$percent_no_mtd)
  expect_lt(safe$total_n_pts, plain$total_n_pts)
  expect_true("lowest_dose_too_toxic" %in% names(safe$stop_reason_percent))
  expect_false("lowest_dose_too_toxic" %in% names(plain$stop_reason_percent))
})

test_that("bound_mtd never selects a higher dose than the unbounded rule", {
  args <- list(target = 0.30, p_true = c(0.10, 0.20, 0.28, 0.36, 0.50),
               n_cohort = 20, cohort_size = 3, n_trials = 1000,
               keep_trials = TRUE, seed = 24)

  plain <- do.call(sim_boin, args)
  bounded <- do.call(sim_boin, c(args, list(bound_mtd = TRUE)))

  both <- !is.na(plain$trials$mtd) & !is.na(bounded$trials$mtd)
  expect_true(all(bounded$trials$mtd[both] <= plain$trials$mtd[both]))
  expect_gte(bounded$percent_no_mtd, plain$percent_no_mtd)
})

test_that("titration reduces the number of patients on safe doses", {
  args <- list(target = 0.30, p_true = c(0.01, 0.02, 0.05, 0.10, 0.30),
               n_cohort = 20, cohort_size = 3, n_trials = 1000, seed = 25)

  plain <- do.call(sim_boin, args)
  titrated <- do.call(sim_boin, c(args, list(titration = TRUE)))

  expect_lt(sum(titrated$n_pts_dose[1:3]), sum(plain$n_pts_dose[1:3]))
})

test_that("the overdosing summary is reported", {
  oc <- sim_boin(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 20, cohort_size = 3, n_trials = 500, seed = 26
  )
  expect_equal(oc$overdose$cutoff, 0.30)
  expect_equal(oc$overdose$doses, c(4L, 5L), ignore_attr = TRUE)
  expect_gte(oc$overdose$pct_trials_over_60, oc$overdose$pct_trials_over_80)

  safe_only <- sim_boin(
    target = 0.30, p_true = c(0.02, 0.05, 0.08, 0.12, 0.20),
    n_cohort = 20, cohort_size = 3, n_trials = 200, seed = 27
  )
  expect_length(safe_only$overdose$doses, 0L)
  expect_equal(safe_only$overdose$pct_patients, 0)
})

test_that("the same seed reproduces the same operating characteristics", {
  args <- list(target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
               n_cohort = 10, cohort_size = 3, n_trials = 300)

  expect_equal(do.call(sim_boin, c(args, list(seed = 5)))$sel_percent,
               do.call(sim_boin, c(args, list(seed = 5)))$sel_percent)
  expect_false(isTRUE(all.equal(
    do.call(sim_boin, c(args, list(seed = 5)))$sel_percent,
    do.call(sim_boin, c(args, list(seed = 6)))$sel_percent
  )))
})

test_that("verbose progress goes to the message stream", {
  args <- list(target = 0.30, p_true = c(0.10, 0.25, 0.40),
               n_cohort = 5, cohort_size = 3, n_trials = 50, seed = 1)

  # capture_messages() takes every message, so none escape to the console.
  reported <- capture_messages(do.call(sim_boin, c(args, list(verbose = TRUE))))
  expect_true(any(grepl("Simulating", reported)))
  expect_true(any(grepl("Selecting the MTD", reported)))

  expect_silent(do.call(sim_boin, args))
})
