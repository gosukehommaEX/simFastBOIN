test_that("waterfall_simulate returns a well formed object", {
  trials <- waterfall_simulate(
    target = 0.30, p_true = zhang2016_p_true(12), cohort_size = 3,
    n_trials = 200, seed = 1
  )

  expect_s3_class(trials, "waterfall_trials")
  expect_equal(dim(trials$n_pts), c(200L, 3L, 5L))
  expect_equal(dim(trials$n_tox), c(200L, 3L, 5L))
  expect_equal(dim(trials$eliminated), c(200L, 3L, 5L))
  expect_equal(dimnames(trials$n_pts)[[2L]], c("A1", "A2", "A3"))
  expect_equal(dimnames(trials$n_pts)[[3L]], paste0("B", 1:5))
  expect_length(trials$n_subtrials, 200L)
  expect_length(trials$n_enrolled, 200L)
  expect_length(trials$stop_reason, 200L)
  expect_type(trials$n_pts, "integer")
  expect_type(trials$eliminated, "logical")
  expect_s3_class(trials$boundary, "boin_boundary")
})

test_that("the default numbers of cohorts are those of the reference", {
  # ceiling(round(4 / cohort_size, 2) * doses of the subtrial): the first
  # subtrial holds J + K - 1 doses and the others K - 1.
  cohorts <- function(p_true, cohort_size) {
    waterfall_simulate(target = 0.30, p_true = p_true,
                       cohort_size = cohort_size, n_trials = 2,
                       seed = 1)$settings$n_cohort
  }
  expect_equal(cohorts(zhang2016_p_true(1), 3), c(6L, 3L))
  expect_equal(cohorts(zhang2016_p_true(12), 3), c(10L, 6L, 6L))
  expect_equal(cohorts(matrix(seq(0.05, 0.80, by = 0.05), 4, 4), 3),
               c(10L, 4L, 4L, 4L))
  expect_equal(cohorts(zhang2016_p_true(1), 1), c(16L, 8L))
  expect_equal(cohorts(zhang2016_p_true(1), 2), c(8L, 4L))
})

test_that("simulated counts respect the design constraints", {
  trials <- waterfall_simulate(
    target = 0.30, p_true = waterfall_steep_p_true(), cohort_size = 3,
    n_trials = 500, seed = 2
  )

  totals <- apply(trials$n_pts, 1L, sum)
  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(totals <= trials$settings$max_total_pts))
  expect_true(all(totals > 0))
  expect_identical(as.integer(totals), trials$n_enrolled)
  expect_true(all(trials$n_subtrials >= 1L & trials$n_subtrials <= 3L))
  expect_true(all(trials$stop_reason %in%
                    c("completed", "first_subtrial_no_mtd",
                      "later_subtrial_no_mtd", "max_sample_size")))
  # Every kind of ending occurs in this matrix.
  expect_true(all(c("completed", "later_subtrial_no_mtd") %in%
                    trials$stop_reason))
})

test_that("the trials agree with a separate implementation in Python", {
  seed <- waterfall_expected_seed()
  check <- function(trials, expected) {
    n_trials <- nrow(expected$n_pts)
    expect_equal(matrix(trials$n_pts, nrow = n_trials), expected$n_pts)
    expect_equal(matrix(trials$n_tox, nrow = n_trials), expected$n_tox)
    expect_equal(matrix(as.integer(trials$eliminated), nrow = n_trials),
                 expected$eliminated)
    expect_equal(trials$n_subtrials, expected$n_subtrials)
    expect_equal(trials$n_enrolled, expected$n_enrolled)
    expect_equal(trials$stop_reason, expected$stop_reason)
  }

  check(
    waterfall_simulate(target = 0.30, p_true = zhang2016_p_true(4),
                       cohort_size = 3, n_trials = 8,
                       boin_compatible = TRUE, seed = seed),
    waterfall_expected_compatible()
  )
  check(
    waterfall_simulate(target = 0.30, p_true = zhang2016_p_true(4),
                       cohort_size = 3, n_trials = 8, seed = seed),
    waterfall_expected_default()
  )
  check(
    waterfall_simulate(target = 0.30, p_true = waterfall_steep_p_true(),
                       cohort_size = 2, n_trials = 6, titration = TRUE,
                       extrasafe = TRUE, bound_mtd = TRUE, seed = 11),
    waterfall_expected_options()
  )
})

test_that("boin_compatible only drops the data of stopped subtrials", {
  args <- list(target = 0.30, p_true = zhang2016_p_true(4), cohort_size = 3,
               n_trials = 400, seed = 5)
  kept <- do.call(waterfall_simulate, args)
  dropped <- do.call(waterfall_simulate, c(args, boin_compatible = TRUE))

  # The same patients are treated, and the trials end for the same reason.
  expect_identical(kept$n_enrolled, dropped$n_enrolled)
  expect_identical(kept$stop_reason, dropped$stop_reason)
  expect_identical(kept$n_subtrials, dropped$n_subtrials)

  completed <- kept$stop_reason == "completed"
  first <- kept$stop_reason == "first_subtrial_no_mtd"
  later <- kept$stop_reason == "later_subtrial_no_mtd"
  expect_true(any(first))
  expect_true(any(later))
  expect_identical(kept$n_pts[completed, , ], dropped$n_pts[completed, , ])
  expect_identical(kept$eliminated[completed, , ],
                   dropped$eliminated[completed, , ])

  kept_total <- apply(kept$n_pts, 1L, sum)
  dropped_total <- apply(dropped$n_pts, 1L, sum)
  expect_identical(as.integer(kept_total), kept$n_enrolled)
  expect_true(all(dropped_total[first] == 0L))
  expect_true(all(dropped_total[later] < kept_total[later]))
  # A subtrial without a candidate has all its doses flagged by default.
  expect_true(all(kept$eliminated[first, 1L, 1L]))
})

test_that("the seed makes the trials reproducible and is restored", {
  set.seed(99)
  before <- .Random.seed
  a <- waterfall_simulate(target = 0.30, p_true = zhang2016_p_true(1),
                          cohort_size = 3, n_trials = 50, seed = 7)
  expect_identical(.Random.seed, before)
  b <- waterfall_simulate(target = 0.30, p_true = zhang2016_p_true(1),
                          cohort_size = 3, n_trials = 50, seed = 7)
  expect_identical(a$n_pts, b$n_pts)
  expect_identical(a$n_tox, b$n_tox)
})

test_that("p_saf and p_tox enter the boundaries", {
  a <- waterfall_simulate(target = 0.30, p_true = zhang2016_p_true(12),
                          cohort_size = 3, n_trials = 300, seed = 3)
  b <- waterfall_simulate(target = 0.30, p_true = zhang2016_p_true(12),
                          cohort_size = 3, n_trials = 300, p_saf = 0.25,
                          p_tox = 0.35, seed = 3)
  expect_equal(a$settings$p_saf, 0.18)
  expect_equal(b$settings$p_saf, 0.25)
  expect_gt(b$boundary$lambda_e, a$boundary$lambda_e)
  expect_lt(b$boundary$lambda_d, a$boundary$lambda_d)
  expect_false(identical(a$n_pts, b$n_pts))
})

test_that("waterfall_simulate validates its arguments", {
  p_true <- zhang2016_p_true(1)
  expect_error(
    waterfall_simulate(target = 0.30, p_true = t(p_true), cohort_size = 3,
                       n_trials = 2),
    "more rows than columns"
  )
  expect_error(
    waterfall_simulate(target = 0.30, p_true = p_true[1L, , drop = FALSE],
                       cohort_size = 3, n_trials = 2),
    "at least two rows"
  )
  expect_error(
    waterfall_simulate(target = 0.30, p_true = p_true, cohort_size = 3,
                       n_cohort = 6, n_trials = 2),
    "'n_cohort' must hold one whole number"
  )
  expect_error(
    waterfall_simulate(target = 0.30, p_true = p_true, cohort_size = 0,
                       n_trials = 2),
    "'cohort_size' must be at least 1"
  )
  expect_error(
    waterfall_simulate(target = 0.30, p_true = p_true, cohort_size = 3,
                       n_trials = 2, boin_compatible = NA),
    "'boin_compatible' must be TRUE or FALSE"
  )
  expect_warning(
    waterfall_simulate(target = 0.30, p_true = p_true, cohort_size = 3,
                       n_trials = 2, n_earlystop = 6),
    "'n_earlystop' is low"
  )
})
