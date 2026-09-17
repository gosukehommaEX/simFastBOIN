test_that("boin_simulate returns a well formed object", {
  trials <- boin_simulate(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 1
  )

  expect_s3_class(trials, "boin_trials")
  expect_equal(dim(trials$n_pts), c(200L, 5L))
  expect_equal(dim(trials$n_tox), c(200L, 5L))
  expect_equal(dim(trials$eliminated), c(200L, 5L))
  expect_length(trials$cohorts_used, 200L)
  expect_length(trials$stop_reason, 200L)
  expect_type(trials$n_pts, "integer")
  expect_type(trials$eliminated, "logical")
})

test_that("simulated counts respect the design constraints", {
  trials <- boin_simulate(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, n_trials = 500, seed = 2
  )

  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(rowSums(trials$n_pts) <= 30))
  expect_true(all(rowSums(trials$n_pts) > 0))
  expect_true(all(trials$cohorts_used >= 1 & trials$cohorts_used <= 10))
  expect_true(all(trials$stop_reason %in%
                    c("lowest_dose_eliminated", "lowest_dose_too_toxic",
                      "n_earlystop", "max_sample_size", "max_cohorts")))
})

test_that("a vector cohort_size sets the maximum sample size to its sum", {
  cohort_size <- c(1, 1, 3, 3, 3, 3)
  trials <- boin_simulate(
    target = 0.30, p_true = c(0.10, 0.25, 0.40, 0.55),
    n_cohort = 6, cohort_size = cohort_size, n_trials = 300, seed = 3
  )

  expect_equal(trials$settings$max_total_pts, sum(cohort_size))
  expect_true(all(rowSums(trials$n_pts) <= sum(cohort_size)))
  expect_true(max(rowSums(trials$n_pts)) == sum(cohort_size))
})

test_that("a vector cohort_size works together with titration", {
  expect_silent(
    boin_simulate(
      target = 0.30, p_true = c(0.10, 0.25, 0.40, 0.55),
      n_cohort = 6, cohort_size = c(3, 3, 3, 3, 3, 3), n_trials = 50,
      titration = TRUE, seed = 4
    )
  )
})

test_that("the same seed reproduces the same trials", {
  args <- list(target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
               n_cohort = 10, cohort_size = 3, n_trials = 100)

  first <- do.call(boin_simulate, c(args, list(seed = 7)))
  second <- do.call(boin_simulate, c(args, list(seed = 7)))
  other <- do.call(boin_simulate, c(args, list(seed = 8)))

  expect_identical(first$n_pts, second$n_pts)
  expect_identical(first$n_tox, second$n_tox)
  expect_false(identical(first$n_tox, other$n_tox))
})

test_that("the state of the random number generator is restored", {
  set.seed(99)
  before <- .Random.seed
  invisible(boin_simulate(
    target = 0.30, p_true = c(0.10, 0.25, 0.40),
    n_cohort = 5, cohort_size = 3, n_trials = 50, seed = 1
  ))
  expect_identical(.Random.seed, before)
})

test_that("seed = NULL draws from the current stream", {
  set.seed(11)
  first <- boin_simulate(
    target = 0.30, p_true = c(0.10, 0.25, 0.40),
    n_cohort = 5, cohort_size = 3, n_trials = 20, seed = NULL
  )
  set.seed(11)
  second <- boin_simulate(
    target = 0.30, p_true = c(0.10, 0.25, 0.40),
    n_cohort = 5, cohort_size = 3, n_trials = 20, seed = NULL
  )
  expect_identical(first$n_tox, second$n_tox)
})

test_that("titration treats one patient per dose until the first DLT", {
  # With a completely safe dose-toxicity curve every trial escalates through the
  # whole range on single patients before the first full cohort.
  trials <- boin_simulate(
    target = 0.30, p_true = rep(0, 5),
    n_cohort = 10, cohort_size = 3, n_trials = 20,
    titration = TRUE, seed = 5
  )

  expect_true(all(trials$n_tox == 0))
  expect_true(all(trials$n_pts[, 1] == 1))
  expect_true(all(trials$n_pts[, 5] >= 3))
})

test_that("titration is switched off when the cohort size is one", {
  trials <- boin_simulate(
    target = 0.30, p_true = c(0.10, 0.25, 0.40),
    n_cohort = 12, cohort_size = 1, n_trials = 20, titration = TRUE, seed = 6
  )
  expect_false(trials$settings$titration)
})

test_that("start_dose moves where the first cohort is treated", {
  # A single cohort, so every trial consists of the starting dose alone.
  one_cohort <- boin_simulate(
    target = 0.30, p_true = c(0.02, 0.05, 0.10, 0.25, 0.40),
    n_cohort = 1, cohort_size = 3, n_trials = 50, start_dose = 3, seed = 8
  )
  expect_true(all(one_cohort$n_pts[, c(1, 2, 4, 5)] == 0))
  expect_true(all(one_cohort$n_pts[, 3] == 3))

  # With two cohorts a trial can de-escalate once, to dose 2 but no further, so
  # the doses below the starting dose minus one are still never reached.
  two_cohorts <- boin_simulate(
    target = 0.30, p_true = c(0.02, 0.05, 0.10, 0.25, 0.40),
    n_cohort = 2, cohort_size = 3, n_trials = 200, start_dose = 3, seed = 8
  )
  expect_true(all(two_cohorts$n_pts[, 1] == 0))
  expect_true(all(two_cohorts$n_pts[, 3] >= 3))
})

test_that("the simple early stopping rule bounds the sample size at a dose", {
  n_earlystop <- 9
  cohort_size <- 3
  trials <- boin_simulate(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 20, cohort_size = cohort_size, n_trials = 300,
    n_earlystop = n_earlystop, n_earlystop_rule = "simple", seed = 9
  )
  expect_true(max(trials$n_pts) <= n_earlystop + cohort_size - 1)
})

test_that("boin_simulate validates its arguments", {
  base <- list(target = 0.30, p_true = c(0.10, 0.25, 0.40),
               n_cohort = 5, cohort_size = 3, n_trials = 10)

  expect_error(do.call(boin_simulate, modifyList(base, list(p_true = c(0.5, 1.5)))),
               "between 0 and 1")
  expect_error(do.call(boin_simulate, modifyList(base, list(n_cohort = 0))),
               "at least 1")
  expect_error(do.call(boin_simulate, modifyList(base, list(cohort_size = 0))),
               "at least 1")
  expect_error(do.call(boin_simulate, modifyList(base, list(start_dose = 4))),
               "must not exceed")
  expect_warning(do.call(boin_simulate, modifyList(base, list(n_earlystop = 3))),
                 "recommended")
  expect_warning(do.call(boin_simulate, modifyList(base, list(p_true = c(0.4, 0.2, 0.3)))),
                 "non-decreasing")
})
