jss_p_true <- function() {
  matrix(c(0.02, 0.04, 0.08, 0.14,
           0.08, 0.25, 0.42, 0.48,
           0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE)
}

test_that("comb_boin_simulate returns a well formed object", {
  trials <- comb_boin_simulate(
    target = 0.25, p_true = jss_p_true(), n_cohort = 10, cohort_size = 3,
    n_trials = 200, seed = 1
  )

  expect_s3_class(trials, "comb_boin_trials")
  expect_equal(dim(trials$n_pts), c(200L, 3L, 4L))
  expect_equal(dim(trials$n_tox), c(200L, 3L, 4L))
  expect_equal(dim(trials$eliminated), c(200L, 3L, 4L))
  expect_equal(dimnames(trials$n_pts)[[2L]], c("A1", "A2", "A3"))
  expect_equal(dimnames(trials$n_pts)[[3L]], c("B1", "B2", "B3", "B4"))
  expect_length(trials$cohorts_used, 200L)
  expect_length(trials$stop_reason, 200L)
  expect_type(trials$n_pts, "integer")
  expect_type(trials$eliminated, "logical")
})

test_that("simulated counts respect the design constraints", {
  trials <- comb_boin_simulate(
    target = 0.25, p_true = jss_p_true(), n_cohort = 10, cohort_size = 3,
    n_trials = 500, seed = 2
  )

  totals <- apply(trials$n_pts, 1L, sum)
  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(totals <= 30))
  expect_true(all(totals > 0))
  expect_true(all(trials$cohorts_used >= 1 & trials$cohorts_used <= 10))
  expect_true(all(trials$stop_reason %in%
                    c("lowest_dose_eliminated", "lowest_dose_too_toxic",
                      "n_earlystop", "max_cohorts")))
})

test_that("the trials agree with a separate implementation in Python", {
  # Reference values from a port of BOIN::get.oc.comb() to Python, run on a
  # replica of R's Mersenne-Twister stream; the first four trials of each
  # setting, with the counts of each trial in column-major order.
  plain <- comb_boin_simulate(
    target = 0.25, p_true = jss_p_true(), n_cohort = 16, cohort_size = 3,
    n_earlystop = 100, n_trials = 4, seed = 6
  )
  expect_equal(as.integer(plain$n_pts[1, , ]),
               c(3, 6, 0, 3, 24, 0, 0, 12, 0, 0, 0, 0))
  expect_equal(as.integer(plain$n_tox[1, , ]),
               c(0, 0, 0, 0, 4, 0, 0, 7, 0, 0, 0, 0))
  expect_equal(as.integer(plain$n_pts[2, , ]),
               c(3, 6, 0, 9, 27, 3, 0, 0, 0, 0, 0, 0))
  expect_equal(as.integer(plain$n_tox[2, , ]),
               c(0, 0, 0, 0, 9, 2, 0, 0, 0, 0, 0, 0))
  expect_equal(as.integer(plain$n_pts[3, , ]),
               c(3, 6, 0, 6, 27, 3, 0, 3, 0, 0, 0, 0))
  expect_equal(as.integer(plain$n_tox[3, , ]),
               c(0, 0, 0, 0, 5, 2, 0, 2, 0, 0, 0, 0))
  expect_equal(as.integer(plain$n_pts[4, , ]),
               c(3, 12, 0, 6, 27, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(as.integer(plain$n_tox[4, , ]),
               c(0, 1, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(plain$stop_reason, rep("max_cohorts", 4L))

  titrated <- comb_boin_simulate(
    target = 0.25, p_true = jss_p_true(), n_cohort = 12, cohort_size = 3,
    n_earlystop = 18, titration = TRUE, n_trials = 4, seed = 11
  )
  expect_equal(as.integer(titrated$n_pts[1, , ]),
               c(1, 4, 3, 0, 18, 3, 0, 3, 0, 0, 0, 0))
  expect_equal(as.integer(titrated$n_tox[1, , ]),
               c(0, 0, 2, 0, 5, 2, 0, 2, 0, 0, 0, 0))
  expect_equal(as.integer(titrated$n_pts[2, , ]),
               c(1, 7, 0, 6, 13, 0, 0, 7, 4, 0, 0, 3))
  expect_equal(as.integer(titrated$n_tox[2, , ]),
               c(0, 0, 0, 1, 5, 0, 0, 3, 3, 0, 0, 2))
  expect_equal(as.integer(titrated$n_pts[3, , ]),
               c(1, 1, 19, 0, 0, 13, 0, 0, 3, 0, 0, 0))
  expect_equal(as.integer(titrated$n_tox[3, , ]),
               c(0, 0, 4, 0, 0, 5, 0, 0, 3, 0, 0, 0))
  expect_equal(as.integer(titrated$n_pts[4, , ]),
               c(1, 3, 0, 7, 13, 4, 0, 0, 10, 0, 0, 3))
  expect_equal(as.integer(titrated$n_tox[4, , ]),
               c(0, 0, 0, 0, 4, 2, 0, 0, 3, 0, 0, 2))
  expect_equal(titrated$stop_reason,
               c("n_earlystop", "max_cohorts", "n_earlystop", "max_cohorts"))
  expect_equal(titrated$cohorts_used, c(10L, 12L, 11L, 12L))
})

test_that("titration adds one patient per combination on the way up", {
  # With no DLT anywhere, titration climbs from A1B1 to the highest
  # combination, one patient at each of the J + K - 1 combinations on the way,
  # and every planned cohort is then treated there.
  p_true <- matrix(0, nrow = 3, ncol = 4)
  trials <- comb_boin_simulate(
    target = 0.25, p_true = p_true, n_cohort = 10, cohort_size = 3,
    n_earlystop = 100, titration = TRUE, n_trials = 50, seed = 3
  )

  totals <- apply(trials$n_pts, 1L, sum)
  expect_equal(totals, rep(3 + 4 - 2 + 10 * 3, 50))
  expect_true(all(trials$n_pts[, 3, 4] == 30L))
  expect_true(all(trials$n_tox == 0L))
})

test_that("without titration every trial uses all cohorts when nothing stops it", {
  p_true <- matrix(c(0.01, 0.02, 0.03,
                     0.02, 0.03, 0.04), nrow = 2, byrow = TRUE)
  cohort_size <- c(1, 1, 3, 3, 3, 3)
  trials <- comb_boin_simulate(
    target = 0.30, p_true = p_true, n_cohort = 6, cohort_size = cohort_size,
    n_earlystop = 100, n_trials = 100, seed = 4
  )

  expect_equal(trials$settings$max_total_pts, sum(cohort_size))
  expect_equal(apply(trials$n_pts, 1L, sum), rep(sum(cohort_size), 100))
  expect_true(all(trials$stop_reason == "max_cohorts"))
})

test_that("a toxic lowest combination stops the trial", {
  p_true <- matrix(c(0.90, 0.95,
                     0.95, 0.99), nrow = 2, byrow = TRUE)
  trials <- comb_boin_simulate(
    target = 0.30, p_true = p_true, n_cohort = 10, cohort_size = 3,
    n_trials = 100, seed = 5
  )

  expect_true(mean(trials$stop_reason == "lowest_dose_eliminated") > 0.9)
  stopped <- trials$stop_reason == "lowest_dose_eliminated"
  expect_true(all(trials$eliminated[stopped, 1, 1]))
  expect_true(all(trials$eliminated[stopped, 2, 2]))
})

test_that("the trial starts at start_dose", {
  trials <- comb_boin_simulate(
    target = 0.25, p_true = jss_p_true(), n_cohort = 1, cohort_size = 3,
    start_dose = c(2, 3), n_trials = 20, seed = 6
  )

  expect_true(all(trials$n_pts[, 2, 3] == 3L))
  expect_equal(apply(trials$n_pts, 1L, sum), rep(3L, 20))
})

test_that("the same seed reproduces the same trials", {
  args <- list(target = 0.25, p_true = jss_p_true(), n_cohort = 10,
               cohort_size = 3, n_trials = 100)

  first <- do.call(comb_boin_simulate, c(args, list(seed = 7)))
  second <- do.call(comb_boin_simulate, c(args, list(seed = 7)))
  other <- do.call(comb_boin_simulate, c(args, list(seed = 8)))

  expect_identical(first$n_pts, second$n_pts)
  expect_identical(first$n_tox, second$n_tox)
  expect_false(identical(first$n_tox, other$n_tox))
})

test_that("the state of the random number generator is restored", {
  set.seed(99)
  before <- .Random.seed
  invisible(comb_boin_simulate(
    target = 0.25, p_true = jss_p_true(), n_cohort = 5, cohort_size = 3,
    n_trials = 50, seed = 1
  ))
  expect_identical(.Random.seed, before)
})

test_that("invalid arguments are rejected", {
  expect_error(
    comb_boin_simulate(target = 0.25, p_true = c(0.1, 0.2, 0.3),
                       n_cohort = 10, cohort_size = 3),
    "p_true"
  )
  expect_error(
    comb_boin_simulate(target = 0.25, p_true = matrix(c(0.1, 1.2), 1),
                       n_cohort = 10, cohort_size = 3),
    "between 0 and 1"
  )
  expect_error(
    comb_boin_simulate(target = 0.25, p_true = jss_p_true(), n_cohort = 10,
                       cohort_size = 3, start_dose = 1),
    "start_dose"
  )
  expect_error(
    comb_boin_simulate(target = 0.25, p_true = jss_p_true(), n_cohort = 10,
                       cohort_size = 3, start_dose = c(4, 1)),
    "start_dose"
  )
  expect_warning(
    comb_boin_simulate(target = 0.25,
                       p_true = matrix(c(0.30, 0.20, 0.40, 0.50), 2),
                       n_cohort = 2, cohort_size = 3, n_trials = 5),
    "non-decreasing"
  )
})
