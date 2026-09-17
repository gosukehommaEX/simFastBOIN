test_that("the renamed functions still work but warn", {
  expect_warning(result <- get_boin_boundary(target = 0.30), "deprecated")
  expect_equal(result, boin_lambda(target = 0.30))

  expect_warning(
    decisions <- get_boin_decision(target = 0.30, max_n = 9), "deprecated"
  )
  expect_equal(decisions, boin_decision_table(target = 0.30, max_n = 9))

  expect_warning(
    bounds <- get_boin_stopping_boundaries(target = 0.30, max_n = 9,
                                           extrasafe = TRUE),
    "deprecated"
  )
  expect_equal(bounds, boin_boundary(target = 0.30, max_n = 9, extrasafe = TRUE))

  expect_warning(
    est <- isotonic_regression(n_pts = c(3, 6, 9), n_tox = c(0, 1, 3)),
    "deprecated"
  )
  expect_equal(est, boin_isotonic(n_pts = c(3, 6, 9), n_tox = c(0, 1, 3)))

  expect_warning(
    mtd <- select_mtd(n_pts = c(3, 6, 9), n_tox = c(0, 1, 3), target = 0.30),
    "deprecated"
  )
  expect_equal(mtd, boin_select_mtd(n_pts = c(3, 6, 9), n_tox = c(0, 1, 3),
                                    target = 0.30))
})

test_that("get_pts_and_tox forwards to boin_simulate", {
  expect_warning(
    trials <- get_pts_and_tox(target = 0.30, p_true = c(0.10, 0.25, 0.40),
                              n_cohort = 5, cohort_size = 3, n_trials = 20,
                              seed = 1),
    "deprecated"
  )
  expect_s3_class(trials, "boin_trials")
})
