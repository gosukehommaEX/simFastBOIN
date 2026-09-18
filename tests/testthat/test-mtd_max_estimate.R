test_that("an explicit cap selects the dose the estimates allow", {
  # Isotonic estimates are 0.0082, 0.1721 and 0.5000 for these data.
  n_pts <- c(6, 6, 6)
  n_tox <- c(0, 1, 3)

  expect_equal(boin_isotonic(n_pts, n_tox)[1, 3], 0.5)

  expect_equal(boin_select_mtd(n_pts, n_tox, target = 0.30)$mtd, 2L)
  expect_equal(
    boin_select_mtd(n_pts, n_tox, target = 0.30, mtd_max_estimate = 0.30)$mtd, 2L
  )
  expect_equal(
    boin_select_mtd(n_pts, n_tox, target = 0.30, mtd_max_estimate = 0.10)$mtd, 1L
  )

  none <- boin_select_mtd(n_pts, n_tox, target = 0.30, mtd_max_estimate = 0.005)
  expect_true(is.na(none$mtd))
  expect_equal(none$reason, "no_dose_below_bound")
})

test_that("a cap above every estimate behaves like no cap", {
  n_pts <- matrix(c(3, 6, 9, 6,
                    6, 6, 6, 6), nrow = 2, byrow = TRUE)
  n_tox <- matrix(c(0, 1, 3, 4,
                    1, 2, 4, 5), nrow = 2, byrow = TRUE)

  expect_equal(
    boin_select_mtd(n_pts, n_tox, target = 0.30, mtd_max_estimate = 0.999)$mtd,
    boin_select_mtd(n_pts, n_tox, target = 0.30)$mtd
  )
})

test_that("the cap can be set at or below the target, which bound_mtd cannot", {
  # The de-escalation boundary always lies above the target, so bound_mtd can
  # never impose a cap of 0.30 when the target is 0.30.
  expect_gt(boin_lambda(target = 0.30)$lambda_d, 0.30)

  args <- list(target = 0.30, p_true = c(0.10, 0.20, 0.30, 0.42, 0.55),
               n_cohort = 20, cohort_size = 3, n_trials = 1000,
               keep_trials = TRUE, seed = 41)

  free <- do.call(sim_boin, args)
  by_lambda <- do.call(sim_boin, c(args, list(bound_mtd = TRUE)))
  by_value <- do.call(sim_boin, c(args, list(mtd_max_estimate = 0.30)))

  # The trials themselves are untouched, only the selection differs.
  expect_identical(free$trials$n_pts, by_value$trials$n_pts)
  expect_identical(free$trials$n_tox, by_value$trials$n_tox)

  # A tighter cap never selects a higher dose.
  both <- !is.na(by_lambda$trials$mtd) & !is.na(by_value$trials$mtd)
  expect_true(all(by_value$trials$mtd[both] <= by_lambda$trials$mtd[both]))
  expect_gte(by_value$percent_no_mtd, by_lambda$percent_no_mtd)
  expect_gte(by_lambda$percent_no_mtd, free$percent_no_mtd)
})

test_that("the cap overrides bound_mtd", {
  n_pts <- c(6, 6, 6)
  n_tox <- c(0, 1, 3)

  expect_equal(
    boin_select_mtd(n_pts, n_tox, target = 0.30, bound_mtd = TRUE,
                    mtd_max_estimate = 0.10)$mtd,
    1L
  )
  expect_equal(
    boin_select_mtd(n_pts, n_tox, target = 0.30, bound_mtd = FALSE,
                    mtd_max_estimate = 0.10)$mtd,
    1L
  )
})

test_that("mtd_max_estimate is validated", {
  n_pts <- c(6, 6, 6)
  n_tox <- c(0, 1, 3)

  expect_error(
    boin_select_mtd(n_pts, n_tox, target = 0.30, mtd_max_estimate = 0),
    "strictly between"
  )
  expect_error(
    boin_select_mtd(n_pts, n_tox, target = 0.30, mtd_max_estimate = "0.3"),
    "single finite number"
  )
})
