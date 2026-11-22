test_that("gboin_select_mtd selects the doses of the rule", {
  n_pts <- rbind(c(3, 6, 12, 9, 0, 0), c(3, 3, 6, 12, 6, 0),
                 c(6, 9, 9, 6, 0, 0), c(3, 3, 3, 3, 6, 12),
                 c(6, 0, 0, 0, 0, 0), c(9, 12, 9, 0, 0, 0))
  ets_total <- rbind(c(0.5, 1.5, 5.5, 6.0, 0, 0), c(0, 0.5, 1.5, 5.0, 5.5, 0),
                     c(1.0, 4.0, 5.5, 6.5, 0, 0), c(0, 0, 0.5, 0.5, 1.5, 4.0),
                     c(7.5, 0, 0, 0, 0, 0), c(3.0, 7.5, 9.0, 0, 0, 0))

  # Expected doses computed apart from the package (Python, scipy)
  res <- gboin_select_mtd(n_pts, ets_total / 1.5, target = 0.47 / 1.5)
  expect_identical(res$mtd, c(3L, 4L, 2L, 6L, NA, 1L))
  expect_identical(res$reason,
                   c(rep("selected", 4), "lowest_dose_eliminated",
                     "selected"))
  expect_identical(res$trial, 1:6)

  # A single trial as vectors
  one <- gboin_select_mtd(n_pts[2, ], ets_total[2, ] / 1.5, 0.47 / 1.5)
  expect_identical(one$mtd, 4L)
})

test_that("whole numbers of DLTs give the selection of the BOIN design", {
  trials <- boin_simulate(
    target = 0.30, p_true = c(0.08, 0.15, 0.31, 0.45, 0.55), n_cohort = 10,
    cohort_size = 3, n_trials = 500, n_earlystop = 12, seed = 77
  )
  for (xs in c(FALSE, TRUE)) {
    expect_identical(
      gboin_select_mtd(trials$n_pts, trials$n_tox, 0.30, extrasafe = xs),
      boin_select_mtd(trials$n_pts, trials$n_tox, 0.30, extrasafe = xs)
    )
  }
  expect_identical(
    gboin_select_mtd(trials$n_pts, trials$n_tox, 0.30, min_mtd_sample = 6),
    boin_select_mtd(trials$n_pts, trials$n_tox, 0.30, min_mtd_sample = 6)
  )
})

test_that("gboin_select_mtd rejects inconsistent data", {
  expect_error(gboin_select_mtd(c(3, 3, 3), c(3.5, 0, 0), 0.3),
               "must not exceed")
  expect_error(gboin_select_mtd(c(3, 3, 3), c(1, 0), 0.3),
               "same dimensions")
  expect_error(gboin_select_mtd(c(3, 3, 3), c(-1, 0, 0), 0.3),
               "non-negative")
  expect_error(gboin_select_mtd(c(3, 3, 3), c(1, 0, 0), 1.3), "target")
  none <- gboin_select_mtd(c(0, 0, 0), c(0, 0, 0), 0.3)
  expect_identical(none$mtd, NA_integer_)
  expect_identical(none$reason, "no_admissible_dose")
})
