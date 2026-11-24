test_that("waterfall_select_mtd selects one MTD per level of drug A", {
  # The estimates (y + 0.05) / (n + 0.1) are 0.016, 0.016, 0.335 in the first
  # row and 0.172, 0.335 in the second, whose third combination is excluded.
  # The expected selections were also computed in Python.
  n_pts <- matrix(c(3, 3, 9,
                    6, 9, 3), nrow = 2, byrow = TRUE)
  n_tox <- matrix(c(0, 0, 3,
                    1, 3, 2), nrow = 2, byrow = TRUE)
  eliminated <- matrix(FALSE, nrow = 2, ncol = 3)
  eliminated[2, 3] <- TRUE

  sel <- waterfall_select_mtd(n_pts, n_tox, eliminated, target = 0.30)

  expect_s3_class(sel, "data.frame")
  expect_named(sel, c("trial", "dose_a", "dose_b"))
  expect_equal(sel$trial, c(1L, 1L))
  expect_equal(sel$dose_a, 1:2)
  expect_equal(sel$dose_b, c(3L, 2L))
})

test_that("an MTD to the left of the MTD of the row above is moved", {
  # Only the first combination of the first row was treated, so it is the
  # closest to the target there; the MTD of the second row is at the second
  # level of drug B, and the first row takes that level.
  n_pts <- matrix(c(3, 0, 0,
                    6, 9, 3), nrow = 2, byrow = TRUE)
  n_tox <- matrix(c(0, 0, 0,
                    1, 3, 2), nrow = 2, byrow = TRUE)
  eliminated <- matrix(FALSE, nrow = 2, ncol = 3)
  eliminated[2, 3] <- TRUE

  sel <- waterfall_select_mtd(n_pts, n_tox, eliminated, target = 0.30)
  expect_equal(sel$dose_b, c(2L, 2L))
})

test_that("rows without admissible data have no MTD", {
  # The second row is excluded from its first combination on.
  n_pts <- matrix(c(3, 6, 0,
                    3, 0, 0), nrow = 2, byrow = TRUE)
  n_tox <- matrix(c(0, 2, 0,
                    3, 0, 0), nrow = 2, byrow = TRUE)
  eliminated <- matrix(c(FALSE, FALSE, FALSE,
                         TRUE, TRUE, TRUE), nrow = 2, byrow = TRUE)

  sel <- waterfall_select_mtd(n_pts, n_tox, eliminated, target = 0.30)
  expect_equal(sel$dose_b, c(2L, NA))

  # A trial without patients selects nothing.
  empty <- matrix(0, nrow = 2, ncol = 3)
  sel <- waterfall_select_mtd(empty, empty, empty == 1, target = 0.30)
  expect_equal(sel$dose_b, c(NA_integer_, NA_integer_))
})

test_that("several trials are selected as each trial alone", {
  trials <- waterfall_simulate(
    target = 0.30, p_true = zhang2016_p_true(12), cohort_size = 3,
    n_trials = 30, seed = 4
  )
  all_trials <- waterfall_select_mtd(trials$n_pts, trials$n_tox,
                                     trials$eliminated, target = 0.30)
  expect_equal(nrow(all_trials), 30L * 3L)
  expect_equal(all_trials$trial, rep(1:30, each = 3L))
  expect_equal(all_trials$dose_a, rep(1:3, times = 30L))

  for (t in c(1L, 7L, 30L)) {
    one <- waterfall_select_mtd(trials$n_pts[t, , ], trials$n_tox[t, , ],
                                trials$eliminated[t, , ], target = 0.30)
    expect_equal(one$dose_b, all_trials$dose_b[all_trials$trial == t])
  }

  # The MTD contour does not turn back: the level of drug B does not
  # increase with the level of drug A.
  mtd <- matrix(all_trials$dose_b, nrow = 30L, byrow = TRUE)
  both <- !is.na(mtd[, 1L]) & !is.na(mtd[, 2L])
  expect_true(all(mtd[both, 1L] >= mtd[both, 2L]))
  both <- !is.na(mtd[, 2L]) & !is.na(mtd[, 3L])
  expect_true(all(mtd[both, 2L] >= mtd[both, 3L]))
})

test_that("waterfall_select_mtd validates its arguments", {
  n_pts <- matrix(3, nrow = 2, ncol = 3)
  n_tox <- matrix(1, nrow = 2, ncol = 3)
  eliminated <- matrix(FALSE, nrow = 2, ncol = 3)

  expect_error(
    waterfall_select_mtd(n_pts, n_tox[, 1:2], eliminated, target = 0.30),
    "must have the same dimensions"
  )
  expect_error(
    waterfall_select_mtd(n_pts, n_pts + 1, eliminated, target = 0.30),
    "'n_tox' must not exceed 'n_pts'"
  )
  expect_error(
    waterfall_select_mtd(n_pts, n_tox, eliminated * 1, target = 0.30),
    "'eliminated' must be logical"
  )
  expect_error(
    waterfall_select_mtd(n_pts, n_tox, eliminated[, 1:2], target = 0.30),
    "'eliminated' must have the same dimensions"
  )
  expect_error(
    waterfall_select_mtd(n_pts, n_tox, eliminated, target = 1.5),
    "'target' must lie strictly between 0 and 1"
  )
})
