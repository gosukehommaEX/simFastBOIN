# The expected selections were obtained independently with a port of the
# selection step of BOIN::get.oc.comb() to Python; in each case the selected
# combination is clear of its competitors, so the result does not depend on
# the last digits of the isotonic regression.

test_that("the combination closest to the target is selected", {
  n_pts <- matrix(c(3, 6, 3,
                    3, 3, 0), nrow = 2, byrow = TRUE)
  n_tox <- matrix(c(0, 1, 1,
                    0, 2, 0), nrow = 2, byrow = TRUE)

  res <- comb_boin_select_mtd(n_pts, n_tox, target = 0.30)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_equal(res$mtd_a, 1L)
  expect_equal(res$mtd_b, 3L)
  expect_equal(res$reason, "selected")
})

test_that("no combination is selected when the lowest one is eliminated", {
  n_pts <- matrix(c(3, 0,
                    0, 0), nrow = 2, byrow = TRUE)
  n_tox <- matrix(c(3, 0,
                    0, 0), nrow = 2, byrow = TRUE)

  res <- comb_boin_select_mtd(n_pts, n_tox, target = 0.30)

  expect_true(is.na(res$mtd_a))
  expect_true(is.na(res$mtd_b))
  expect_equal(res$reason, "lowest_dose_eliminated")
})

test_that("bound_mtd and extrasafe can rule out every combination", {
  # Two DLTs out of three: the isotonic estimate is about 0.65, above the
  # de-escalation boundary 0.358, and Pr(p > 0.30) = 0.916 exceeds the extra
  # safety cutoff 0.95 - 0.05 but not the elimination cutoff 0.95.
  n_pts <- matrix(c(3, 0,
                    0, 0), nrow = 2, byrow = TRUE)
  n_tox <- matrix(c(2, 0,
                    0, 0), nrow = 2, byrow = TRUE)

  plain <- comb_boin_select_mtd(n_pts, n_tox, target = 0.30)
  expect_equal(c(plain$mtd_a, plain$mtd_b), c(1L, 1L))

  bounded <- comb_boin_select_mtd(n_pts, n_tox, target = 0.30,
                                  bound_mtd = TRUE)
  expect_true(is.na(bounded$mtd_a))
  expect_equal(bounded$reason, "no_dose_below_bound")

  safe <- comb_boin_select_mtd(n_pts, n_tox, target = 0.30, extrasafe = TRUE)
  expect_true(is.na(safe$mtd_a))
  expect_equal(safe$reason, "lowest_dose_eliminated")
})

test_that("elimination at selection takes the rest of the row and the column", {
  n_pts <- as_comb_count_array(matrix(3, nrow = 2, ncol = 3), "n_pts")

  one <- as_comb_count_array(matrix(c(0, 3, 0,
                                      0, 0, 0), nrow = 2, byrow = TRUE), "n_tox")
  elim <- comb_boin_eliminate_cpp(n_pts, one, 0.30, 0.95, FALSE, 0.05)
  expect_equal(elim[1, , ], matrix(c(FALSE, TRUE, TRUE,
                                     FALSE, TRUE, FALSE), nrow = 2, byrow = TRUE))

  # A2B3 survives although A1B2 is eliminated, and, pooled with its
  # neighbors, it is the combination selected, as in the reference
  # implementation.
  res <- comb_boin_select_mtd(matrix(3, nrow = 2, ncol = 3),
                              matrix(c(0, 3, 0,
                                       0, 0, 0), nrow = 2, byrow = TRUE),
                              target = 0.30)
  expect_equal(c(res$mtd_a, res$mtd_b), c(2L, 3L))

  # The scan of each row stops at its first eliminated combination.
  two <- as_comb_count_array(matrix(c(0, 3, 0,
                                      3, 0, 0), nrow = 2, byrow = TRUE), "n_tox")
  elim <- comb_boin_eliminate_cpp(n_pts, two, 0.30, 0.95, FALSE, 0.05)
  expect_equal(elim[1, , ], matrix(c(FALSE, TRUE, TRUE,
                                     TRUE, TRUE, TRUE), nrow = 2, byrow = TRUE))
})

test_that("an array of trials gives the same selections as single trials", {
  trials <- comb_boin_simulate(
    target = 0.30,
    p_true = matrix(c(0.05, 0.15, 0.30,
                      0.15, 0.30, 0.45), nrow = 2, byrow = TRUE),
    n_cohort = 8, cohort_size = 3, n_trials = 20, seed = 10
  )

  together <- comb_boin_select_mtd(trials$n_pts, trials$n_tox, target = 0.30)
  expect_equal(nrow(together), 20L)

  for (t in c(1, 7, 20)) {
    alone <- comb_boin_select_mtd(trials$n_pts[t, , ], trials$n_tox[t, , ],
                                  target = 0.30)
    expect_equal(alone$mtd_a, together$mtd_a[t])
    expect_equal(alone$mtd_b, together$mtd_b[t])
  }
})

test_that("invalid counts are rejected", {
  expect_error(
    comb_boin_select_mtd(matrix(3, 2, 2), matrix(0, 2, 3), target = 0.30),
    "same dimensions"
  )
  expect_error(
    comb_boin_select_mtd(matrix(3, 2, 2), matrix(4, 2, 2), target = 0.30),
    "must not exceed"
  )
  expect_error(
    comb_boin_select_mtd(c(3, 3), c(0, 0), target = 0.30),
    "matrix for one trial"
  )
})
