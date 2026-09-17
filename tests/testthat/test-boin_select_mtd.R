test_that("boin_select_mtd picks the dose closest to the target", {
  n_pts <- matrix(c(3, 6, 9, 3,
                    3, 6, 9, 3), nrow = 2, byrow = TRUE)
  n_tox <- matrix(c(0, 1, 3, 2,
                    0, 1, 2, 1), nrow = 2, byrow = TRUE)

  out <- boin_select_mtd(n_pts, n_tox, target = 0.30)

  expect_s3_class(out, "data.frame")
  expect_equal(names(out), c("trial", "mtd", "reason"))
  expect_equal(out$mtd, c(3L, 4L))
  expect_equal(out$reason, c("selected", "selected"))
})

test_that("no dose is selected once the lowest dose is eliminated", {
  # Nine DLTs out of nine patients at the lowest dose.
  out <- boin_select_mtd(c(9, 3, 0), c(9, 3, 0), target = 0.30)

  expect_true(is.na(out$mtd))
  expect_equal(out$reason, "lowest_dose_eliminated")
})

test_that("untreated doses cannot be selected", {
  out <- boin_select_mtd(c(3, 3, 0, 0), c(1, 1, 0, 0), target = 0.30)
  expect_equal(out$mtd, 1L)

  none <- boin_select_mtd(c(0, 0, 0), c(0, 0, 0), target = 0.30)
  expect_true(is.na(none$mtd))
  expect_equal(none$reason, "no_admissible_dose")
})

test_that("min_mtd_sample removes sparsely treated doses", {
  n_pts <- matrix(c(12, 12, 3), nrow = 1)
  n_tox <- matrix(c(0, 1, 1), nrow = 1)

  relaxed <- boin_select_mtd(n_pts, n_tox, target = 0.30, min_mtd_sample = 1)
  strict <- boin_select_mtd(n_pts, n_tox, target = 0.30, min_mtd_sample = 6)

  expect_equal(relaxed$mtd, 3L)
  expect_equal(strict$mtd, 2L)
})

test_that("bound_mtd refuses doses above the de-escalation boundary", {
  # Every treated dose is well above lambda_d but none is eliminated.
  n_pts <- matrix(c(3, 3), nrow = 1)
  n_tox <- matrix(c(2, 2), nrow = 1)

  free <- boin_select_mtd(n_pts, n_tox, target = 0.30, bound_mtd = FALSE)
  bounded <- boin_select_mtd(n_pts, n_tox, target = 0.30, bound_mtd = TRUE)

  expect_false(is.na(free$mtd))
  expect_true(is.na(bounded$mtd))
  expect_equal(bounded$reason, "no_dose_below_lambda_d")
})

test_that("extrasafe can withhold the MTD when the lowest dose looks toxic", {
  # Pr(p > 0.30 | 6 of 12) is 0.938, which clears the relaxed safety cutoff of
  # 0.90 but not the elimination cutoff of 0.95.
  n_pts <- matrix(c(12, 12, 12), nrow = 1)
  n_tox <- matrix(c(6, 6, 6), nrow = 1)

  expect_equal(
    boin_select_mtd(n_pts, n_tox, target = 0.30, extrasafe = FALSE)$mtd, 1L
  )
  expect_true(
    is.na(boin_select_mtd(n_pts, n_tox, target = 0.30, extrasafe = TRUE)$mtd)
  )
  expect_equal(
    boin_select_mtd(n_pts, n_tox, target = 0.30, extrasafe = TRUE)$reason,
    "lowest_dose_eliminated"
  )
})

test_that("the MTD is chosen from the isotonic fit over admissible doses only", {
  # Six DLTs out of six patients eliminates dose 3 and every dose above it, so
  # the selection uses the fit restricted to doses 1 and 2.
  n_pts <- matrix(c(6, 6, 6), nrow = 1)
  n_tox <- matrix(c(0, 2, 6), nrow = 1)

  out <- boin_select_mtd(n_pts, n_tox, target = 0.30)
  est <- boin_isotonic(n_pts, n_tox,
                       admissible = matrix(c(TRUE, TRUE, FALSE), nrow = 1))

  expect_equal(out$mtd, 2L)
  expect_true(is.na(est[1, 3]))
  expect_equal(unname(which.min(abs(est[1, 1:2] - 0.30))), 2L)
})

test_that("boin_select_mtd validates its arguments", {
  expect_error(boin_select_mtd(c(3, 6), c(0, 1, 2), target = 0.3), "same dimensions")
  expect_error(boin_select_mtd(c(3, 6), c(0, 9), target = 0.3), "must not exceed")
  expect_error(boin_select_mtd(c(3, 6), c(0, 1), target = 0.3, min_mtd_sample = 0),
               "at least 1")
})
