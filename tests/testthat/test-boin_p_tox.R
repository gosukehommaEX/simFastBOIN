test_that("boin_p_tox inverts boin_lambda", {
  for (target in c(0.15, 0.20, 0.25, 0.30, 0.40)) {
    for (p_tox in c(1.2, 1.4, 2.0) * target) {
      if (p_tox >= 1) next
      lambda_d <- boin_lambda(target, p_tox = p_tox)$lambda_d
      expect_equal(boin_p_tox(target, lambda_d), p_tox, tolerance = 1e-6)
    }
  }
})

test_that("the returned threshold produces the requested boundary", {
  for (target in c(0.20, 0.25, 0.30)) {
    for (lambda_d in c(target + 0.03, target + 0.06, target + 0.10)) {
      p_tox <- boin_p_tox(target, lambda_d)
      expect_equal(boin_lambda(target, p_tox = p_tox)$lambda_d, lambda_d,
                   tolerance = 1e-8)
      expect_gt(p_tox, target)
      expect_lt(p_tox, 1)
    }
  }
})

test_that("the de-escalation boundary lies between the target and the threshold", {
  p_tox <- boin_p_tox(target = 0.30, lambda_d = 0.33)

  expect_gt(p_tox, 0.30)
  expect_lt(p_tox, 0.42)
  expect_equal(round(p_tox, 4), 0.3610)
})

test_that("tightening the boundary lowers the de-escalation counts", {
  p_tox <- boin_p_tox(target = 0.30, lambda_d = 0.33)

  default <- boin_boundary(0.30, 18)$b_deesc[seq(3, 18, by = 3)]
  tighter <- boin_boundary(0.30, 18, p_tox = p_tox)$b_deesc[seq(3, 18, by = 3)]

  expect_equal(default, c(2L, 3L, 4L, 5L, 6L, 7L))
  expect_equal(tighter, c(1L, 2L, 3L, 4L, 5L, 6L))
})

test_that("a boundary too close to the target is flagged", {
  # The other functions reject 'p_tox' within ten percent of 'target', and a
  # boundary just above the target lands inside that band.
  expect_warning(p <- boin_p_tox(0.30, 0.31), "too close to 'target'")
  expect_lt(p - 0.30, 0.1 * 0.30)
  expect_warning(boin_p_tox(0.30, 0.31), "smallest usable boundary")

  expect_warning(boin_p_tox(0.25, 0.26), "too close to 'target'")

  # A boundary far enough above the target is returned without complaint.
  expect_silent(usable <- boin_p_tox(0.30, 0.33))
  expect_gte(usable - 0.30, 0.1 * 0.30)
  expect_silent(usable <- boin_p_tox(0.25, 0.28))
  expect_gte(usable - 0.25, 0.1 * 0.25)
})

test_that("the warning threshold matches where the other functions give way", {
  # Just below the smallest usable boundary the value warns and is rejected
  # downstream; just above it the whole chain goes through.
  smallest <- 0.3149

  expect_warning(too_tight <- boin_p_tox(0.30, smallest - 0.001))
  expect_error(boin_lambda(0.30, p_tox = too_tight), "clearly above")

  expect_silent(usable <- boin_p_tox(0.30, smallest + 0.001))
  expect_equal(boin_lambda(0.30, p_tox = usable)$lambda_d, smallest + 0.001,
               tolerance = 1e-8)
})

test_that("boin_p_tox validates its arguments", {
  expect_error(boin_p_tox(0.30, 0.30), "greater than 'target'")
  expect_error(boin_p_tox(0.30, 0.20), "greater than 'target'")
  expect_error(boin_p_tox(0.30, 1.2), "strictly between")
  expect_error(boin_p_tox(c(0.2, 0.3), 0.35), "single finite number")
})
