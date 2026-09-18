test_that("boin_lambda reproduces the published boundaries for a 30% target", {
  # Liu and Yuan (2015) report lambda_e = 0.236 and lambda_d = 0.359 for a
  # target DLT rate of 0.30 with the default p_saf and p_tox.
  lambda <- boin_lambda(target = 0.30)

  expect_equal(round(lambda$lambda_e, 3), 0.236)
  expect_equal(round(lambda$lambda_d, 3), 0.359)
})

test_that("boin_lambda is stable across targets", {
  # Regression values for the default thresholds. They are cross-checked against
  # BOIN::get.boundary() in test-boin_reference.R when that package is installed.
  expect_equal(round(boin_lambda(target = 0.20)$lambda_e, 3), 0.157)
  expect_equal(round(boin_lambda(target = 0.20)$lambda_d, 3), 0.238)
  expect_equal(round(boin_lambda(target = 0.25)$lambda_e, 3), 0.197)
  expect_equal(round(boin_lambda(target = 0.25)$lambda_d, 3), 0.298)
  expect_equal(round(boin_lambda(target = 0.40)$lambda_e, 3), 0.316)
  expect_equal(round(boin_lambda(target = 0.40)$lambda_d, 3), 0.480)
  expect_equal(round(boin_lambda(target = 0.50)$lambda_e, 3), 0.397)
  expect_equal(round(boin_lambda(target = 0.50)$lambda_d, 3), 0.603)
})

test_that("boin_lambda orders the boundaries around the target", {
  for (target in c(0.10, 0.20, 0.25, 0.30, 0.40, 0.50)) {
    lambda <- boin_lambda(target)
    expect_lt(lambda$lambda_e, target)
    expect_gt(lambda$lambda_d, target)
  }
})

test_that("boin_lambda honors custom thresholds", {
  wide <- boin_lambda(target = 0.30, p_saf = 0.10, p_tox = 0.60)
  narrow <- boin_lambda(target = 0.30, p_saf = 0.20, p_tox = 0.40)

  expect_lt(wide$lambda_e, narrow$lambda_e)
  expect_gt(wide$lambda_d, narrow$lambda_d)
})

test_that("boin_lambda rejects unusable thresholds", {
  expect_error(boin_lambda(target = 0.03), "too low")
  expect_error(boin_lambda(target = 0.70), "too high")
  expect_error(boin_lambda(target = 0.30, p_saf = 0.29), "clearly below")
  expect_error(boin_lambda(target = 0.30, p_tox = 0.31), "clearly above")
  expect_error(boin_lambda(target = c(0.2, 0.3)), "single finite number")
})
