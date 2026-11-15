test_that("the boundaries of Takeda et al. (2018) are reproduced", {
  bd <- boin_et_boundary(target = 0.30, target_resp = 0.60)

  expect_s3_class(bd, "boin_et_boundary")
  expect_equal(bd$lambda1, 0.14)
  expect_equal(bd$lambda2, 0.35)
  expect_equal(bd$eta1, 0.48)
  expect_equal(c(bd$p_saf, bd$p_tox, bd$p_saf_resp), c(0.03, 0.42, 0.36))
  expect_gt(bd$prob_incorrect, 0)
  expect_lt(bd$prob_incorrect, 1)
})

test_that("other settings agree with an independent implementation", {
  # Expected values from a separate Python implementation of the grid search.
  cases <- list(
    list(args = list(target = 0.25, target_resp = 0.40),
         expected = c(0.115, 0.29, 0.31)),
    list(args = list(target = 0.30, target_resp = 0.30),
         expected = c(0.13, 0.35, 0.23)),
    list(args = list(target = 0.20, target_resp = 0.50),
         expected = c(0.09, 0.23, 0.40)),
    list(args = list(target = 0.30, target_resp = 0.60, p_saf = 0.18),
         expected = c(0.29, 0.35, 0.47)),
    list(args = list(target = 0.33, target_resp = 0.45, p_tox = 0.5,
                     p_saf_resp = 0.2),
         expected = c(0.153, 0.41, 0.31))
  )
  for (case in cases) {
    bd <- do.call(boin_et_boundary, case$args)
    expect_equal(c(bd$lambda1, bd$lambda2, bd$eta1), case$expected)
  }
})

test_that("the boundaries agree with boinet::gridoptim()", {
  skip_if_not_installed("boinet")
  settings <- list(c(0.30, 0.60), c(0.25, 0.40), c(0.20, 0.35), c(0.35, 0.50))
  for (s in settings) {
    ref <- boinet::gridoptim(phi = s[1], phi1 = s[1] * 0.1, phi2 = s[1] * 1.4,
                             delta = s[2], delta1 = s[2] * 0.6)
    bd <- boin_et_boundary(target = s[1], target_resp = s[2])
    expect_identical(c(bd$lambda1, bd$lambda2, bd$eta1),
                     c(ref$lambda1, ref$lambda2, ref$eta1))
  }
})

test_that("the prior changes the boundaries", {
  flat <- boin_et_boundary(target = 0.30, target_resp = 0.60)
  skewed <- boin_et_boundary(target = 0.30, target_resp = 0.60,
                             prior = c(0.05, 0.05, 0.1, 0.1, 0.35, 0.35))
  expect_false(identical(c(flat$lambda1, flat$lambda2, flat$eta1),
                         c(skewed$lambda1, skewed$lambda2, skewed$eta1)))
})

test_that("invalid design parameters are rejected", {
  expect_error(boin_et_boundary(target = 0.30, target_resp = 0.60, p_saf = 0.3),
               "p_saf < target")
  expect_error(boin_et_boundary(target = 0.30, target_resp = 0.60, p_tox = 0.25),
               "p_saf < target")
  expect_error(boin_et_boundary(target = 0.30, target_resp = 0.60,
                                p_saf_resp = 0.7), "p_saf_resp")
  expect_error(boin_et_boundary(target = 1.2, target_resp = 0.60), "target")
  expect_error(boin_et_boundary(target = 0.30, target_resp = 0.60,
                                prior = rep(1 / 5, 5)), "prior")
})
