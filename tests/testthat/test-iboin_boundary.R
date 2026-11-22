zhou_skeleton <- function() c(0.10, 0.19, 0.30, 0.42, 0.54)

test_that("iboin_boundary reproduces Table 1 of Zhou et al. (2021)", {
  bd <- iboin_boundary(target = 0.30, skeleton = zhou_skeleton(),
                       prior_ess = 3, max_n = 30)
  at <- seq(3, 30, by = 3)

  # Escalate if the number of DLTs is at most (doses in rows)
  escalate <- rbind(
    c(1, 1, 2, 3, 4, 4, 5, 6, 6, 7),
    c(0, 1, 2, 3, 3, 4, 5, 5, 6, 7),
    c(0, 1, 2, 2, 3, 4, 4, 5, 6, 7),
    c(0, 1, 1, 2, 3, 3, 4, 5, 6, 6),
    c(0, 0, 1, 2, 2, 3, 4, 5, 5, 6)
  )
  # De-escalate if the number of DLTs is at least
  deescalate <- rbind(
    c(2, 3, 4, 5, 7, 8, 9, 10, 11, 12),
    c(2, 3, 4, 5, 6, 7, 8, 9, 11, 12),
    c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
    c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11),
    c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11)
  )

  expect_s3_class(bd, "iboin_boundary")
  expect_equal(unname(bd$b_esc[, at]), escalate)
  expect_equal(unname(bd$b_deesc[, at]), deescalate)
  expect_false(any(bd$searched))
  expect_equal(dim(bd$b_esc), c(5L, 30L))
})

test_that("the prior probabilities follow equation (5)", {
  bd <- iboin_boundary(0.30, zhou_skeleton(), 3, 9)

  # Computed apart from the package (Python, scipy)
  expected <- rbind(
    c(0.3258969787, 0.4445807483, 0.2295222730),
    c(0.3321544593, 0.3920553963, 0.2757901445),
    c(0.3352362784, 0.3316281176, 0.3331356041),
    c(0.3331356041, 0.2706155062, 0.3962488897),
    c(0.3256381943, 0.2149274412, 0.4594343644)
  )
  expect_lt(max(abs(unname(bd$prior) - expected)), 1e-9)
  expect_equal(unname(rowSums(bd$prior)), rep(1, 5), tolerance = 1e-12)
  expect_equal(colnames(bd$prior), c("target", "p_saf", "p_tox"))

  # lambda_e and lambda_d of the lowest dose after three patients
  expect_lt(abs(bd$lambda_e[1, 3] - 0.3912128317), 1e-9)
  expect_lt(abs(bd$lambda_d[1, 3] - 0.5813125288), 1e-9)
})

test_that("a prior effective sample size of 0 gives the BOIN boundaries", {
  for (target in c(0.10, 0.25, 0.30)) {
    bd <- iboin_boundary(target, zhou_skeleton(), 0, 30, extrasafe = TRUE)
    ref <- boin_boundary(target, 30, extrasafe = TRUE)
    for (j in 1:5) {
      expect_identical(unname(bd$b_esc[j, ]), ref$b_esc)
      expect_identical(unname(bd$b_deesc[j, ]), ref$b_deesc)
    }
    expect_identical(bd$b_elim, ref$b_elim)
    expect_identical(bd$b_stop, ref$b_stop)
    expect_equal(unname(bd$prior), matrix(1 / 3, 5, 3), tolerance = 1e-15)
  }
})

test_that("the prior effective sample size can differ between doses", {
  mixed <- iboin_boundary(0.30, zhou_skeleton(), c(3, 0, 6, 1, 0), 12)
  expect_identical(mixed$prior_ess, c(3L, 0L, 6L, 1L, 0L))
  for (j in 1:5) {
    single <- iboin_boundary(0.30, zhou_skeleton(), c(3, 0, 6, 1, 0)[j], 12)
    expect_identical(mixed$b_esc[j, ], single$b_esc[j, ])
    expect_identical(mixed$b_deesc[j, ], single$b_deesc[j, ])
  }
})

test_that("the robust prior drops the information above a high prior MTD", {
  # Prior MTD 3 of 5 doses: the doses above it lose their prior information
  bd <- iboin_boundary(0.30, zhou_skeleton(), 3, 9, robust = TRUE)
  expect_identical(bd$prior_mtd, 3L)
  expect_identical(bd$prior_ess, c(3L, 3L, 3L, 0L, 0L))
  expect_identical(bd$prior_ess_input, rep(3L, 5))
  ref <- boin_boundary(0.30, 9)
  plain <- iboin_boundary(0.30, zhou_skeleton(), 3, 9)
  for (j in 4:5) {
    expect_identical(unname(bd$b_esc[j, ]), ref$b_esc)
    expect_identical(unname(bd$b_deesc[j, ]), ref$b_deesc)
  }
  expect_identical(bd$b_esc[1:3, ], plain$b_esc[1:3, ])

  # Prior MTD 1 of 5 doses: unchanged
  low <- iboin_boundary(0.30, c(0.30, 0.42, 0.54, 0.64, 0.73), 3, 9,
                        robust = TRUE)
  expect_identical(low$prior_mtd, 1L)
  expect_identical(low$prior_ess, rep(3L, 5))

  # Prior MTD at the highest dose: nothing above it
  high <- iboin_boundary(0.30, c(0.01, 0.04, 0.10, 0.19, 0.30), 3, 9,
                         robust = TRUE)
  expect_identical(high$prior_mtd, 5L)
  expect_identical(high$prior_ess, rep(3L, 5))

  # Prior MTD 2 of 4 doses counts as the upper half
  four <- iboin_boundary(0.30, c(0.2, 0.3, 0.4, 0.5), 3, 9, robust = TRUE)
  expect_identical(four$prior_ess, c(3L, 3L, 0L, 0L))
})

test_that("boundaries are searched where lambda_e exceeds lambda_d", {
  bd <- iboin_boundary(0.10, c(0.04, 0.13, 0.16, 0.25, 0.40), 8, 6)

  # Computed apart from the package (Python, scipy)
  searched <- rbind(
    c(0, 0, 0, 0, 0, 0),
    c(1, 0, 0, 0, 0, 0),
    c(1, 1, 0, 0, 0, 0),
    c(1, 1, 1, 1, 1, 1),
    c(1, 1, 1, 1, 1, 1)
  ) == 1
  escalate <- rbind(
    c(0, 0, 0, 0, 0, 0),
    c(-1, -1, 0, 0, 0, 0),
    c(-1, -1, -1, -1, 0, 0),
    c(-1, -1, -1, -1, -1, -1),
    c(-1, -1, -1, -1, -1, -1)
  )
  deescalate <- rbind(
    c(1, 1, 1, 2, 2, 2),
    c(0, 0, 1, 1, 1, 1),
    c(0, 0, 0, 0, 1, 1),
    c(0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0)
  )
  expect_identical(unname(bd$searched), searched)
  expect_equal(unname(bd$b_esc), escalate)
  expect_equal(unname(bd$b_deesc), deescalate)
  expect_identical(bd$b_elim, c(NA, NA, 2L, 2L, 2L, 2L))
  expect_true(all((bd$lambda_e > bd$lambda_d) == bd$searched))
})

test_that("iboin_search minimizes the probability of an incorrect decision", {
  # Expected pairs computed apart from the package (Python, scipy)
  expect_equal(iboin_search(3, c(0.02, 0.49, 0.49), 0.3, 0.18, 0.42),
               c(0, 1))
  expect_equal(iboin_search(6, c(0.02, 0.49, 0.49), 0.3, 0.18, 0.42),
               c(1, 2))
  expect_equal(iboin_search(9, c(0.05, 0.60, 0.35), 0.3, 0.18, 0.42),
               c(3, 4))
  expect_equal(iboin_search(12, c(0.01, 0.30, 0.69), 0.3, 0.18, 0.42),
               c(2, 3))

  # With equal prior probabilities the search finds the BOIN boundaries
  ref <- boin_boundary(0.30, 12)
  for (n in c(3, 6, 9, 12)) {
    found <- iboin_search(n, rep(1 / 3, 3), 0.3, 0.18, 0.42)
    expect_equal(found, c(ref$b_esc[n], floor(ref$lambda_d * n) + 1))
  }
})

test_that("the boundaries are always ordered and capped by the elimination", {
  settings <- list(
    list(0.30, zhou_skeleton(), 3),
    list(0.30, c(0.01, 0.04, 0.10, 0.19, 0.30), 20),
    list(0.30, c(0.30, 0.42, 0.54, 0.64, 0.73), 20),
    list(0.10, c(0.04, 0.13, 0.16, 0.25, 0.40), 12),
    list(0.20, c(0.60, 0.70, 0.80, 0.90, 0.95), 9)
  )
  for (s in settings) {
    bd <- iboin_boundary(s[[1]], s[[2]], s[[3]], 30)
    expect_true(all(bd$b_esc < bd$b_deesc))
    expect_true(all(bd$b_esc >= -1L))
    elim <- matrix(bd$b_elim, 5, 30, byrow = TRUE)
    has <- !is.na(elim)
    expect_true(all(bd$b_deesc[has] <= elim[has]))
    n <- matrix(1:30, 5, 30, byrow = TRUE)
    expect_true(all(bd$b_deesc <= n + 1L))
  }
})

test_that("iboin_boundary rejects invalid prior information", {
  expect_error(iboin_boundary(0.30, c(0.1, 0.3, 1), 3, 9), "skeleton")
  expect_error(iboin_boundary(0.30, c(0.1, 0.3, 0), 3, 9), "skeleton")
  expect_error(iboin_boundary(0.30, "a", 3, 9), "skeleton")
  expect_error(iboin_boundary(0.30, zhou_skeleton(), 2.5, 9), "prior_ess")
  expect_error(iboin_boundary(0.30, zhou_skeleton(), -1, 9), "prior_ess")
  expect_error(iboin_boundary(0.30, zhou_skeleton(), c(3, 3), 9),
               "prior_ess")
  expect_error(iboin_boundary(0.30, zhou_skeleton(), 3, 9, robust = NA),
               "robust")
  expect_error(iboin_boundary(0.30, zhou_skeleton(), 3, 0), "max_n")
})
