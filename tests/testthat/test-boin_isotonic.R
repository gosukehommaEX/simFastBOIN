test_that("boin_isotonic returns monotone estimates", {
  n_pts <- matrix(c(3, 6, 9, 12,
                    3, 6, 9, 12,
                    3, 3, 3, 3), nrow = 3, byrow = TRUE)
  n_tox <- matrix(c(0, 1, 3, 4,
                    0, 0, 2, 3,
                    0, 3, 0, 1), nrow = 3, byrow = TRUE)

  est <- boin_isotonic(n_pts, n_tox)

  expect_equal(dim(est), dim(n_pts))
  expect_true(all(est >= 0 & est <= 1))
  for (i in seq_len(nrow(est))) {
    expect_true(all(diff(est[i, ]) >= -1e-12))
  }
})

test_that("boin_isotonic reproduces known pooled estimates", {
  expect_equal(
    round(as.vector(boin_isotonic(c(3, 6, 9, 12), c(0, 1, 3, 4))), 6),
    c(0.016129, 0.172131, 0.334908, 0.334908)
  )
  expect_equal(
    round(as.vector(boin_isotonic(c(3, 6, 9, 12), c(0, 0, 2, 3))), 6),
    c(0.010008, 0.010008, 0.225275, 0.252066)
  )
  expect_equal(
    round(as.vector(boin_isotonic(c(3, 6, 9, 12), c(1, 2, 4, 6))), 6),
    c(0.337031, 0.337031, 0.445055, 0.5)
  )
  # Equal weights, so a reversal of 0 of 3 and 3 of 3 pools to exactly one half.
  expect_equal(
    as.vector(boin_isotonic(c(3, 3, 3), c(0, 3, 0))),
    c(0.05 / 3.1, 0.5, 0.5)
  )
})

test_that("untreated and inadmissible doses are returned as NA", {
  est <- boin_isotonic(c(3, 6, 0, 0), c(0, 1, 0, 0))
  expect_true(all(is.na(est[1, 3:4])))
  expect_false(any(is.na(est[1, 1:2])))

  admissible <- matrix(c(TRUE, TRUE, FALSE, FALSE), nrow = 1)
  restricted <- boin_isotonic(c(3, 6, 9, 12), c(0, 1, 3, 4), admissible = admissible)
  expect_true(all(is.na(restricted[1, 3:4])))
  expect_equal(restricted[1, 1:2], est[1, 1:2], ignore_attr = TRUE)
})

test_that("excluding a dose can change the estimates of the doses kept", {
  # Dose 3 is a reversal relative to dose 2, so dropping it removes the pooling.
  full <- boin_isotonic(c(6, 6, 6), c(1, 4, 0))
  dropped <- boin_isotonic(c(6, 6, 6), c(1, 4, 0),
                           admissible = matrix(c(TRUE, TRUE, FALSE), nrow = 1))

  expect_equal(full[1, 2], full[1, 3])
  expect_false(isTRUE(all.equal(full[1, 2], dropped[1, 2])))
})

test_that("boin_isotonic agrees with Iso::pava", {
  skip_if_not_installed("Iso")
  expect_true(have_package("Iso"))

  set.seed(42)
  n_doses <- 6
  for (i in 1:50) {
    n <- sample(0:12, n_doses, replace = TRUE)
    y <- vapply(n, function(k) if (k == 0) 0L else sample(0:k, 1), integer(1))
    treated <- n > 0
    if (!any(treated)) next

    phat <- (y[treated] + 0.05) / (n[treated] + 0.1)
    variance <- (y[treated] + 0.05) * (n[treated] - y[treated] + 0.05) /
      ((n[treated] + 0.1)^2 * (n[treated] + 0.1 + 1))
    reference <- Iso::pava(phat, w = 1 / variance)

    expect_equal(as.vector(boin_isotonic(n, y))[treated], reference,
                 tolerance = 1e-10)
  }
})

test_that("boin_isotonic validates its arguments", {
  expect_error(boin_isotonic(c(3, 6), c(0, 1, 2)), "same dimensions")
  expect_error(boin_isotonic(c(3, 6), c(0, 7)), "must not exceed")
  expect_error(boin_isotonic(c(3, -1), c(0, 0)), "non-negative whole numbers")
  expect_error(boin_isotonic(c(3, 6), c(0, 1), admissible = c(TRUE, TRUE, TRUE)),
               "same dimensions")
})
