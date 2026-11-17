# Expected values from the Python implementation (Validation/comb_boin12/
# test_expect_cb12.py of the development folder); the raw DLT estimates of
# example E2 are already monotone, so that the isotonic regression returns
# them unchanged.
e2 <- list(n = matrix(c(3, 6, 9, 3, 3, 0), 3),
           y = matrix(c(0, 1, 3, 0, 1, 0), 3),
           x = matrix(c(0, 2, 4, 1, 2, 0), 3))

test_that("the rectangle and the contour bounds give the expected OBDC", {
  # MTDC (3, 1): the rectangle keeps column 1 only, the contour also keeps
  # (1, 2) and (2, 2) of column 2, where (2, 2) has the highest utility.
  args <- list(n_pts = e2$n, n_tox = e2$y, n_resp = e2$x, target = 0.30,
               target_resp = 0.20, est_method = "observed")
  rect <- do.call(comb_boin12_select_obd, args)
  cont <- do.call(comb_boin12_select_obd,
                  c(args, list(obd_bound = "contour")))
  expect_identical(c(rect$obd_a, rect$obd_b), c(3L, 1L))
  expect_identical(c(cont$obd_a, cont$obd_b), c(2L, 2L))
  expect_identical(c(rect$mtd_a, rect$mtd_b), c(3L, 1L))
  expect_identical(rect$reason, "selected")
  raw <- (e2$y + 0.05) / (e2$n + 0.1)
  expect_equal(unname(rect$tox_est[1, , ]), raw, tolerance = 1e-12)
  util <- c(39.354839, 53.114754, 53.260073, 59.354839, 66.451613, NA)
  expect_equal(as.vector(rect$utility[1, , ]), util, tolerance = 1e-6)
  expect_true(all(is.na(rect$resp_est[1, , ][e2$n == 0])))
  expect_equal(rect$resp_est[1, , ][e2$n > 0], (e2$x / e2$n)[e2$n > 0])
})

test_that("the final admissibility uses the cutoff of extrasafe", {
  # Two DLTs in three patients at (1, 1): Pr(pT > 0.30) = 0.9163, below 0.95
  # but not below 0.85.
  n <- matrix(c(3, 3, 3, 0), 2)
  y <- matrix(c(2, 2, 3, 0), 2)
  x <- matrix(c(1, 2, 2, 0), 2)
  args <- list(n_pts = n, n_tox = y, n_resp = x, target = 0.30,
               target_resp = 0.20, est_method = "observed")
  plain <- do.call(comb_boin12_select_obd, args)
  safe <- do.call(comb_boin12_select_obd,
                  c(args, list(extrasafe = TRUE, offset = 0.10)))
  expect_identical(c(plain$obd_a, plain$obd_b), c(1L, 1L))
  expect_true(is.na(safe$obd_a))
  expect_identical(safe$reason, "no_admissible_dose")
})

test_that("the logistic estimates agree with a separate implementation", {
  # Python with the same Halton points (4,000, two refits); long Markov
  # chain Monte Carlo runs differ from these values by at most 0.0022.
  n <- matrix(c(3, 6, 3, 0, 6, 12, 3, 0, 3, 3, 0, 0, 0, 0, 0, 0), 4)
  x <- matrix(c(0, 1, 1, 0, 1, 5, 2, 0, 1, 2, 0, 0, 0, 0, 0, 0), 4)
  y <- matrix(0, 4, 4)
  res <- comb_boin12_select_obd(n, y, x, target = 0.35, target_resp = 0.2,
                                dose_a = c(80, 160, 240, 320),
                                dose_b = c(75, 150, 225, 300))
  expect_equal(res$resp_est[1, , ][n > 0],
               c(0.08726675, 0.19621046, 0.32772446, 0.19566430,
                 0.40641758, 0.57382698, 0.32726284, 0.57450728),
               tolerance = 1e-6)
  n <- matrix(c(3, 6, 9, 12, 6, 6), 2)
  x <- matrix(c(0, 1, 3, 5, 3, 4), 2)
  res <- comb_boin12_select_obd(n, matrix(0, 2, 3), x, target = 0.33,
                                target_resp = 0.2, dose_a = c(100, 200),
                                dose_b = c(50, 100, 200))
  expect_equal(as.vector(res$resp_est[1, , ]),
               c(0.16357160, 0.22485910, 0.29540055, 0.39334188,
                 0.52284631, 0.62699096), tolerance = 1e-6)
  # No random number is used.
  set.seed(1)
  before <- .Random.seed
  again <- comb_boin12_select_obd(n, matrix(0, 2, 3), x, target = 0.33,
                                  target_resp = 0.2, dose_a = c(100, 200),
                                  dose_b = c(50, 100, 200))
  expect_identical(.Random.seed, before)
  expect_identical(again$resp_est, res$resp_est)
})

test_that("several trials, early stops and a single row are handled", {
  n <- array(0L, c(3L, 3L, 2L))
  n[1, , ] <- e2$n
  n[2, , ] <- e2$n
  n[3, , ] <- e2$n
  y <- n
  y[1, , ] <- e2$y
  y[2, , ] <- e2$y
  y[3, , ] <- e2$y
  x <- n
  x[1, , ] <- e2$x
  x[2, , ] <- e2$x
  x[3, , ] <- e2$x
  res <- comb_boin12_select_obd(n, y, x, target = 0.30, target_resp = 0.20,
                                est_method = "observed",
                                early_stop = c(FALSE, TRUE, FALSE))
  expect_identical(res$obd_a, c(3L, NA, 3L))
  expect_identical(res$reason, c("selected", "early_stop", "selected"))
  expect_true(all(is.na(res$utility[2, , ])))
  # One row: the isotonic regression is that of a single drug.
  one <- comb_boin12_select_obd(matrix(c(3, 6, 3), 1), matrix(c(0, 1, 2), 1),
                                matrix(c(0, 3, 2), 1), target = 0.30,
                                target_resp = 0.20, est_method = "observed")
  expect_equal(as.vector(one$tox_est[1, , ]),
               as.vector(Iso::pava(c(0.05, 1.05, 2.05) / c(3.1, 6.1, 3.1),
                                   c(3.1, 6.1, 3.1))))
  expect_false(is.na(one$obd_b))
})

test_that("invalid arguments are rejected", {
  args <- list(n_pts = e2$n, n_tox = e2$y, n_resp = e2$x, target = 0.30,
               target_resp = 0.20)
  expect_error(do.call(comb_boin12_select_obd,
                       modifyList(args, list(n_tox = e2$n + 1))), "n_tox")
  expect_error(do.call(comb_boin12_select_obd,
                       modifyList(args, list(n_resp = e2$x[1:2, ]))),
               "dimensions")
  expect_error(do.call(comb_boin12_select_obd,
                       modifyList(args, list(dose_a = c(3, 2, 1)))),
               "dose_a")
  expect_error(do.call(comb_boin12_select_obd,
                       modifyList(args, list(logistic_scale = 1))),
               "logistic_scale")
  expect_error(do.call(comb_boin12_select_obd,
                       modifyList(args, list(n_draws = 10))), "n_draws")
  expect_error(do.call(comb_boin12_select_obd,
                       modifyList(args, list(early_stop = NA))),
               "early_stop")
  expect_error(do.call(comb_boin12_select_obd,
                       modifyList(args, list(obd_bound = "box"))),
               "should be one of")
})
