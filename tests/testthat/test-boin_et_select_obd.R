# Final data of single trials and the OBD under the four selection rules, from
# a separate Python implementation of the rules of boinet::obd.select().
obd_cases <- list(
  list(n = c(3, 6, 12, 3), y = c(0, 1, 3, 2), x = c(0, 2, 7, 2),
       obd = c(3, 3, 3, 3)),
  list(n = c(3, 3, 9, 15, 6, 0), y = c(0, 0, 1, 3, 3, 0),
       x = c(0, 1, 4, 9, 3, 0), obd = c(4, 4, 4, 4)),
  list(n = c(3, 9, 6), y = c(0, 2, 5), x = c(1, 5, 4), obd = c(2, 2, 2, 2)),
  list(n = c(3, 12, 3), y = c(0, 1, 0), x = c(0, 0, 2), obd = c(3, 3, 3, 3)),
  list(n = c(9, 0, 0), y = c(2, 0, 0), x = c(4, 0, 0), obd = c(1, 1, 1, 1)),
  list(n = c(6, 6, 6), y = c(0, 0, 0), x = c(3, 3, 3), obd = c(1, 1, 1, 1)),
  list(n = c(3, 0, 9, 3), y = c(0, 0, 2, 1), x = c(2, 0, 6, 1),
       obd = c(1, 1, 1, 1)),
  list(n = c(3, 6, 15, 9, 3), y = c(0, 0, 2, 4, 2), x = c(1, 3, 9, 6, 1),
       obd = c(4, 3, 2, 3)),
  list(n = c(3, 6, 0, 12, 6), y = c(0, 2, 0, 5, 2), x = c(1, 4, 0, 6, 5),
       obd = c(2, 5, 1, 5)),
  list(n = c(6, 0, 6, 9, 3), y = c(0, 0, 2, 4, 1), x = c(3, 0, 4, 7, 2),
       obd = c(3, 4, 1, 4))
)
obd_methods <- c("max_resp", "utility_weighted", "utility_truncated_linear",
                 "utility_scoring")

test_that("the OBD of single trials agrees with an independent implementation", {
  for (case in obd_cases) {
    for (m in seq_along(obd_methods)) {
      sel <- boin_et_select_obd(case$n, case$y, case$x, target = 0.30,
                                target_resp = 0.60, obd_method = obd_methods[m])
      expect_identical(sel$obd, as.integer(case$obd[m]))
      expect_identical(sel$reason, "selected")
    }
  }
})

test_that("a matrix of trials gives the same OBDs as the trials one by one", {
  n_doses <- 6
  pad <- function(v) c(v, rep(0, n_doses - length(v)))
  n <- t(sapply(obd_cases, function(z) pad(z$n)))
  y <- t(sapply(obd_cases, function(z) pad(z$y)))
  x <- t(sapply(obd_cases, function(z) pad(z$x)))
  for (m in seq_along(obd_methods)) {
    sel <- boin_et_select_obd(n, y, x, target = 0.30, target_resp = 0.60,
                              obd_method = obd_methods[m])
    expect_identical(sel$obd,
                     as.integer(vapply(obd_cases, function(z) z$obd[m], 1)))
  }
  expect_equal(dim(sel$probe), c(length(obd_cases), n_doses))
  expect_true(all(is.na(sel$probe[n == 0])))
  expect_equal(sel$probe[n > 0], x[n > 0] / n[n > 0])
})

test_that("early stops and inadmissible doses select no OBD", {
  none <- boin_et_select_obd(c(6, 0), c(5, 0), c(0, 0), target = 0.30,
                             target_resp = 0.60)
  expect_identical(none$obd, NA_integer_)
  expect_identical(none$reason, "no_admissible_dose")

  stopped <- boin_et_select_obd(c(3, 6, 12, 3), c(0, 1, 3, 2), c(0, 2, 7, 2),
                                target = 0.30, target_resp = 0.60,
                                early_stop = TRUE)
  expect_identical(stopped$obd, NA_integer_)
  expect_identical(stopped$reason, "early_stop")
})

test_that("the OBD agrees with boinet::obd.select() on simulated trials", {
  skip_if_not_installed("boinet")
  boinet_names <- c(max_resp = "max.effprob",
                    utility_weighted = "utility.weighted",
                    utility_truncated_linear = "utility.truncated.linear",
                    utility_scoring = "utility.scoring")
  # The end of boinet::boinet(), applied to final data.
  reference <- function(n, y, x, est, method) {
    phi <- 0.30
    delta <- 0.60
    eva <- which(n != 0)
    tterm <- stats::pbeta(phi, 1 + y[eva], 1 + n[eva] - y[eva])
    eterm <- 1 - stats::pbeta(delta * 0.6, 1 + x[eva], 1 + n[eva] - x[eva])
    ok <- (tterm >= (1 - 0.95)) & (eterm >= (1 - 0.99))
    if (length(eva) == 1L) return(if (ok) eva else NA_integer_)
    if (!any(ok)) return(NA_integer_)
    estpt <- Iso::pava(y[eva] / n[eva])
    estpe <- if (est == "multi_iso") {
      boinet::multi.iso(obs = x[eva], n = n[eva])
    } else {
      x[eva] / n[eva]
    }
    eva[boinet::obd.select(
      probt = estpt, probe = estpe, method = boinet_names[[method]],
      phi = phi, phi1 = phi * 0.1, phi2 = phi * 1.4, delta = delta,
      delta1 = delta * 0.6, tterm = tterm, eterm = eterm, stopT = 0.95,
      stopE = 0.99, w1 = 0.33, w2 = 1.09, plow.ast = phi * 0.1,
      pupp.ast = phi * 1.4, qlow.ast = delta * 0.6 / 2, qupp.ast = delta,
      psi00 = 40, psi11 = 60)]
  }
  scenarios <- list(
    list(p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
         p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70)),
    list(p_true = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
         p_resp = c(0.30, 0.70, 0.60, 0.50, 0.40, 0.30)),
    list(p_true = c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45),
         p_resp = c(0.65, 0.60, 0.55, 0.50, 0.45, 0.40))
  )
  checked <- 0L
  for (s in seq_along(scenarios)) {
    tr <- boin_et_simulate(target = 0.30, target_resp = 0.60,
                           p_true = scenarios[[s]]$p_true,
                           p_resp = scenarios[[s]]$p_resp, n_cohort = 12,
                           cohort_size = 3, n_trials = 100, seed = s)
    keep <- !tr$early_stop
    for (est in c("observed", "multi_iso")) {
      for (method in obd_methods) {
        sel <- boin_et_select_obd(tr$n_pts[keep, ], tr$n_tox[keep, ],
                                  tr$n_resp[keep, ], target = 0.30,
                                  target_resp = 0.60, est_method = est,
                                  obd_method = method)
        expected <- vapply(which(keep), function(i) {
          as.integer(reference(tr$n_pts[i, ], tr$n_tox[i, ], tr$n_resp[i, ],
                               est, method))
        }, integer(1))
        expect_identical(sel$obd, expected)
        checked <- checked + length(expected)
      }
    }
  }
  expect_gt(checked, 2000L)
})

test_that("the estimators agree with those of boinet", {
  skip_if_not_installed("boinet")
  obs <- c(1, 4, 6, 3, 1)
  n <- c(3, 6, 12, 9, 3)
  expect_identical(boin_et_multi_iso(obs, n), boinet::multi.iso(obs = obs, n = n))
  skip_if_not_installed("mfp")
  expect_identical(boin_et_fp_logit(obs, n, dose = 1:5),
                   boinet::fp.logit(obs = obs, n = n, dose = 1:5))
  expect_identical(boin_et_fp_logit(obs[-2], n[-2], dose = c(1, 3, 4, 5)),
                   boinet::fp.logit(obs = obs[-2], n = n[-2], dose = c(1, 3, 4, 5)))
})

test_that("the fractional polynomial estimate enters the selection", {
  skip_if_not_installed("mfp")
  n <- rbind(c(3, 6, 15, 9, 3), c(3, 6, 12, 3, 0))
  y <- rbind(c(0, 0, 2, 4, 2), c(0, 1, 3, 2, 0))
  x <- rbind(c(1, 3, 9, 6, 1), c(0, 2, 7, 2, 0))
  sel <- boin_et_select_obd(n, y, x, target = 0.30, target_resp = 0.60,
                            est_method = "fp_logistic")
  expect_identical(sel$reason, c("selected", "selected"))
  expect_equal(sel$probe[1, ], boin_et_fp_logit(x[1, ], n[1, ], 1:5),
               ignore_attr = TRUE)
  expect_equal(sel$probe[2, 1:4], boin_et_fp_logit(x[2, 1:4], n[2, 1:4], 1:4),
               ignore_attr = TRUE)
})

test_that("invalid input is rejected", {
  expect_error(boin_et_select_obd(c(3, 3), c(4, 0), c(0, 0), 0.3, 0.6),
               "must not exceed")
  expect_error(boin_et_select_obd(c(3, 3), c(0, 0), c(0, 0, 0), 0.3, 0.6),
               "same dimensions")
  expect_error(boin_et_select_obd(c(0, 0), c(0, 0), c(0, 0), 0.3, 0.6),
               "at least one patient")
  expect_error(boin_et_select_obd(c(3, 3), c(0, 0), c(0, 0), 0.3, 0.6,
                                  early_stop = c(TRUE, FALSE)), "early_stop")
  expect_error(boin_et_select_obd(c(3, 3), c(0, 0), c(0, 0), 0.3, 0.6,
                                  w1 = NA), "w1")
  expect_error(boin_et_select_obd(c(3, 3), c(0, 0), c(0, 0), 0.3, 0.6,
                                  obd_method = "other"))
})
