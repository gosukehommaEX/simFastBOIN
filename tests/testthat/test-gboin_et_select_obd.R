# Final data of eight trials: patients, and totals of the toxicity and
# efficacy scores on the raw scale (largest score 1.5).
obd_cases <- function() {
  list(
    n = rbind(c(3, 6, 12, 9, 0, 0), c(3, 6, 12, 9, 0, 0),
              c(6, 6, 6, 6, 6, 0), c(6, 3, 6, 9, 0, 0),
              c(6, 12, 9, 6, 0, 0), c(6, 0, 0, 0, 0, 0),
              c(3, 0, 0, 0, 0, 0), c(3, 3, 3, 3, 6, 12)),
    ets = rbind(c(0.5, 1.5, 5.5, 6.0, 0, 0), c(0.0, 1.0, 5.0, 5.5, 0, 0),
                c(0.5, 1.0, 1.5, 2.0, 2.5, 0), c(1.0, 4.5, 1.5, 3.0, 0, 0),
                c(1.0, 2.0, 3.5, 3.0, 0, 0), c(7.5, 0, 0, 0, 0, 0),
                c(0.5, 0, 0, 0, 0, 0), c(0, 0, 0.5, 0.5, 1.5, 4.0)),
    ees = rbind(c(1.0, 4.5, 12.0, 10.5, 0, 0), c(4.0, 3.0, 8.0, 6.5, 0, 0),
                c(1.0, 7.0, 8.0, 1.5, 2.5, 0), c(3.0, 2.0, 6.0, 9.5, 0, 0),
                c(2.0, 0.0, 9.0, 3.0, 0, 0), c(4.0, 0, 0, 0, 0, 0),
                c(2.0, 0, 0, 0, 0, 0), c(0, 1.0, 1.5, 3.0, 7.0, 12.5))
  )
}

select_cases <- function(...) {
  cases <- obd_cases()
  gboin_et_select_obd(
    n_pts = cases$n, n_tox = cases$ets / 1.5, n_resp = cases$ees / 1.5,
    target = 0.47 / 1.5, target_resp = 0.95 / 1.5, ...
  )
}

test_that("gboin_et_select_obd selects the doses of the rules", {
  # Expected doses from an implementation of the rules in Python, which
  # agrees with the C++ selection trial by trial on simulated trials
  res <- select_cases()
  expect_identical(res$obd, c(3L, 1L, 3L, 1L, 3L, NA, 1L, 5L))
  expect_identical(res$reason,
                   c(rep("selected", 5), "no_admissible_dose", "selected",
                     "selected"))
  # The estimates are the observed mean scores
  cases <- obd_cases()
  expect_equal(unname(res$probe[1, 1:4]), cases$ees[1, 1:4] / 1.5 /
                 cases$n[1, 1:4], tolerance = 1e-8)
  expect_true(all(is.na(res$probe[1, 5:6])))
  expect_identical(colnames(res$probe), paste0("DL", 1:6))
})

test_that("the linear fractional polynomial changes the selection", {
  # Expected doses and estimates as in the first test (numpy least squares)
  res <- select_cases(est_method = "fp_linear")
  expect_identical(res$obd, c(3L, 1L, 2L, 1L, 4L, NA, 1L, 5L))
  expect_equal(unname(res$probe[3, 1:5]),
               c(0.1052916799, 0.8484200670, 0.6709381608, 0.4171257717,
                 0.1804465428), tolerance = 1e-8)
  expect_equal(unname(res$probe[5, 1:4]),
               c(0.2122970762, 0.0764154663, 0.4012818607, 0.5885047557),
               tolerance = 1e-8)
  expect_equal(unname(res$probe[8, ]),
               c(0.0230716325, 0.1637772834, 0.3929729321, 0.6355618188,
                 0.7813612888, 0.6943623278), tolerance = 1e-8)
})

test_that("three doses or fewer keep the observed means", {
  n <- c(6, 9, 12, 0)
  ees <- c(2, 6.5, 5)
  obs <- gboin_et_select_obd(n, c(0.5, 1, 2, 0), c(ees, 0) / 1.5,
                             0.47 / 1.5, 0.95 / 1.5)
  fp <- gboin_et_select_obd(n, c(0.5, 1, 2, 0), c(ees, 0) / 1.5,
                            0.47 / 1.5, 0.95 / 1.5, est_method = "fp_linear")
  expect_identical(fp$probe, obs$probe)
  expect_identical(fp$obd, obs$obd)
})

test_that("min_obd_sample keeps doses with few patients from the selection", {
  # Expected doses as in the first test
  res <- select_cases(min_obd_sample = 6)
  expect_identical(res$obd, c(3L, 3L, 3L, 1L, 3L, NA, NA, 5L))
  expect_identical(res$reason[7], "no_admissible_dose")
  res <- select_cases(min_obd_sample = 6, est_method = "fp_linear")
  expect_identical(res$obd, c(3L, 3L, 2L, 1L, 4L, NA, NA, 5L))
})

test_that("boinet_compatible lets a dose above an excluded one be selected", {
  # Trial 4: dose 2 is excluded by its toxicity score. Expected doses as in
  # the first test.
  expect_identical(select_cases()$obd[4], 1L)
  expect_identical(select_cases(boinet_compatible = TRUE)$obd,
                   c(3L, 1L, 3L, 4L, 3L, NA, 1L, 5L))
})

test_that("the utilities select their doses", {
  # Expected doses as in the first test
  expect_identical(select_cases(obd_method = "utility_weighted")$obd,
                   c(3L, 1L, 3L, 1L, 3L, NA, 1L, 5L))
  expect_identical(select_cases(obd_method = "utility_truncated_linear")$obd,
                   c(2L, 1L, 2L, 1L, 3L, NA, 1L, 4L))
  expect_identical(select_cases(obd_method = "utility_scoring")$obd,
                   c(4L, 1L, 3L, 1L, 3L, NA, 1L, 5L))
  expect_identical(
    select_cases(obd_method = "utility_scoring", min_obd_sample = 6)$obd,
    c(4L, 2L, 3L, 1L, 3L, NA, NA, 5L)
  )
})

test_that("doses with the same utility as fractions tie", {
  # Doses 1 and 3 have the scoring utility 200 / 3: dose 1 from mean scores
  # of 0 and 4 / 9, dose 3 from 1 / 18 (after pooling) and 13 / 27. In
  # floating point the two utilities differ in the last digits, in a
  # direction that depends on the platform; rounded, they tie and the lowest
  # dose is selected. Expected doses as in the first test.
  n <- c(3, 6, 18, 3)
  tox <- c(0, 2 / 3, 0, 1 / 3)
  resp <- c(4 / 3, 7 / 3, 26 / 3, 1)
  sel <- function(method) {
    gboin_et_select_obd(n, tox, resp, 0.47 / 1.5, 0.95 / 1.5,
                        obd_method = method)$obd
  }
  expect_identical(sel("utility_scoring"), 1L)
  expect_identical(sel("max_resp"), 3L)
  expect_identical(sel("utility_weighted"), 3L)
  expect_identical(sel("utility_truncated_linear"), 3L)
})

test_that("trials stopped early select no dose", {
  res <- select_cases(early_stop = c(TRUE, rep(FALSE, 7)))
  expect_identical(res$obd[1], NA_integer_)
  expect_identical(res$reason[1], "early_stop")
  expect_identical(res$obd[-1], select_cases()$obd[-1])
})

test_that("whole numbers of DLTs and responses give the BOIN-ET selection", {
  set.seed(5)
  n <- matrix(sample(c(0, 3, 6, 9, 12), 200 * 5, replace = TRUE), ncol = 5)
  n[, 1] <- pmax(n[, 1], 3)
  y <- matrix(rbinom(length(n), n, 0.3), ncol = 5)
  x <- matrix(rbinom(length(n), n, 0.5), ncol = 5)
  # Under every rule: both functions compare the response estimates and the
  # utilities after rounding to ten decimals.
  for (compat in c(FALSE, TRUE)) {
    for (method in c("max_resp", "utility_weighted",
                     "utility_truncated_linear", "utility_scoring")) {
      ref <- boin_et_select_obd(n, y, x, 0.30, 0.60, obd_method = method,
                                boinet_compatible = compat)
      new <- gboin_et_select_obd(n, y, x, 0.30, 0.60, obd_method = method,
                                 boinet_compatible = compat)
      expect_identical(new$obd, ref$obd)
      expect_identical(new$reason, ref$reason)
      expect_equal(new$probe, ref$probe)
    }
  }
})

test_that("a vector is one trial", {
  res <- gboin_et_select_obd(c(3, 6, 12, 9), c(0.5, 1.5, 5.5, 6.0) / 1.5,
                             c(1.0, 4.5, 12.0, 10.5) / 1.5, 0.47 / 1.5,
                             0.95 / 1.5)
  expect_identical(res$obd, 3L)
  expect_identical(dim(res$probe), c(1L, 4L))
})

test_that("gboin_et_select_obd rejects invalid input", {
  n <- c(3, 6, 12, 9)
  expect_error(gboin_et_select_obd(n, c(1, 2, 3), c(1, 2, 3, 4), 0.3, 0.6),
               "same dimensions")
  expect_error(gboin_et_select_obd(n, c(4, 2, 3, 1), c(1, 2, 3, 4), 0.3, 0.6),
               "must not exceed")
  expect_error(gboin_et_select_obd(n, c(-1, 2, 3, 1), c(1, 2, 3, 4), 0.3,
                                   0.6), "n_tox")
  expect_error(gboin_et_select_obd(n, c(1, 2, 3, 1), c(1, NA, 3, 4), 0.3,
                                   0.6), "n_resp")
  expect_error(gboin_et_select_obd(c(0, 0, 0, 0), c(0, 0, 0, 0),
                                   c(0, 0, 0, 0), 0.3, 0.6), "at least one")
  expect_error(gboin_et_select_obd(n, c(1, 2, 3, 1), c(1, 2, 3, 4), 0.3, 0.6,
                                   min_obd_sample = 0), "min_obd_sample")
  expect_error(gboin_et_select_obd(n, c(1, 2, 3, 1), c(1, 2, 3, 4), 0.3, 0.6,
                                   est_method = "multi_iso"),
               "should be one of")
  expect_error(gboin_et_select_obd(n, c(1, 2, 3, 1), c(1, 2, 3, 4), 0.3, 0.6,
                                   early_stop = NA), "early_stop")
})
