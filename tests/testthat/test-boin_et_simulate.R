test_that("boin_et_simulate returns consistent trial data", {
  tr <- boin_et_simulate(
    target = 0.30, target_resp = 0.60,
    p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
    p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
    n_cohort = 12, cohort_size = 3, n_trials = 500, seed = 1
  )

  expect_s3_class(tr, "boin_et_trials")
  expect_equal(dim(tr$n_pts), c(500L, 6L))
  expect_identical(colnames(tr$n_resp), paste0("DL", 1:6))
  expect_true(all(tr$n_tox <= tr$n_pts))
  expect_true(all(tr$n_resp <= tr$n_pts))
  expect_identical(rowSums(tr$n_pts), 3 * tr$cohorts_used)
  expect_true(all(tr$stop_reason %in% c("no_admissible_dose",
                                        "no_lower_admissible_dose",
                                        "n_earlystop", "max_cohorts",
                                        "top_dose_excluded")))
  expect_identical(tr$design, "BOIN-ET")
  expect_false(tr$settings$boinet_compatible)
  expect_identical(tr$early_stop,
                   tr$stop_reason %in% c("no_admissible_dose",
                                         "no_lower_admissible_dose"))
  expect_equal(c(tr$boundary$lambda1, tr$boundary$lambda2, tr$boundary$eta1),
               c(0.14, 0.35, 0.48))
  expect_identical(tr$settings$n_earlystop, 36L)
})

test_that("the seed reproduces the trials and the session state is kept", {
  args <- list(target = 0.30, target_resp = 0.60,
               p_true = c(0.10, 0.20, 0.30), p_resp = c(0.2, 0.4, 0.6),
               n_cohort = 6, cohort_size = 3, n_trials = 200)
  set.seed(99)
  before <- .Random.seed
  a <- do.call(boin_et_simulate, c(args, seed = 7))
  expect_identical(.Random.seed, before)
  b <- do.call(boin_et_simulate, c(args, seed = 7))
  expect_identical(a$n_pts, b$n_pts)
  expect_identical(a$n_resp, b$n_resp)
  c2 <- do.call(boin_et_simulate, c(args, seed = 8))
  expect_false(identical(a$n_pts, c2$n_pts))
})

test_that("a design without DLTs or responses escalates to the top dose", {
  # No DLT and no response: escalate every cohort, then stay at the top dose
  # until 12 patients without a response exclude it, which ends the trial.
  # Under the rules of boinet the trial stays there to the end.
  args <- list(target = 0.30, target_resp = 0.60, p_true = c(0, 0, 0),
               p_resp = c(0, 0, 0), n_cohort = 12, cohort_size = 3,
               n_trials = 50, seed = 3)
  tr <- do.call(boin_et_simulate, args)
  expect_true(all(tr$n_pts[, 1] == 3L & tr$n_pts[, 2] == 3L &
                    tr$n_pts[, 3] == 12L))
  expect_true(all(tr$stop_reason == "top_dose_excluded"))
  expect_false(any(tr$early_stop))
  tr <- do.call(boin_et_simulate, c(args, boinet_compatible = TRUE))
  expect_true(all(tr$n_pts[, 1] == 3L & tr$n_pts[, 2] == 3L &
                    tr$n_pts[, 3] == 30L))
  expect_true(all(tr$stop_reason == "max_cohorts"))
})

test_that("an escalation does not skip a dose excluded by its DLT data", {
  # Dose 2 has 3 DLTs in 3 patients, which excludes it, and the trial returns
  # to dose 1, where escalation is decided again. The trial stays at dose 1;
  # under the rules of boinet it skips to the untried dose 3 (independent
  # Python implementation of both rules).
  args <- list(target = 0.30, target_resp = 0.60, p_true = c(0, 1, 1),
               p_resp = c(0, 0, 0), n_cohort = 4, cohort_size = 3,
               n_trials = 20, seed = 4)
  tr <- do.call(boin_et_simulate, args)
  expect_true(all(tr$n_pts[, 1] == 9L & tr$n_pts[, 2] == 3L &
                    tr$n_pts[, 3] == 0L))
  tr <- do.call(boin_et_simulate, c(args, boinet_compatible = TRUE))
  expect_true(all(tr$n_pts[, 1] == 6L & tr$n_pts[, 2] == 3L &
                    tr$n_pts[, 3] == 3L))
  expect_true(all(tr$n_tox[, 3] == 3L))
})

test_that("a dose that always responds keeps the trial there", {
  tr <- boin_et_simulate(
    target = 0.30, target_resp = 0.60, p_true = c(0, 0.1, 0.2),
    p_resp = c(1, 1, 1), n_cohort = 12, cohort_size = 3, n_trials = 50,
    n_earlystop = 12, seed = 3
  )
  expect_true(all(tr$n_pts[, 1] == 12L))
  expect_true(all(tr$stop_reason == "n_earlystop"))
})

test_that("a toxic lowest dose stops the trial early", {
  # Excluding dose 1 by its DLT data excludes every dose. Under the rules of
  # boinet the untried doses stay admissible, and the trial stops because no
  # dose at or below the intended one is admissible.
  args <- list(target = 0.30, target_resp = 0.60, p_true = c(1, 1, 1),
               p_resp = c(0.5, 0.5, 0.5), n_cohort = 6, cohort_size = 3,
               n_trials = 20, seed = 3)
  tr <- do.call(boin_et_simulate, args)
  expect_true(all(tr$stop_reason == "no_admissible_dose"))
  expect_true(all(tr$early_stop))
  expect_true(all(tr$n_pts[, 1] == 3L))
  tr <- do.call(boin_et_simulate, c(args, boinet_compatible = TRUE))
  expect_true(all(tr$stop_reason == "no_lower_admissible_dose"))
  expect_true(all(tr$early_stop))
})

test_that("the copula and the binary link give the intended correlation", {
  # One dose and no exclusion, so that every trial treats 36 patients there
  # and the covariance of the DLT and response counts is 36 times that of a
  # single patient. Expected values from the bivariate normal distribution;
  # over 4000 trials the standard deviation of the estimate is about 0.1.
  base <- list(target = 0.30, target_resp = 0.60, p_true = 0.2, p_resp = 0.5,
               n_cohort = 12, cohort_size = 3, n_trials = 4000,
               cutoff_eli = 0.999999, cutoff_eli_resp = 0.999999)
  tr <- do.call(boin_et_simulate, c(base, resp_cor = 0.5, seed = 11))
  expect_true(all(tr$n_pts == 36L))
  expect_lt(abs(stats::cov(tr$n_tox[, 1], tr$n_resp[, 1]) - 2.0313), 0.5)
  expect_lt(abs(mean(tr$n_tox) / 36 - 0.2), 0.01)
  expect_lt(abs(mean(tr$n_resp) / 36 - 0.5), 0.01)

  tr <- do.call(boin_et_simulate, c(base, resp_cor = -0.4, seed = 12))
  expect_lt(abs(stats::cov(tr$n_tox[, 1], tr$n_resp[, 1]) + 1.6198), 0.5)

  # Binary link: covariance 36 * 0.5 * sqrt(0.2 * 0.8 * 0.5 * 0.5) = 3.6.
  tr <- do.call(boin_et_simulate, c(base, resp_cor = 0.5,
                                    resp_cor_type = "binary", seed = 13))
  expect_lt(abs(stats::cov(tr$n_tox[, 1], tr$n_resp[, 1]) - 3.6), 0.5)
  expect_lt(abs(mean(tr$n_resp) / 36 - 0.5), 0.01)
})

test_that("the conditional response probabilities of the binary link", {
  link <- boin_et_resp_link(p_true = c(0.2, 0.4), p_resp = c(0.5, 0.3),
                            resp_cor = 0.2)
  s <- sqrt(c(0.5 * 0.5 * 0.2 * 0.8, 0.3 * 0.7 * 0.4 * 0.6))
  p11 <- c(0.5 * 0.2, 0.3 * 0.4) + 0.2 * s
  expect_equal(link$q1, p11 / c(0.2, 0.4))
  expect_equal(link$q0, (c(0.5, 0.3) - p11) / c(0.8, 0.6))
  expect_error(boin_et_resp_link(0.05, 0.95, 0.9), "not attainable")
})

test_that("cohort sizes, the start dose and the stopping size are honoured", {
  tr <- boin_et_simulate(
    target = 0.30, target_resp = 0.60, p_true = c(0.05, 0.1, 0.2, 0.3),
    p_resp = c(0.1, 0.2, 0.4, 0.6), n_cohort = 5, cohort_size = c(1, 2, 3),
    start_dose = 2, n_earlystop = 6, n_trials = 300, seed = 5
  )
  expect_identical(tr$settings$cohort_size, c(1L, 2L, 3L, 3L, 3L))
  expect_true(all(tr$n_pts[, 2] >= 1L))
  sizes <- cumsum(c(1L, 2L, 3L, 3L, 3L))
  expect_true(all(rowSums(tr$n_pts) == sizes[tr$cohorts_used]))
  stopped <- tr$stop_reason == "n_earlystop"
  expect_true(any(stopped))
  expect_true(all(apply(tr$n_pts[stopped, , drop = FALSE], 1, max) >= 6L))
})

test_that("invalid input is rejected", {
  args <- list(target = 0.30, target_resp = 0.60, p_true = c(0.1, 0.2),
               p_resp = c(0.2, 0.4), n_cohort = 4, cohort_size = 3,
               n_trials = 10)
  expect_error(do.call(boin_et_simulate, modifyList(args, list(p_resp = 0.2))),
               "p_resp")
  expect_error(do.call(boin_et_simulate, modifyList(args, list(start_dose = 3))),
               "start_dose")
  expect_error(do.call(boin_et_simulate, modifyList(args, list(resp_cor = 1))),
               "resp_cor")
  expect_error(do.call(boin_et_simulate,
                       modifyList(args, list(resp_cor_type = "other"))))
  expect_error(do.call(boin_et_simulate,
                       modifyList(args, list(boinet_compatible = NA))),
               "boinet_compatible")
})
