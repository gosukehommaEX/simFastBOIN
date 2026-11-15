sc1 <- list(p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
            p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70))

test_that("tite_boin_et_simulate returns consistent trial data", {
  tr <- tite_boin_et_simulate(
    target = 0.30, target_resp = 0.60, p_true = sc1$p_true,
    p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3, window = 30,
    resp_window = 45, accrual_rate = 1 / 10, n_trials = 500, seed = 1
  )
  expect_s3_class(tr, "boin_et_trials")
  expect_identical(tr$design, "TITE-BOIN-ET")
  expect_equal(dim(tr$n_pts), c(500L, 6L))
  expect_identical(colnames(tr$n_tox), paste0("DL", 1:6))
  expect_true(all(tr$n_tox <= tr$n_pts))
  expect_true(all(tr$n_resp <= tr$n_pts))
  expect_identical(rowSums(tr$n_pts), 3 * tr$cohorts_used)
  expect_true(all(tr$stop_reason %in% boin_et_stop_reasons()))
  expect_identical(tr$early_stop,
                   tr$stop_reason %in% c("no_admissible_dose",
                                         "no_lower_admissible_dose"))
  expect_true(all(tr$duration > 0))
  # A full trial lasts at least the arrivals of 36 patients and the longer
  # window after the last one.
  full <- tr$cohorts_used == 12L
  expect_true(all(tr$duration[full] >= 45))
  expect_equal(c(tr$boundary$lambda1, tr$boundary$lambda2, tr$boundary$eta1),
               c(0.14, 0.35, 0.48))
})

test_that("the seed reproduces the trials and the session state is kept", {
  args <- list(target = 0.30, target_resp = 0.60, p_true = c(0.10, 0.20, 0.30),
               p_resp = c(0.2, 0.4, 0.6), n_cohort = 6, cohort_size = 3,
               window = 28, accrual_rate = 0.2, n_trials = 200)
  set.seed(99)
  before <- .Random.seed
  a <- do.call(tite_boin_et_simulate, c(args, seed = 7))
  expect_identical(.Random.seed, before)
  b <- do.call(tite_boin_et_simulate, c(args, seed = 7))
  expect_identical(a$n_pts, b$n_pts)
  expect_identical(a$duration, b$duration)
  c2 <- do.call(tite_boin_et_simulate, c(args, seed = 8))
  expect_false(identical(a$duration, c2$duration))
})

test_that("fixed arrivals without events give the expected timeline", {
  # Arrivals every 10 days, a 30-day window and no events. Cohort 1 arrives
  # at 0, 10 and 20; the decision waits from day 30 until two of the three
  # patients have finished (day 40) and escalates. Cohort 2 arrives at 40,
  # 50 and 60 and the trial ends 30 days after the last arrival.
  tr <- tite_boin_et_simulate(
    target = 0.30, target_resp = 0.60, p_true = c(0, 0), p_resp = c(0, 0),
    n_cohort = 2, cohort_size = 3, window = 30, accrual_rate = 0.1,
    accrual = "fixed", event_time = "uniform", n_trials = 5, seed = 1
  )
  expect_true(all(tr$n_pts[, 1] == 3L & tr$n_pts[, 2] == 3L))
  expect_equal(tr$duration, rep(90, 5))
})

test_that("the dose rules are those of boin_et_simulate()", {
  # With arrivals 10 days apart, a 30-day DLT window and a 60-day response
  # window, every DLT at a dose is known when it is decided on, so that the
  # trials of the deterministic example of test-boin_et_simulate.R recur:
  # no skipping over dose 2, excluded by 3 DLTs, unless boinet_compatible.
  args <- list(target = 0.30, target_resp = 0.60, p_true = c(0, 1, 1),
               p_resp = c(0, 0, 0), n_cohort = 4, cohort_size = 3,
               window = 30, resp_window = 60, accrual_rate = 0.1,
               accrual = "fixed", event_time = "uniform", n_trials = 20,
               seed = 4)
  tr <- do.call(tite_boin_et_simulate, args)
  expect_true(all(tr$n_pts[, 1] == 9L & tr$n_pts[, 2] == 3L &
                    tr$n_pts[, 3] == 0L))
  tr <- do.call(tite_boin_et_simulate, c(args, boinet_compatible = TRUE))
  expect_true(all(tr$n_pts[, 1] == 6L & tr$n_pts[, 2] == 3L &
                    tr$n_pts[, 3] == 3L))

  # No DLT and no response: the excluded highest dose ends the trial.
  args <- list(target = 0.30, target_resp = 0.60, p_true = c(0, 0, 0),
               p_resp = c(0, 0, 0), n_cohort = 12, cohort_size = 3,
               window = 30, accrual_rate = 0.2, n_trials = 20, seed = 5)
  tr <- do.call(tite_boin_et_simulate, args)
  expect_true(all(tr$n_pts[, 3] == 12L))
  expect_true(all(tr$stop_reason == "top_dose_excluded"))
  tr <- do.call(tite_boin_et_simulate, c(args, boinet_compatible = TRUE))
  expect_true(all(tr$n_pts[, 3] == 30L))
  expect_true(all(tr$stop_reason == "max_cohorts"))
})

test_that("the event probabilities and their correlation are honoured", {
  # One dose and no exclusion: 36 patients per trial, all followed to the
  # end. The DLT and response indicators are built from the same copula as
  # in boin_et_simulate(), so the covariance of the counts is again 2.0313
  # for a correlation of 0.5 (standard deviation of the estimate about 0.1).
  base <- list(target = 0.30, target_resp = 0.60, p_true = 0.2, p_resp = 0.5,
               n_cohort = 12, cohort_size = 3, window = 30, resp_window = 60,
               accrual_rate = 0.5, resp_cor = 0.5, n_trials = 4000,
               cutoff_eli = 0.999999, cutoff_eli_resp = 0.999999)
  for (et in c("weibull", "uniform")) {
    tr <- do.call(tite_boin_et_simulate, c(base, event_time = et, seed = 21))
    expect_true(all(tr$n_pts == 36L))
    expect_lt(abs(mean(tr$n_tox) / 36 - 0.2), 0.01)
    expect_lt(abs(mean(tr$n_resp) / 36 - 0.5), 0.01)
    expect_lt(abs(stats::cov(tr$n_tox[, 1], tr$n_resp[, 1]) - 2.0313), 0.5)
  }
})

test_that("invalid input is rejected", {
  args <- list(target = 0.30, target_resp = 0.60, p_true = c(0.1, 0.2),
               p_resp = c(0.2, 0.4), n_cohort = 4, cohort_size = 3,
               window = 30, accrual_rate = 0.1, n_trials = 10)
  expect_error(do.call(tite_boin_et_simulate,
                       modifyList(args, list(window = 0))), "window")
  expect_error(do.call(tite_boin_et_simulate,
                       modifyList(args, list(accrual_rate = -1))),
               "accrual_rate")
  expect_error(do.call(tite_boin_et_simulate,
                       modifyList(args, list(late_fraction = 1))),
               "late_fraction")
  expect_error(do.call(tite_boin_et_simulate,
                       modifyList(args, list(p_resp = c(0.2, 1)))),
               "below 1")
  expect_s3_class(do.call(tite_boin_et_simulate,
                          modifyList(args, list(p_resp = c(0.2, 1),
                                                event_time = "uniform"))),
                  "boin_et_trials")
  expect_error(do.call(tite_boin_et_simulate,
                       modifyList(args, list(accrual = "other"))))
  expect_error(do.call(tite_boin_et_simulate,
                       modifyList(args, list(p_resp = 0.2))), "p_resp")
})
