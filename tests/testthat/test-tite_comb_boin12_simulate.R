# Scenario 4 of Lu et al. (2025), with their utility (wT = 0.11)
sc4 <- list(
  p_true = matrix(c(0.10, 0.20, 0.28, 0.32, 0.18, 0.22, 0.35, 0.45,
                    0.20, 0.24, 0.48, 0.56, 0.28, 0.50, 0.55, 0.60), 4),
  p_resp = matrix(c(0.20, 0.24, 0.28, 0.32, 0.20, 0.25, 0.30, 0.35,
                    0.40, 0.55, 0.58, 0.62, 0.42, 0.58, 0.62, 0.65), 4)
)
lu_psi00 <- 100 * 0.11 / 1.11
lu_psi11 <- 100 / 1.11

test_that("tite_comb_boin12_simulate returns a well formed object", {
  trials <- tite_comb_boin12_simulate(
    target = 0.35, target_resp = 0.2, p_true = sc4$p_true,
    p_resp = sc4$p_resp, n_cohort = 8, cohort_size = 3, window = 1,
    resp_window = 2, accrual_rate = 2, event_time = "uniform",
    n_trials = 200, seed = 1)
  expect_s3_class(trials, "comb_boin12_trials")
  expect_equal(trials$design, "TITE-Comb-BOIN12")
  for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
    expect_equal(dim(trials[[nm]]), c(200L, 4L, 4L))
    expect_type(trials[[nm]], "integer")
  }
  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(trials$n_tox_resp <= pmin(trials$n_tox, trials$n_resp)))
  expect_length(trials$duration, 200L)
  expect_true(all(trials$duration > 0))
  expect_equal(trials$settings$window, 1)
  expect_equal(trials$settings$resp_window, 2)
  expect_equal(trials$settings$min_completed_share, 0.5)
  expect_equal(trials$settings$pending_estimate, "conditional")
  out <- utils::capture.output(res <- print(trials))
  expect_identical(res, trials)
  expect_true(any(grepl("Simulated TITE-Comb-BOIN12 trials", out,
                        fixed = TRUE)))
  expect_true(any(grepl("Average trial duration", out, fixed = TRUE)))
})

test_that("the trials agree with the Python implementation", {
  # Trials of the fixture tite_comb_boin12_trials.csv: the Python
  # implementation fed with R's random numbers after set.seed(20261118) and
  # checked against the C++ core (Validation/comb_boin12/
  # make_tcb12_fixture.py of the development folder).
  fx <- utils::read.csv(test_path("fixtures", "tite_comb_boin12_trials.csv"),
                        stringsAsFactors = FALSE)
  groups <- list(
    G1 = list(tie_break = "higher", no_candidate = "nearest",
              extrasafe = TRUE, accrual = "uniform", event_time = "uniform",
              window = 1, resp_window = 2, accrual_rate = 2,
              min_completed_share = 0.5, pending_estimate = "conditional"),
    G2 = list(tie_break = "random", no_candidate = "stop", extrasafe = FALSE,
              accrual = "exponential", event_time = "weibull", window = 2,
              resp_window = 3, accrual_rate = 2, min_completed_share = 0.5,
              pending_estimate = "marginal"),
    G3 = list(tie_break = "higher", no_candidate = "nearest",
              extrasafe = FALSE, accrual = "uniform", event_time = "uniform",
              window = 1, resp_window = 2, accrual_rate = 1,
              min_completed_share = 0, pending_estimate = "conditional",
              n_earlystop = 9)
  )
  for (g in names(groups)) {
    trials <- do.call(tite_comb_boin12_simulate, c(list(
      target = 0.35, target_resp = 0.2, p_true = sc4$p_true,
      p_resp = sc4$p_resp, n_cohort = 8, cohort_size = 3, n_trials = 30,
      psi00 = lu_psi00, psi11 = lu_psi11, benchmark_weight = 0.6,
      offset = 0.10, seed = 20261118), groups[[g]]))
    exp <- fx[fx$group == g, , drop = FALSE]
    for (v in c("n", "y", "x", "z")) {
      nm <- c(n = "n_pts", y = "n_tox", x = "n_resp", z = "n_tox_resp")[[v]]
      cols <- paste0(v, 1:16)
      expect_equal(unname(matrix(trials[[nm]], 30L)),
                   unname(as.matrix(exp[, cols])), info = paste(g, nm))
    }
    expect_equal(trials$cohorts_used, exp$cohorts_used, info = g)
    expect_equal(trials$stop_reason, exp$stop_reason, info = g)
    expect_equal(trials$duration, exp$duration, tolerance = 1e-12, info = g)
  }
})

test_that("with every outcome observed the trials are those of Comb-BOIN12", {
  # Fixed arrivals draw no variate for the times, and with
  # min_completed_share = 1 no decision sees a pending outcome.
  for (tie in c("higher", "random")) {
    args <- list(target = 0.35, target_resp = 0.2, p_true = sc4$p_true,
                 p_resp = sc4$p_resp, n_cohort = 10, cohort_size = 3,
                 n_trials = 200, psi00 = lu_psi00, psi11 = lu_psi11,
                 benchmark_weight = 0.6, extrasafe = TRUE, offset = 0.10,
                 tie_break = tie, seed = 7)
    comb <- do.call(comb_boin12_simulate, args)
    tite <- do.call(tite_comb_boin12_simulate,
                    c(args, list(window = 2, resp_window = 3,
                                 accrual_rate = 2, accrual = "fixed",
                                 min_completed_share = 1)))
    for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp", "cohorts_used",
                 "stop_reason", "early_stop")) {
      expect_identical(tite[[nm]], comb[[nm]], info = paste(tie, nm))
    }
  }
})

test_that("one row or one column gives the trials of TITE-BOIN12", {
  p_true <- c(0.05, 0.12, 0.27, 0.35, 0.50)
  p_resp <- c(0.20, 0.35, 0.36, 0.37, 0.38)
  args <- list(target = 0.35, target_resp = 0.25, n_cohort = 12,
               cohort_size = 3, window = 45, resp_window = 60,
               accrual_rate = 0.1, n_trials = 200, tie_break = "random",
               no_candidate = "stop", seed = 11)
  single <- do.call(tite_boin12_simulate,
                    c(args, list(p_true = p_true, p_resp = p_resp,
                                 n_explore = 100)))
  for (dims in list(c(1L, 5L), c(5L, 1L))) {
    comb <- do.call(tite_comb_boin12_simulate,
                    c(args, list(p_true = matrix(p_true, dims[1L]),
                                 p_resp = matrix(p_resp, dims[1L]))))
    for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
      expect_identical(unname(matrix(comb[[nm]], 200L)),
                       unname(single[[nm]]))
    }
    expect_identical(comb$stop_reason, single$stop_reason)
    expect_equal(comb$duration, single$duration, tolerance = 1e-12)
  }
})

test_that("invalid arguments are rejected", {
  args <- list(target = 0.35, target_resp = 0.2, p_true = sc4$p_true,
               p_resp = sc4$p_resp, n_cohort = 4, cohort_size = 3,
               window = 1, accrual_rate = 2, n_trials = 5)
  expect_error(do.call(tite_comb_boin12_simulate,
                       c(args, list(min_completed_share = 1.5))),
               "min_completed_share")
  expect_error(do.call(tite_comb_boin12_simulate,
                       c(args[-4], list(p_resp = sc4$p_resp[, 1:3]))))
  expect_error(do.call(tite_comb_boin12_simulate,
                       c(args[-7], list(window = -1))))
  expect_error(do.call(tite_comb_boin12_simulate,
                       c(args, list(start_dose = c(5, 1)))))
})
