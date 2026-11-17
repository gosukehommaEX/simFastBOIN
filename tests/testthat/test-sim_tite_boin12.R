sc1 <- list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
            p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))

test_that("sim_tite_boin12 summarizes the trials, their OBDs and durations", {
  args <- list(target = 0.35, target_resp = 0.25, p_true = sc1$p_true,
               p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
               window = 45, resp_window = 60, accrual_rate = 0.1,
               n_trials = 300, seed = 1)
  oc <- do.call(sim_tite_boin12, c(args, list(keep_trials = TRUE)))
  expect_s3_class(oc, "boin12_oc")
  expect_equal(oc$design, "TITE-BOIN12")
  expect_lt(abs(sum(oc$sel_percent) + oc$percent_no_obd - 100), 1e-10)
  expect_lt(abs(sum(oc$stop_reason_percent) - 100), 1e-10)
  expect_lt(abs(oc$total_n_pts - mean(rowSums(oc$trials$n_pts))), 1e-10)
  expect_equal(oc$duration, mean(oc$trials$duration))
  trials <- do.call(tite_boin12_simulate, args)
  sel <- boin12_select_obd(trials$n_pts, trials$n_tox, trials$n_resp,
                           trials$n_tox_resp, target = 0.35,
                           target_resp = 0.25, early_stop = trials$early_stop)
  expect_identical(oc$obd, sel$obd)
  expect_identical(oc$trials$n_pts, trials$n_pts)
  expect_equal(oc$settings$pending_estimate, "conditional")
  expect_equal(oc$est_method, "posterior_mean")
})

test_that("waiting for every outcome gives the BOIN12 results", {
  common <- list(target = 0.35, target_resp = 0.25, p_true = sc1$p_true,
                 p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
                 n_trials = 300, seed = 2)
  a <- do.call(sim_boin12, common)
  b <- do.call(sim_tite_boin12,
               c(common, list(window = 45, resp_window = 60,
                              accrual_rate = 0.1, accrual = "fixed",
                              min_completed_share = 1)))
  expect_equal(b$sel_percent, a$sel_percent)
  expect_equal(b$n_pts_dose, a$n_pts_dose)
  expect_equal(b$percent_no_obd, a$percent_no_obd)
  expect_null(a$duration)
  expect_gt(b$duration, 0)
})

test_that("the output and the messages", {
  oc <- sim_tite_boin12(target = 0.35, target_resp = 0.25,
                        p_true = sc1$p_true, p_resp = sc1$p_resp,
                        n_cohort = 12, cohort_size = 3, window = 45,
                        resp_window = 60, accrual_rate = 0.1, n_trials = 50,
                        seed = 3)
  expect_output(print(oc), "TITE-BOIN12 operating characteristics")
  expect_output(print(oc), "Average trial duration")
  msgs <- capture_messages(
    sim_tite_boin12(target = 0.35, target_resp = 0.25, p_true = sc1$p_true,
                    p_resp = sc1$p_resp, n_cohort = 4, cohort_size = 3,
                    window = 45, accrual_rate = 0.1, n_trials = 5,
                    verbose = TRUE)
  )
  expect_true(any(grepl("Simulating", msgs)))
  expect_true(any(grepl("Selecting the OBD", msgs)))
  expect_error(sim_tite_boin12(0.35, 0.25, sc1$p_true, sc1$p_resp, 4, 3,
                               window = 45, accrual_rate = 0.1,
                               keep_trials = NA),
               "keep_trials")
})

test_that("the rules for the next dose reach the trials", {
  args <- list(target = 0.35, target_resp = 0.25,
               p_true = c(0.10, 0.20, 0.40, 0.55, 0.70),
               p_resp = c(0.05, 0.10, 0.30, 0.35, 0.40), n_cohort = 12,
               cohort_size = 3, window = 45, resp_window = 60,
               accrual_rate = 0.1, n_trials = 200, tie_break = "random",
               no_candidate = "stop", seed = 3)
  oc <- do.call(sim_tite_boin12, c(args, list(keep_trials = TRUE)))
  trials <- do.call(tite_boin12_simulate, args)
  expect_identical(oc$trials$n_pts, trials$n_pts)
  expect_identical(oc$trials$duration, trials$duration)
  expect_equal(oc$settings$tie_break, "random")
  expect_equal(oc$settings$no_candidate, "stop")
  expect_named(oc$stop_reason_percent, boin12_stop_reasons())
  expect_equal(oc$percent_early_stop, 100 * mean(trials$early_stop))
  expect_gt(unname(oc$stop_reason_percent["no_admissible_candidate"]), 0)
  expect_error(do.call(sim_tite_boin12,
                       modifyList(args, list(no_candidate = "x"))),
               "should be one of")
})
