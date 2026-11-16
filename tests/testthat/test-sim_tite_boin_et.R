sc1 <- list(p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
            p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70))

test_that("sim_tite_boin_et summarizes the trials and their OBDs", {
  args <- list(target = 0.30, target_resp = 0.60, p_true = sc1$p_true,
               p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
               window = 30, resp_window = 45, accrual_rate = 1 / 10,
               n_trials = 400, seed = 2)
  oc <- do.call(sim_tite_boin_et, c(args, keep_trials = TRUE,
                                    obd_method = "utility_scoring"))
  tr <- do.call(tite_boin_et_simulate, args)
  sel <- boin_et_select_obd(tr$n_pts, tr$n_tox, tr$n_resp, target = 0.30,
                            target_resp = 0.60, obd_method = "utility_scoring",
                            early_stop = tr$early_stop)

  expect_s3_class(oc, "boin_et_oc")
  expect_identical(oc$design, "TITE-BOIN-ET")
  expect_identical(oc$trials$n_pts, tr$n_pts)
  expect_identical(oc$obd, sel$obd)
  expected <- as.numeric(table(factor(sel$obd, levels = 1:6))) / 400 * 100
  expect_equal(unname(oc$sel_percent), expected)
  expect_equal(sum(oc$sel_percent) + oc$percent_no_obd, 100)
  expect_equal(sum(oc$stop_reason_percent), 100)
  expect_equal(oc$duration, mean(tr$duration))
  expect_identical(oc$settings$window, 30)
  expect_identical(oc$settings$resp_window, 45)
  expect_identical(oc$settings$min_completed_share, 0.5)

  # The share reaches the simulation.
  oc0 <- do.call(sim_tite_boin_et, c(args, min_completed_share = 0))
  tr0 <- do.call(tite_boin_et_simulate, c(args, min_completed_share = 0))
  expect_equal(oc0$duration, mean(tr0$duration))
  expect_lt(oc0$duration, oc$duration)
})

test_that("boinet_compatible reaches both the simulation and the selection", {
  args <- list(target = 0.30, target_resp = 0.60, p_true = c(0, 0, 0),
               p_resp = c(0, 0, 0), n_cohort = 12, cohort_size = 3,
               window = 30, accrual_rate = 0.2, n_trials = 30, seed = 5)
  oc <- do.call(sim_tite_boin_et, args)
  expect_equal(unname(oc$n_pts_dose), c(3, 3, 12))
  expect_equal(unname(oc$sel_percent), c(100, 0, 0))
  oc <- do.call(sim_tite_boin_et, c(args, boinet_compatible = TRUE))
  expect_equal(unname(oc$n_pts_dose), c(3, 3, 30))
  expect_true(oc$settings$boinet_compatible)
})

test_that("verbose reports the steps", {
  msgs <- capture_messages(
    sim_tite_boin_et(target = 0.30, target_resp = 0.60,
                     p_true = c(0.1, 0.3), p_resp = c(0.3, 0.6), n_cohort = 4,
                     cohort_size = 3, window = 30, accrual_rate = 0.1,
                     n_trials = 20, verbose = TRUE, seed = 1)
  )
  expect_true(any(grepl("Simulating", msgs)))
  expect_true(any(grepl("Selecting the OBD", msgs)))
})
