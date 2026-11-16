ex <- list(p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
           p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70))

test_that("sim_bf_boin_et returns consistent operating characteristics", {
  oc <- sim_bf_boin_et(
    target = 0.30, target_resp = 0.25, p_true = ex$p_true, p_resp = ex$p_resp,
    n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
    n_trials = 300, keep_trials = TRUE, seed = 1
  )
  expect_s3_class(oc, "bf_boin_et_oc")
  expect_identical(oc$design, "BF-BOIN-ET")
  expect_identical(names(oc$sel_percent), paste0("DL", 1:5))
  expect_equal(sum(oc$sel_percent) + oc$percent_no_obd, 100)
  expect_equal(sum(oc$stop_reason_percent), 100)
  expect_lte(oc$percent_early_stop, 100)
  # Boundaries of BOIN for a target of 0.30 and the same formula for 0.25,
  # values from an independent Python implementation.
  expect_equal(c(oc$lambda1, oc$lambda2, oc$eta1),
               c(0.23649068523646805, 0.35851946464092954,
                 0.19680087055548712))

  tr <- oc$trials
  expect_true(all(tr$n_bf <= tr$n_pts))
  expect_true(all(tr$n_tox <= tr$n_pts) && all(tr$n_resp <= tr$n_pts))
  expect_true(all(tr$n_tox_bf <= tr$n_tox) && all(tr$n_resp_bf <= tr$n_resp))
  expect_identical(rowSums(tr$n_pts) - rowSums(tr$n_bf),
                   3 * tr$cohorts_used)
  expect_gt(mean(rowSums(tr$n_bf)), 0)
  expect_equal(oc$total_n_bf, mean(rowSums(tr$n_bf)))
  expect_equal(unname(oc$n_pts_dose), unname(colMeans(tr$n_pts)))
  expect_equal(oc$duration, mean(tr$duration))
  expect_equal(unname(oc$sel_percent),
               as.numeric(table(factor(tr$obd, levels = 1:5))) / 300 * 100)
  expect_identical(tr$early_stop, tr$stop_reason == "no_admissible_dose")
})

test_that("degenerate scenarios give the expected trials", {
  # No DLT and no response: escalation to dose 3, which is excluded once it
  # has 9 patients without a response, which ends the trial; nobody is
  # backfilled and dose 1 is selected (independent Python implementation).
  oc <- sim_bf_boin_et(
    target = 0.30, target_resp = 0.25, p_true = c(0, 0, 0),
    p_resp = c(0, 0, 0), n_cohort = 10, cohort_size = 3, window = 1,
    accrual_rate = 3, n_trials = 40, keep_trials = TRUE, seed = 2
  )
  expect_true(all(oc$trials$n_pts[, 1] == 3L & oc$trials$n_pts[, 2] == 3L &
                    oc$trials$n_pts[, 3] == 9L))
  expect_true(all(oc$trials$n_bf == 0L))
  expect_equal(unname(oc$stop_reason_percent["top_dose_excluded"]), 100)
  expect_equal(unname(oc$sel_percent), c(100, 0, 0))

  # Every patient has a DLT: dose 1 is excluded, with every dose above it,
  # once it has 6 patients, an early stop without an OBD.
  oc <- sim_bf_boin_et(
    target = 0.30, target_resp = 0.25, p_true = c(1, 1, 1),
    p_resp = c(0.3, 0.3, 0.3), n_cohort = 10, cohort_size = 3, window = 1,
    accrual_rate = 3, event_time = "uniform", n_trials = 40, seed = 3
  )
  expect_equal(unname(oc$n_pts_dose), c(6, 0, 0))
  expect_equal(oc$percent_early_stop, 100)
  expect_equal(oc$percent_no_obd, 100)
})

test_that("the seed reproduces the trials and the session state is kept", {
  args <- list(target = 0.30, target_resp = 0.25, p_true = ex$p_true,
               p_resp = ex$p_resp, n_cohort = 8, cohort_size = 3, window = 1,
               accrual_rate = 3, n_trials = 100, keep_trials = TRUE)
  set.seed(99)
  before <- .Random.seed
  a <- do.call(sim_bf_boin_et, c(args, seed = 7))
  expect_identical(.Random.seed, before)
  b <- do.call(sim_bf_boin_et, c(args, seed = 7))
  expect_identical(a$trials$n_pts, b$trials$n_pts)
  expect_identical(a$trials$duration, b$trials$duration)
})

test_that("verbose reports the simulation", {
  msgs <- capture_messages(
    sim_bf_boin_et(target = 0.30, target_resp = 0.25, p_true = c(0.1, 0.3),
                   p_resp = c(0.2, 0.5), n_cohort = 4, cohort_size = 3,
                   window = 1, accrual_rate = 3, n_trials = 10,
                   verbose = TRUE, seed = 1)
  )
  expect_true(any(grepl("Simulating", msgs)))
})

test_that("invalid input is rejected", {
  args <- list(target = 0.30, target_resp = 0.25, p_true = c(0.1, 0.2),
               p_resp = c(0.2, 0.4), n_cohort = 4, cohort_size = 3,
               window = 1, accrual_rate = 3, n_trials = 10)
  expect_error(do.call(sim_bf_boin_et, modifyList(args, list(p_saf = 0.35))),
               "p_saf")
  expect_error(do.call(sim_bf_boin_et,
                       modifyList(args, list(p_saf_resp = 0.3))),
               "p_saf_resp")
  expect_error(do.call(sim_bf_boin_et, modifyList(args, list(suspend = 2))),
               "suspend")
  expect_error(do.call(sim_bf_boin_et,
                       modifyList(args, list(backfill_extension = -1))),
               "backfill_extension")
  expect_error(do.call(sim_bf_boin_et,
                       modifyList(args, list(n_elimination = 0))),
               "n_elimination")
  expect_error(do.call(sim_bf_boin_et,
                       modifyList(args, list(psi_resp_no_dlt = NA))),
               "psi_resp_no_dlt")
  expect_error(do.call(sim_bf_boin_et, modifyList(args, list(window = 0))),
               "window")
})
