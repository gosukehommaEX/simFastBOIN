test_that("without any response the results are those of sim_tite_boin", {
  common <- list(target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
                 n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
                 max_pending_ratio = 0.49, min_follow_up = 0.25,
                 n_earlystop = 9, bound_mtd = TRUE, overdose_cutoff = 0.3,
                 n_trials = 1000, seed = 8)
  be <- do.call(sim_be_boin, c(common, list(p_resp = rep(0, 5))))
  tite <- do.call(sim_tite_boin, c(common, list(method = "imputation")))

  for (component in c("sel_percent", "percent_no_mtd", "n_pts_dose",
                      "n_tox_dose", "total_n_pts", "total_n_tox", "overdose",
                      "stop_reason_percent", "duration_mean", "duration_sd",
                      "pct_trials_suspended", "avg_n_suspensions",
                      "avg_time_suspended")) {
    expect_equal(be[[component]], tite[[component]], info = component)
  }
  expect_equal(be$total_n_bf, 0)
})

test_that("the summaries agree with the trial by trial data", {
  oc <- sim_be_boin(
    target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
    p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58), n_cohort = 10, cohort_size = 3,
    window = 3, accrual_rate = 2, n_trials = 500, keep_trials = TRUE,
    no_slot = "leave", seed = 9
  )
  trials <- oc$trials

  expect_s3_class(oc, "be_boin_oc")
  expect_s3_class(oc, "backfill_oc")
  expect_s3_class(oc, "boin_oc")
  expect_gt(oc$total_n_bf, 0)
  expect_equal(unname(oc$n_bf_dose), unname(colMeans(trials$n_bf)))
  expect_equal(oc$avg_n_turned_away, mean(trials$n_turned_away))
  expect_equal(oc$avg_time_suspended, mean(trials$time_suspended))
  expect_identical(oc$settings$conflict_dose, "lowest")
  expect_equal(oc$settings$min_follow_up, 0.25)
  expect_true(is.call(oc$call))
})
