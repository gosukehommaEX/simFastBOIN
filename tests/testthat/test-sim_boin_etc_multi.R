test_that("sim_boin_etc_multi gives the results of sim_boin_etc per scenario", {
  scenarios <- list(
    "Scenario 1" = list(
      p_true = matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2),
      p_resp = matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)),
    "Scenario 2" = list(
      p_true = matrix(c(0.05, 0.15, 0.15, 0.35, 0.30, 0.55), nrow = 2),
      p_resp = matrix(c(0.25, 0.60, 0.60, 0.65, 0.65, 0.70), nrow = 2))
  )
  res <- sim_boin_etc_multi(target = 0.30, target_resp = 0.60,
                            scenarios = scenarios, n_cohort = 8,
                            cohort_size = 3, n_trials = 200,
                            subtrials = "etc3", seed = 4)
  expect_s3_class(res, "boin_etc_oc_multi")
  expect_equal(res$scenario_names, names(scenarios))
  expect_equal(nrow(res$summary_table), 2L)
  expect_equal(c(res$n_rows, res$n_cols), c(2L, 3L))
  single <- sim_boin_etc(target = 0.30, target_resp = 0.60,
                         p_true = scenarios[[2]]$p_true,
                         p_resp = scenarios[[2]]$p_resp, n_cohort = 8,
                         cohort_size = 3, n_trials = 200, subtrials = "etc3",
                         seed = 4)
  expect_equal(res$results[[2]]$sel_percent, single$sel_percent)
  expect_equal(res$results[[2]]$n_pts_dose, single$n_pts_dose)
  expect_equal(res$summary_table$Patients[2],
               round(single$total_n_pts, 1))
  bad <- list(list(p_true = matrix(0.1, 2, 2)))
  expect_error(sim_boin_etc_multi(target = 0.30, target_resp = 0.60,
                                  scenarios = bad, n_cohort = 8,
                                  cohort_size = 3, n_trials = 10))
})
