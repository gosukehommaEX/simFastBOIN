scenarios <- list(
  "Scenario 16" = list(
    p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2),
    p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)),
  "Scenario 17" = list(
    p_true = matrix(c(0.05, 0.10, 0.15, 0.20, 0.40, 0.50), nrow = 2),
    p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.55), nrow = 2))
)

test_that("sim_comb_boin12_multi matches sim_comb_boin12 per scenario", {
  multi <- sim_comb_boin12_multi(target = 0.33, target_resp = 0.2,
                                 scenarios = scenarios, n_cohort = 14,
                                 cohort_size = 3, n_trials = 100,
                                 est_method = "observed", seed = 1)
  expect_s3_class(multi, "comb_boin12_oc_multi")
  expect_equal(multi$scenario_names, names(scenarios))
  expect_equal(c(multi$n_rows, multi$n_cols), c(2L, 3L))
  single <- sim_comb_boin12(target = 0.33, target_resp = 0.2,
                            p_true = scenarios[[2]]$p_true,
                            p_resp = scenarios[[2]]$p_resp, n_cohort = 14,
                            cohort_size = 3, n_trials = 100,
                            est_method = "observed", seed = 1)
  expect_identical(multi$results[[2]]$sel_percent, single$sel_percent)
  expect_identical(multi$results[[2]]$n_pts_dose, single$n_pts_dose)
  tab <- multi$summary_table
  expect_equal(nrow(tab), 2L)
  expect_equal(tab[["True OBD"]], c("A1B3", "A2B2"))
  expect_equal(tab[["OBD selected (%)"]][2],
               round(single$pct_true_obd_selected, 1))
  expect_equal(tab[["No OBD (%)"]][2], round(single$percent_no_obd, 1))
})

test_that("scenarios without response probabilities are rejected", {
  expect_error(
    sim_comb_boin12_multi(target = 0.33, target_resp = 0.2,
                          scenarios = list(scenarios[[1]]$p_true),
                          n_cohort = 14, cohort_size = 3, n_trials = 10),
    "p_resp"
  )
  bad <- scenarios
  bad[[2]]$p_resp <- matrix(0.5, 3, 3)
  expect_error(
    sim_comb_boin12_multi(target = 0.33, target_resp = 0.2, scenarios = bad,
                          n_cohort = 14, cohort_size = 3, n_trials = 10),
    "p_resp"
  )
})
