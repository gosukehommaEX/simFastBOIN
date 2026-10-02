test_that("a scenario matches the same scenario run on its own", {
  scenarios <- list(
    list(name = "Rising", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
         p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58)),
    list(name = "Plateau", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
         p_resp = c(0.30, 0.32, 0.35, 0.36, 0.36))
  )
  multi <- sim_be_boin_multi(
    target = 0.25, scenarios = scenarios, n_cohort = 10, cohort_size = 3,
    window = 3, accrual_rate = 2, n_trials = 300, seed = 42
  )

  expect_s3_class(multi, "be_boin_oc_multi")
  expect_s3_class(multi, "backfill_oc_multi")
  expect_equal(nrow(multi$summary_table), 14L)
  for (i in seq_along(scenarios)) {
    single <- sim_be_boin(
      target = 0.25, p_true = scenarios[[i]]$p_true,
      p_resp = scenarios[[i]]$p_resp, n_cohort = 10, cohort_size = 3,
      window = 3, accrual_rate = 2, n_trials = 300, seed = 42
    )
    expect_equal(multi$results[[i]]$sel_percent, single$sel_percent)
    expect_equal(multi$results[[i]]$n_bf_dose, single$n_bf_dose)
    expect_equal(multi$results[[i]]$duration_mean, single$duration_mean)
  }
})
