scenarios <- list(
  "Scenario 1" = list(p_true = c(0.05, 0.15, 0.25, 0.30),
                      p_resp = c(0.05, 0.35, 0.55, 0.65)),
  "Scenario 2" = list(p_true = c(0.05, 0.15, 0.25, 0.35),
                      p_resp = c(0.30, 0.70, 0.60, 0.50))
)

test_that("sim_tite_boin_et_multi runs every scenario with the same design", {
  res <- sim_tite_boin_et_multi(
    target = 0.30, target_resp = 0.60, scenarios = scenarios, n_cohort = 8,
    cohort_size = 3, window = 30, accrual_rate = 0.1, n_trials = 100,
    seed = 3
  )
  expect_s3_class(res, "boin_et_oc_multi")
  expect_identical(res$scenario_names, names(scenarios))
  expect_equal(nrow(res$summary_table), 12L)
  one <- sim_tite_boin_et(
    target = 0.30, target_resp = 0.60, p_true = scenarios[[2]]$p_true,
    p_resp = scenarios[[2]]$p_resp, n_cohort = 8, cohort_size = 3,
    window = 30, accrual_rate = 0.1, n_trials = 100, seed = 3
  )
  expect_identical(res$results[[2]]$sel_percent, one$sel_percent)
  expect_identical(res$results[[2]]$design, "TITE-BOIN-ET")
  expect_output(print(res), "TITE-BOIN-ET operating characteristics across 2")
})

test_that("verbose reports each scenario", {
  msgs <- capture_messages(
    sim_tite_boin_et_multi(
      target = 0.30, target_resp = 0.60, scenarios = scenarios, n_cohort = 4,
      cohort_size = 3, window = 30, accrual_rate = 0.1, n_trials = 10,
      verbose = TRUE, seed = 3
    )
  )
  expect_length(grep("^Scenario", msgs), 2L)
})
