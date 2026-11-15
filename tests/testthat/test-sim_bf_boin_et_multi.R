scenarios <- list(
  "Example" = list(p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
                   p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70)),
  "Plateau" = list(p_true = c(0.05, 0.10, 0.15, 0.30, 0.45),
                   p_resp = c(0.10, 0.30, 0.45, 0.50, 0.52))
)

test_that("sim_bf_boin_et_multi runs every scenario with the same design", {
  res <- sim_bf_boin_et_multi(
    target = 0.30, target_resp = 0.25, scenarios = scenarios, n_cohort = 8,
    cohort_size = 3, window = 1, accrual_rate = 3, n_trials = 100, seed = 4
  )
  expect_s3_class(res, "bf_boin_et_oc_multi")
  expect_identical(res$scenario_names, names(scenarios))
  expect_equal(nrow(res$summary_table), 14L)
  expect_identical(res$summary_table$Scenario[c(1, 8)], names(scenarios))
  one <- sim_bf_boin_et(
    target = 0.30, target_resp = 0.25, p_true = scenarios[[2]]$p_true,
    p_resp = scenarios[[2]]$p_resp, n_cohort = 8, cohort_size = 3,
    window = 1, accrual_rate = 3, n_trials = 100, seed = 4
  )
  expect_identical(res$results[[2]]$sel_percent, one$sel_percent)
})

test_that("scenarios must carry DLT and response probabilities", {
  expect_error(
    sim_bf_boin_et_multi(target = 0.30, target_resp = 0.25,
                         scenarios = list(list(p_true = c(0.1, 0.2))),
                         n_cohort = 4, cohort_size = 3, window = 1,
                         accrual_rate = 3, n_trials = 10),
    "p_resp"
  )
})

test_that("verbose reports each scenario", {
  msgs <- capture_messages(
    sim_bf_boin_et_multi(target = 0.30, target_resp = 0.25,
                         scenarios = scenarios, n_cohort = 4, cohort_size = 3,
                         window = 1, accrual_rate = 3, n_trials = 10,
                         verbose = TRUE, seed = 4)
  )
  expect_length(grep("^Scenario", msgs), 2L)
})
