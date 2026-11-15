scenarios <- list(
  list(name = "Monotone", p_true = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.40),
       p_resp = c(0.05, 0.10, 0.25, 0.40, 0.65, 0.70)),
  list(name = "Peaked", p_true = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
       p_resp = c(0.30, 0.70, 0.60, 0.50, 0.40, 0.30))
)

test_that("sim_boin_et_multi runs every scenario with the same design", {
  res <- sim_boin_et_multi(target = 0.30, target_resp = 0.60,
                           scenarios = scenarios, n_cohort = 12,
                           cohort_size = 3, n_trials = 300, seed = 5)

  expect_s3_class(res, "boin_et_oc_multi")
  expect_identical(res$scenario_names, c("Monotone", "Peaked"))
  expect_identical(res$n_doses, 6L)
  expect_equal(nrow(res$summary_table), 12)
  expect_identical(res$summary_table$Scenario[c(1, 7)], c("Monotone", "Peaked"))

  single <- sim_boin_et(target = 0.30, target_resp = 0.60,
                        p_true = scenarios[[2]]$p_true,
                        p_resp = scenarios[[2]]$p_resp, n_cohort = 12,
                        cohort_size = 3, n_trials = 300, seed = 5)
  expect_identical(res$results$Peaked$sel_percent, single$sel_percent)
  expect_identical(res$results$Peaked$n_pts_dose, single$n_pts_dose)
})

test_that("scenarios must carry DLT and response probabilities", {
  expect_error(sim_boin_et_multi(target = 0.30, target_resp = 0.60,
                                 scenarios = list(c(0.1, 0.2)), n_cohort = 4,
                                 cohort_size = 3, n_trials = 10),
               "p_resp")
})

test_that("verbose reports each scenario", {
  msgs <- capture_messages(
    sim_boin_et_multi(target = 0.30, target_resp = 0.60, scenarios = scenarios,
                      n_cohort = 4, cohort_size = 3, n_trials = 20,
                      verbose = TRUE, seed = 1)
  )
  expect_identical(sum(grepl("Scenario [12] of 2", msgs)), 2L)
})
