test_that("print.comb_boin12_oc_multi shows one row per scenario", {
  scenarios <- list(
    "Scenario 16" = list(
      p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2),
      p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)),
    "Scenario 17" = list(
      p_true = matrix(c(0.05, 0.10, 0.15, 0.20, 0.40, 0.50), nrow = 2),
      p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.55), nrow = 2))
  )
  oc <- sim_comb_boin12_multi(target = 0.33, target_resp = 0.2,
                              scenarios = scenarios, n_cohort = 14,
                              cohort_size = 3, n_trials = 50,
                              est_method = "observed", seed = 1)
  out <- capture.output(res <- print(oc))
  expect_identical(res, oc)
  expect_true(any(grepl("across 2 scenarios", out)))
  expect_true(any(grepl("Scenario 17", out)))
  expect_true(any(grepl("A2B2", out)))
  skip_if_not_installed("knitr")
  md <- capture.output(print(oc, kable = TRUE))
  expect_true(any(grepl("|Scenario", md, fixed = TRUE)))
})
