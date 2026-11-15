test_that("print.bf_boin_et_oc_multi shows one block per scenario", {
  scenarios <- list(
    "Example" = list(p_true = c(0.03, 0.05, 0.20), p_resp = c(0.05, 0.40, 0.60)),
    "Flat" = list(p_true = c(0.05, 0.10, 0.15), p_resp = c(0.30, 0.30, 0.30))
  )
  res <- sim_bf_boin_et_multi(
    target = 0.30, target_resp = 0.25, scenarios = scenarios, n_cohort = 6,
    cohort_size = 3, window = 1, accrual_rate = 3, n_trials = 30, seed = 2
  )
  out <- capture.output(val <- print(res, digits = 2))
  expect_identical(val, res)
  expect_true(any(grepl("across 2 scenarios", out)))
  expect_true(any(grepl("Example", out)) && any(grepl("Flat", out)))
  out <- capture.output(print(res, percent = TRUE))
  expect_true(any(grepl("(%)", out, fixed = TRUE)))
  skip_if_not_installed("knitr")
  out <- capture.output(print(res, kable = TRUE))
  expect_true(any(grepl("^\\|", out)))
})
