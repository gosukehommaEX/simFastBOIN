test_that("print.boin_et_trials summarizes the trials", {
  tr <- boin_et_simulate(
    target = 0.30, target_resp = 0.60, p_true = c(0.1, 0.2, 0.3),
    p_resp = c(0.2, 0.4, 0.6), n_cohort = 6, cohort_size = 3, n_trials = 100,
    seed = 1
  )
  out <- capture.output(res <- print(tr))

  expect_identical(res, tr)
  expect_true(any(grepl("Simulated BOIN-ET trials", out)))
  expect_true(any(grepl("Responses observed", out)))
  expect_true(any(grepl("Stopping reason", out)))
})
