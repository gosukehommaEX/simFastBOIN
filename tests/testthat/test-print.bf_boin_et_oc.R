test_that("print.bf_boin_et_oc shows the table and the settings", {
  oc <- sim_bf_boin_et(
    target = 0.30, target_resp = 0.25,
    p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
    p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70), n_cohort = 6, cohort_size = 3,
    window = 1, accrual_rate = 3, n_trials = 50, seed = 1
  )
  out <- capture.output(res <- print(oc))
  expect_identical(res, oc)
  expect_true(any(grepl("BF-BOIN-ET operating characteristics", out)))
  expect_true(any(grepl("Patients backfilled", out)))
  expect_true(any(grepl("Average trial duration", out)))
  out <- capture.output(print(oc, percent = TRUE))
  expect_true(any(grepl("Patients backfilled \\(%\\)", out)))
  expect_error(print(oc, percent = NA), "percent")
  skip_if_not_installed("knitr")
  out <- capture.output(print(oc, kable = TRUE))
  expect_true(any(grepl("^\\|", out)))
})
