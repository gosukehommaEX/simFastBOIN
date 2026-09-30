test_that("print.tite_boin_oc_multi shows the design and every scenario", {
  scenarios <- list(
    list(name = "MTD at dose 3", p_true = c(0.05, 0.15, 0.30, 0.45, 0.60)),
    list(name = "All toxic",     p_true = c(0.35, 0.45, 0.55, 0.65, 0.75))
  )
  oc <- sim_tite_boin_multi(
    target = 0.30, scenarios = scenarios, n_cohort = 10, cohort_size = 3,
    window = 3, accrual_rate = 2, method = "ess", n_trials = 100, seed = 123
  )
  out <- capture.output(print(oc))

  expect_identical(out[1L],
                   "TITE-BOIN operating characteristics across 2 scenarios")
  expect_true(any(grepl("effective sample size", out, fixed = TRUE)))
  expect_true(any(grepl("MTD at dose 3", out, fixed = TRUE)))
  expect_true(any(grepl("All toxic", out, fixed = TRUE)))
  expect_true(any(grepl("Trial duration (mean)", out, fixed = TRUE)))
  expect_true(any(grepl("Accrual suspended (% trials)", out, fixed = TRUE)))

  shares <- capture.output(print(oc, percent = TRUE))
  expect_true(any(grepl("Patients treated (%)", shares, fixed = TRUE)))
  expect_output(expect_error(print(oc, percent = "yes"), "percent"), NA)

  result <- quiet_print(oc)
  expect_false(result$visible)
  expect_identical(result$value, oc)
})
