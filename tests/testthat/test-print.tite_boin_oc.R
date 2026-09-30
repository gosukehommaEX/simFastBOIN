test_that("print.tite_boin_oc shows the design, the table and the timing", {
  oc <- sim_tite_boin(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 200, seed = 123
  )
  out <- capture.output(print(oc))

  expect_identical(out[1L], "TITE-BOIN operating characteristics")
  expect_true(any(grepl("single mean imputation", out, fixed = TRUE)))
  expect_true(any(grepl("MTD selected (%)", out, fixed = TRUE)))
  expect_true(any(grepl("Patients treated there", out, fixed = TRUE)))
  expect_true(any(grepl(paste0("Average trial duration              : ",
                               format(round(oc$duration_mean, 1))),
                        out, fixed = TRUE)))
  expect_true(any(grepl(paste0("Trials with accrual suspended       : ",
                               format(round(oc$pct_trials_suspended, 1)), "%"),
                        out, fixed = TRUE)))
})

test_that("percent changes the table and is validated before printing", {
  oc <- sim_tite_boin(
    target = 0.30, p_true = c(0.05, 0.15, 0.30), n_cohort = 5,
    cohort_size = 3, window = 3, accrual_rate = 2, n_trials = 100, seed = 1
  )

  counts <- capture.output(print(oc))
  shares <- capture.output(print(oc, percent = TRUE))
  expect_true(any(grepl("Patients treated (%)", shares, fixed = TRUE)))
  expect_false(any(grepl("Patients treated (%)", counts, fixed = TRUE)))
  expect_output(expect_error(print(oc, percent = NA), "percent"), NA)
})

test_that("print.tite_boin_oc returns its input invisibly", {
  oc <- sim_tite_boin(
    target = 0.30, p_true = c(0.10, 0.25, 0.40), n_cohort = 5,
    cohort_size = 3, window = 3, accrual_rate = 2, n_trials = 20, seed = 2
  )
  result <- quiet_print(oc)

  expect_false(result$visible)
  expect_identical(result$value, oc)
})

test_that("kable output is available", {
  skip_if_not_installed("knitr")
  expect_true(have_package("knitr"))

  oc <- sim_tite_boin(
    target = 0.30, p_true = c(0.10, 0.25, 0.40), n_cohort = 5,
    cohort_size = 3, window = 3, accrual_rate = 2, n_trials = 20, seed = 2
  )
  out <- capture.output(print(oc, kable = TRUE))

  expect_true(any(grepl("|", out, fixed = TRUE)))
  expect_false(any(grepl("TITE-BOIN operating characteristics", out,
                         fixed = TRUE)))
})
