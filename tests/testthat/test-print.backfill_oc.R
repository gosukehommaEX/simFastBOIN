test_that("print.backfill_oc shows the design, the table and the timing", {
  oc <- sim_bf_boin(
    target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
    p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58), n_cohort = 10, cohort_size = 3,
    window = 1, accrual_rate = 3, n_trials = 200, seed = 1
  )
  out <- capture.output(print(oc))

  expect_match(out[1L], "BF-BOIN operating characteristics", fixed = TRUE)
  expect_true(any(grepl("True response rate (%)", out, fixed = TRUE)))
  expect_true(any(grepl("of which backfilled", out, fixed = TRUE)))
  expect_true(any(grepl("Trials with backfilled patients", out, fixed = TRUE)))
  expect_false(any(grepl("turned away", out, fixed = TRUE)))

  pct <- capture.output(print(oc, percent = TRUE))
  expect_true(any(grepl("of which backfilled (%)", pct, fixed = TRUE)))
  expect_error(print(oc, percent = NA), "percent")

  be <- sim_be_boin(
    target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
    p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58), n_cohort = 10, cohort_size = 3,
    window = 3, accrual_rate = 2, no_slot = "leave", n_trials = 200, seed = 1
  )
  out <- capture.output(print(be))
  expect_match(out[1L], "BE-BOIN operating characteristics", fixed = TRUE)
  expect_true(any(grepl("Average patients turned away", out, fixed = TRUE)))
})

test_that("print.backfill_oc returns its input invisibly and can use kable", {
  skip_if_not_installed("knitr")
  expect_true(have_package("knitr"))

  oc <- sim_bf_boin(
    target = 0.25, p_true = c(0.10, 0.25, 0.40), p_resp = c(0.2, 0.3, 0.4),
    n_cohort = 5, cohort_size = 3, window = 1, accrual_rate = 3,
    n_trials = 50, seed = 2
  )
  result <- quiet_print(oc)
  expect_false(result$visible)
  expect_identical(result$value, oc)

  out <- capture.output(print(oc, kable = TRUE))
  expect_true(any(grepl("True response rate (%)", out, fixed = TRUE)))
  expect_false(any(grepl("operating characteristics", out, fixed = TRUE)))
})
