test_that("print.tite_boin_trials summarizes the trials", {
  trials <- tite_boin_simulate(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 100, seed = 1
  )
  out <- capture.output(print(trials))

  expect_identical(out[1L], "Simulated TITE-BOIN trials")
  expect_true(any(grepl("single mean imputation", out, fixed = TRUE)))
  expect_true(any(grepl("window          : 3", out, fixed = TRUE)))
  expect_true(any(grepl("2 patients per unit of time (exponential)", out,
                        fixed = TRUE)))
  expect_true(any(grepl(paste0("mean ", format(round(mean(trials$duration), 2))),
                        out, fixed = TRUE)))
  expect_true(any(grepl("Trials with suspension", out, fixed = TRUE)))
  expect_true(any(grepl("Stopping reason", out, fixed = TRUE)))

  ess <- tite_boin_simulate(
    target = 0.30, p_true = c(0.05, 0.15, 0.30), n_cohort = 5,
    cohort_size = 3, window = 3, accrual_rate = 2, method = "ess",
    n_trials = 20, seed = 1
  )
  expect_true(any(grepl("effective sample size",
                        capture.output(print(ess)), fixed = TRUE)))
})

test_that("print.tite_boin_trials returns its input invisibly", {
  trials <- tite_boin_simulate(
    target = 0.30, p_true = c(0.10, 0.25, 0.40), n_cohort = 5,
    cohort_size = 3, window = 3, accrual_rate = 2, n_trials = 20, seed = 2
  )
  result <- quiet_print(trials)

  expect_false(result$visible)
  expect_identical(result$value, trials)
})
