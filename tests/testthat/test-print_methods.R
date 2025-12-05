# tests/testthat/test-print_methods.R

# Test for print methods
test_that("print.boin_summary works without errors", {
  result <- sim_boin(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  expect_output(print(result$summary))
  expect_output(print(result$summary, percent = TRUE))
  expect_output(print(result$summary, kable = TRUE))
})

test_that("print.boin_multi_summary works without errors", {
  scenarios <- list(
    list(name = "Scenario 1", p_true = c(0.05, 0.10, 0.20, 0.30, 0.45)),
    list(name = "Scenario 2", p_true = c(0.10, 0.15, 0.30, 0.45, 0.60))
  )

  result <- sim_boin_multi(
    scenarios = scenarios,
    target = 0.30,
    n_trials = 10,
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  expect_output(print(result$summary))
  expect_output(print(result$summary, kable = TRUE))
})

test_that("print methods handle different kable_format options", {
  result <- sim_boin(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  formats <- c("pipe", "simple", "html", "latex")

  for (format in formats) {
    expect_output(print(result$summary, kable = TRUE, kable_format = format))
  }
})
