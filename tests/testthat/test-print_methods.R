#' Test for Print Methods
#'
#' @description
#'   Test suite for print methods for BOIN summary objects, including
#'   print.boin_summary and print.boin_multi_summary.
#'
#' @details
#'   Tests include:
#'   - Basic print output without errors
#'   - percent parameter functionality
#'   - kable format output
#'   - Different kable_format options (pipe, simple, html, latex)
#'
#' @importFrom testthat test_that expect_output

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

  formats <- c("pipe", "simple")

  for (format in formats) {
    expect_output(print(result$summary, kable = TRUE, kable_format = format))
  }
})

test_that("print.boin_summary percent parameter affects output", {
  result <- sim_boin(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  # Capture output with percent = FALSE (default)
  output_absolute <- capture.output(print(result$summary, percent = FALSE))

  # Capture output with percent = TRUE
  output_percent <- capture.output(print(result$summary, percent = TRUE))

  # Outputs should differ
  expect_false(identical(output_absolute, output_percent))
})

test_that("print methods handle kable parameter correctly", {
  result <- sim_boin(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  # Test with kable = FALSE (default)
  output_no_kable <- capture.output(print(result$summary, kable = FALSE))

  # Test with kable = TRUE
  output_with_kable <- capture.output(print(result$summary, kable = TRUE))

  # Outputs should differ
  expect_false(identical(output_no_kable, output_with_kable))
})
