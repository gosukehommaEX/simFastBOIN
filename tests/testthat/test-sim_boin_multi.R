# tests/testthat/test-sim_boin_multi.R

# Test for sim_boin_multi function
test_that("sim_boin_multi runs successfully", {
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

  expect_type(result, "list")
  expect_named(result, c("results", "summary"))
  expect_length(result$results, 2)
  expect_s3_class(result$summary, "boin_multi_summary")
})

test_that("sim_boin_multi passes parameters correctly", {
  scenarios <- list(
    list(name = "Test", p_true = c(0.10, 0.25, 0.40, 0.55, 0.70))
  )

  result <- sim_boin_multi(
    scenarios = scenarios,
    target = 0.30,
    n_trials = 10,
    n_cohort = 10,
    cohort_size = 3,
    extrasafe = TRUE,
    boundMTD = TRUE,
    seed = 123
  )

  expect_s3_class(result$summary, "boin_multi_summary")
})

test_that("sim_boin_multi summary has correct structure", {
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

  expect_true("scenario_name" %in% names(result$summary))
  expect_equal(nrow(result$summary$scenario_name), 2)
})
