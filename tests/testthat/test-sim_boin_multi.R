#' Test for sim_boin_multi Function
#'
#' @description
#'   Test suite for the sim_boin_multi function which runs BOIN simulations
#'   across multiple dose-toxicity scenarios.
#'
#' @details
#'   Tests include:
#'   - Successful execution with multiple scenarios
#'   - Correct parameter passing
#'   - Proper summary structure for multi-scenario results
#'
#' @importFrom testthat test_that expect_type expect_named expect_length
#'   expect_s3_class expect_true expect_equal

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
  expect_named(result, c("results_by_scenario", "combined_summary_df",
                         "scenario_names", "n_doses", "call"))
  expect_length(result$results_by_scenario, 2)
  expect_s3_class(result, "boin_multi_summary")
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

  expect_s3_class(result, "boin_multi_summary")
  expect_true("combined_summary_df" %in% names(result))
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

  # Check that combined_summary_df exists and has "Scenario" column
  expect_true("combined_summary_df" %in% names(result))
  expect_true(is.data.frame(result$combined_summary_df))
  expect_true("Scenario" %in% colnames(result$combined_summary_df))

  # Check that we have results for both scenarios
  expect_equal(length(result$scenario_names), 2)
})

test_that("sim_boin_multi handles different numbers of scenarios", {
  # Test with 3 scenarios
  scenarios <- list(
    list(name = "Scenario 1", p_true = c(0.05, 0.10, 0.20, 0.30, 0.45)),
    list(name = "Scenario 2", p_true = c(0.10, 0.15, 0.30, 0.45, 0.60)),
    list(name = "Scenario 3", p_true = c(0.15, 0.25, 0.35, 0.50, 0.70))
  )

  result <- sim_boin_multi(
    scenarios = scenarios,
    target = 0.30,
    n_trials = 10,
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  expect_length(result$results_by_scenario, 3)
  expect_equal(length(result$scenario_names), 3)
})

test_that("sim_boin_multi respects seed for reproducibility", {
  scenarios <- list(
    list(name = "Scenario 1", p_true = c(0.05, 0.10, 0.20, 0.30, 0.45))
  )

  result1 <- sim_boin_multi(
    scenarios = scenarios,
    target = 0.30,
    n_trials = 100,
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  result2 <- sim_boin_multi(
    scenarios = scenarios,
    target = 0.30,
    n_trials = 100,
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  # Results should be identical
  expect_equal(result1$results_by_scenario[[1]]$summary$mtd_selection_percent,
               result2$results_by_scenario[[1]]$summary$mtd_selection_percent)
})

test_that("sim_boin_multi verbose parameter controls output", {
  scenarios <- list(
    list(name = "Test Scenario", p_true = c(0.10, 0.25, 0.40, 0.55, 0.70))
  )

  # Test verbose = FALSE (default, silent mode)
  expect_silent({
    result_silent <- sim_boin_multi(
      scenarios = scenarios,
      target = 0.30,
      n_trials = 10,
      n_cohort = 10,
      cohort_size = 3,
      verbose = FALSE,
      seed = 123
    )
  })

  # Test verbose = TRUE (with output)
  expect_output({
    result_verbose <- sim_boin_multi(
      scenarios = scenarios,
      target = 0.30,
      n_trials = 10,
      n_cohort = 10,
      cohort_size = 3,
      verbose = TRUE,
      seed = 123
    )
  }, "Multi-Scenario Simulation")

  # Results should be identical regardless of verbose setting
  result_silent <- sim_boin_multi(
    scenarios = scenarios,
    target = 0.30,
    n_trials = 10,
    n_cohort = 10,
    cohort_size = 3,
    verbose = FALSE,
    seed = 123
  )

  result_verbose_test <- sim_boin_multi(
    scenarios = scenarios,
    target = 0.30,
    n_trials = 10,
    n_cohort = 10,
    cohort_size = 3,
    verbose = TRUE,
    seed = 123
  )

  expect_equal(
    result_silent$results_by_scenario[[1]]$summary$mtd_selection_percent,
    result_verbose_test$results_by_scenario[[1]]$summary$mtd_selection_percent
  )
})
