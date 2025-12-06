#' Test for sim_boin Function
#'
#' @description
#'   Test suite for the sim_boin function which runs BOIN trial simulations
#'   and computes operating characteristics.
#'
#' @details
#'   Tests include:
#'   - Successful execution with default parameters
#'   - return_details parameter functionality
#'   - Safety feature handling (extrasafe, boundMTD)
#'   - Reproducibility with seed
#'   - Valid summary statistics
#'   - Custom p_saf and p_tox parameters
#'
#' @importFrom testthat test_that expect_type expect_named expect_s3_class
#'   expect_true expect_equal expect_null expect_length

# Test for sim_boin function
test_that("sim_boin runs successfully with default parameters", {
  result <- sim_boin(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  expect_type(result, "list")
  expect_named(result, c("detailed_results", "summary"))
  expect_s3_class(result$summary, "boin_summary")
})

test_that("sim_boin return_details works correctly", {
  # With details
  result_with <- sim_boin(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    return_details = TRUE,
    seed = 123
  )
  expect_type(result_with$detailed_results, "list")
  expect_length(result_with$detailed_results, 10)

  # Without details
  result_without <- sim_boin(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    return_details = FALSE,
    seed = 123
  )
  expect_null(result_without$detailed_results)
})

test_that("sim_boin handles safety features correctly", {
  result <- sim_boin(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    extrasafe = TRUE,
    offset = 0.05,
    boundMTD = TRUE,
    seed = 123
  )

  expect_s3_class(result$summary, "boin_summary")
  expect_true("avg_n_pts" %in% names(result$summary))
})

test_that("sim_boin respects seed for reproducibility", {
  result1 <- sim_boin(
    n_trials = 100,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  result2 <- sim_boin(
    n_trials = 100,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  expect_equal(result1$summary$mtd_selection_percent,
               result2$summary$mtd_selection_percent)
})

test_that("sim_boin summary statistics are valid", {
  result <- sim_boin(
    n_trials = 100,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  summary <- result$summary

  # Check selection percentages sum to 100
  expect_equal(sum(summary$mtd_selection_percent), 100, tolerance = 0.01)

  # Check average values are non-negative
  expect_true(all(summary$avg_n_pts >= 0))
  expect_true(all(summary$avg_n_tox >= 0))
  expect_true(summary$avg_total_n_pts > 0)
  expect_true(summary$avg_total_n_tox >= 0)
})

test_that("sim_boin handles n_earlystop_rule parameter", {
  result_with_stay <- sim_boin(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    n_earlystop_rule = "with_stay",
    seed = 123
  )

  result_simple <- sim_boin(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    n_earlystop_rule = "simple",
    seed = 456
  )

  # Both should complete successfully
  expect_s3_class(result_with_stay$summary, "boin_summary")
  expect_s3_class(result_simple$summary, "boin_summary")
})

test_that("sim_boin handles titration parameter", {
  result <- sim_boin(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.05, 0.10, 0.20, 0.30, 0.45),
    n_cohort = 20,
    cohort_size = 3,
    titration = TRUE,
    seed = 123
  )

  expect_s3_class(result$summary, "boin_summary")
})

test_that("sim_boin verbose parameter controls output", {
  # Test verbose = FALSE (default, silent mode)
  expect_silent({
    result_silent <- sim_boin(
      n_trials = 10,
      target = 0.30,
      p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
      n_cohort = 10,
      cohort_size = 3,
      verbose = FALSE,
      seed = 123
    )
  })

  # Test verbose = TRUE (with output)
  expect_output({
    result_verbose <- sim_boin(
      n_trials = 10,
      target = 0.30,
      p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
      n_cohort = 10,
      cohort_size = 3,
      verbose = TRUE,
      seed = 123
    )
  }, "BOIN Simulation")

  # Results should be identical regardless of verbose setting
  result_silent <- sim_boin(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    verbose = FALSE,
    seed = 123
  )

  result_verbose_test <- sim_boin(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    verbose = TRUE,
    seed = 123
  )

  expect_equal(result_silent$summary$mtd_selection_percent,
               result_verbose_test$summary$mtd_selection_percent)
})

test_that("sim_boin handles custom p_saf and p_tox parameters", {
  # Test with custom p_saf and p_tox
  result_custom <- sim_boin(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    p_saf = 0.15,
    p_tox = 0.45,
    seed = 123
  )

  expect_s3_class(result_custom$summary, "boin_summary")
  expect_true("p_saf" %in% names(result_custom$summary))
  expect_true("p_tox" %in% names(result_custom$summary))
  expect_equal(result_custom$summary$p_saf, 0.15)
  expect_equal(result_custom$summary$p_tox, 0.45)
})

test_that("sim_boin uses default p_saf and p_tox when not specified", {
  result_default <- sim_boin(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  # Check that default values are used
  expect_equal(result_default$summary$p_saf, 0.6 * 0.30)
  expect_equal(result_default$summary$p_tox, 1.4 * 0.30)
})

test_that("sim_boin p_saf and p_tox are stored correctly in summary", {
  # Run with custom p_saf and p_tox
  result_custom <- sim_boin(
    n_trials = 100,
    target = 0.25,
    p_true = c(0.05, 0.15, 0.25, 0.35, 0.50),
    n_cohort = 10,
    cohort_size = 3,
    p_saf = 0.12,
    p_tox = 0.40,
    seed = 123
  )

  # Verify that custom values are stored correctly
  expect_equal(result_custom$summary$p_saf, 0.12)
  expect_equal(result_custom$summary$p_tox, 0.40)

  # Run with default p_saf and p_tox
  result_default <- sim_boin(
    n_trials = 100,
    target = 0.25,
    p_true = c(0.05, 0.15, 0.25, 0.35, 0.50),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  # Verify that default values are calculated correctly
  expect_equal(result_default$summary$p_saf, 0.6 * 0.25)
  expect_equal(result_default$summary$p_tox, 1.4 * 0.25)
})
