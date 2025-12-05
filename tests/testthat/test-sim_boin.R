# tests/testthat/test-sim_boin.R

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
