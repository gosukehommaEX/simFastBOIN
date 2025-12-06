#' Test for get_pts_and_tox Function
#'
#' @description
#'   Test suite for the get_pts_and_tox function which generates patient
#'   enrollment and toxicity data for BOIN simulations.
#'
#' @details
#'   Tests include:
#'   - Correct data structure generation
#'   - Reproducibility with seed
#'   - Valid toxicity counts (n_tox <= n_pts)
#'   - Matrix dimensions consistency
#'   - Custom p_saf and p_tox parameters
#'
#' @importFrom testthat test_that expect_type expect_true expect_equal

# Test for get_pts_and_tox function
test_that("get_pts_and_tox generates correct data structure", {
  result <- get_pts_and_tox(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    n_earlystop = 18,
    cutoff_eli = 0.95,
    extrasafe = FALSE,
    offset = 0.05,
    n_earlystop_rule = "with_stay",
    titration = FALSE,
    seed = 123
  )

  expect_type(result, "list")
  expect_true(all(c("n_pts_all", "n_tox_all", "eliminated_mat",
                    "cohorts_completed", "stop_reason") %in% names(result)))
  expect_equal(dim(result$n_pts_all), c(10, 5))  # n_trials x n_doses
  expect_equal(dim(result$n_tox_all), c(10, 5))
  expect_equal(dim(result$eliminated_mat), c(10, 5))
})

test_that("get_pts_and_tox respects seed for reproducibility", {
  result1 <- get_pts_and_tox(
    n_trials = 100,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  result2 <- get_pts_and_tox(
    n_trials = 100,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  expect_equal(result1$n_pts_all, result2$n_pts_all)
  expect_equal(result1$n_tox_all, result2$n_tox_all)
  expect_equal(result1$stop_reason, result2$stop_reason)
})

test_that("get_pts_and_tox generates valid toxicity counts", {
  result <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  # Toxicities should not exceed number of patients
  expect_true(all(result$n_tox_all <= result$n_pts_all))

  # Both should be non-negative
  expect_true(all(result$n_pts_all >= 0))
  expect_true(all(result$n_tox_all >= 0))
})

test_that("get_pts_and_tox handles extrasafe parameter", {
  result_safe <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.60, 0.70, 0.80, 0.90, 0.95),  # All toxic scenario
    n_cohort = 10,
    cohort_size = 3,
    extrasafe = TRUE,
    offset = 0.05,
    seed = 123
  )

  result_nosafe <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.60, 0.70, 0.80, 0.90, 0.95),
    n_cohort = 10,
    cohort_size = 3,
    extrasafe = FALSE,
    seed = 123
  )

  # With extrasafe, trials should stop earlier on average
  expect_true(mean(result_safe$cohorts_completed) <=
                mean(result_nosafe$cohorts_completed))
})

test_that("get_pts_and_tox handles titration parameter", {
  result_titration <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.05, 0.10, 0.20, 0.30, 0.45),
    n_cohort = 20,
    cohort_size = 3,
    titration = TRUE,
    seed = 123
  )

  result_no_titration <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.05, 0.10, 0.20, 0.30, 0.45),
    n_cohort = 20,
    cohort_size = 3,
    titration = FALSE,
    seed = 123
  )

  # Both should complete successfully
  expect_type(result_titration, "list")
  expect_type(result_no_titration, "list")
})

test_that("get_pts_and_tox handles n_earlystop_rule parameter", {
  result_with_stay <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    n_earlystop = 9,
    n_earlystop_rule = "with_stay",
    seed = 123
  )

  result_simple <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    n_earlystop = 9,
    n_earlystop_rule = "simple",
    seed = 456
  )

  # Both should complete successfully
  expect_type(result_with_stay, "list")
  expect_type(result_simple, "list")
})

test_that("get_pts_and_tox handles custom p_saf and p_tox parameters", {
  # Test with custom p_saf and p_tox
  result_custom <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    p_saf = 0.15,
    p_tox = 0.45,
    seed = 123
  )

  # Test with default p_saf and p_tox
  result_default <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  # Both should complete successfully
  expect_type(result_custom, "list")
  expect_type(result_default, "list")

  # Results should be valid
  expect_true(all(result_custom$n_tox_all <= result_custom$n_pts_all))
  expect_true(all(result_default$n_tox_all <= result_default$n_pts_all))
})

test_that("get_pts_and_tox uses default p_saf and p_tox when not specified", {
  result <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  # Should complete successfully with default values
  expect_type(result, "list")
  expect_true(all(c("n_pts_all", "n_tox_all", "eliminated_mat") %in% names(result)))
})

test_that("get_pts_and_tox p_saf and p_tox parameters are accepted", {
  # Run with different p_saf and p_tox values
  result1 <- get_pts_and_tox(
    n_trials = 50,
    target = 0.25,
    p_true = c(0.05, 0.15, 0.25, 0.35, 0.50),
    n_cohort = 10,
    cohort_size = 3,
    p_saf = 0.10,
    p_tox = 0.35,
    seed = 123
  )

  result2 <- get_pts_and_tox(
    n_trials = 50,
    target = 0.25,
    p_true = c(0.05, 0.15, 0.25, 0.35, 0.50),
    n_cohort = 10,
    cohort_size = 3,
    p_saf = 0.20,
    p_tox = 0.50,
    seed = 456
  )

  # Both should complete successfully
  expect_type(result1, "list")
  expect_type(result2, "list")

  # Verify valid results
  expect_true(all(result1$n_tox_all <= result1$n_pts_all))
  expect_true(all(result2$n_tox_all <= result2$n_pts_all))
})

test_that("get_pts_and_tox handles edge cases with p_saf and p_tox", {
  # Test with p_saf close to target
  result_close <- get_pts_and_tox(
    n_trials = 20,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    p_saf = 0.25,
    p_tox = 0.35,
    seed = 123
  )

  expect_type(result_close, "list")
  expect_true(all(result_close$n_tox_all <= result_close$n_pts_all))
})
