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
    p_true = c(0.60, 0.70, 0.80, 0.90, 0.95),  # All toxic scenario
    n_cohort = 10,
    cohort_size = 3,
    extrasafe = FALSE,
    seed = 123
  )

  # With extrasafe, trials should stop earlier in toxic scenarios
  expect_true(mean(result_safe$cohorts_completed) <=
                mean(result_nosafe$cohorts_completed))
})

test_that("get_pts_and_tox handles different cohort sizes", {
  # Variable cohort size
  result <- get_pts_and_tox(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = c(1, 3, 3, 3, 3, 3, 3, 3, 3, 3),  # First cohort size 1
    titration = FALSE,
    seed = 123
  )

  expect_equal(dim(result$n_pts_all), c(10, 5))
})

test_that("get_pts_and_tox stop_reason values are valid", {
  result <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  valid_reasons <- c("max_cohorts_reached", "lowest_dose_eliminated",
                     "lowest_dose_too_toxic", "all_doses_eliminated",
                     "n_earlystop_with_stay", "n_earlystop_simple",
                     "max_sample_size_reached")

  expect_true(all(result$stop_reason %in% valid_reasons))
})
