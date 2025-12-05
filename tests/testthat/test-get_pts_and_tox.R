# tests/testthat/test-get_pts_and_tox.R

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
  expect_true(all(c("n_pts", "n_tox", "reason", "decision_table", "stop_bound") %in% names(result)))
  expect_equal(dim(result$n_pts), c(10, 5))  # n_trials x n_doses
  expect_equal(dim(result$n_tox), c(10, 5))
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

  expect_equal(result1$n_pts, result2$n_pts)
  expect_equal(result1$n_tox, result2$n_tox)
  expect_equal(result1$reason, result2$reason)
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
  expect_true(all(result$n_tox <= result$n_pts))

  # Both should be non-negative
  expect_true(all(result$n_pts >= 0))
  expect_true(all(result$n_tox >= 0))
})

test_that("get_pts_and_tox handles titration correctly", {
  result_titration <- get_pts_and_tox(
    n_trials = 20,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    titration = TRUE,
    seed = 123
  )

  expect_type(result_titration, "list")
  expect_equal(dim(result_titration$n_pts), c(20, 5))
})

test_that("get_pts_and_tox handles extrasafe stopping", {
  # High toxicity at all doses to trigger extrasafe
  result <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.50, 0.60, 0.70, 0.80, 0.90),
    n_cohort = 10,
    cohort_size = 3,
    extrasafe = TRUE,
    offset = 0.05,
    seed = 123
  )

  # Check that some trials stopped due to safety
  safety_stops <- sum(grepl("lowest_dose_too_toxic|lowest_dose_eliminated", result$reason))
  expect_true(safety_stops > 0)
})

test_that("get_pts_and_tox handles early stopping rules", {
  # Test "simple" rule
  result_simple <- get_pts_and_tox(
    n_trials = 20,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    n_earlystop = 12,
    n_earlystop_rule = "simple",
    seed = 123
  )

  # Test "with_stay" rule
  result_with_stay <- get_pts_and_tox(
    n_trials = 20,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    n_earlystop = 12,
    n_earlystop_rule = "with_stay",
    seed = 123
  )

  expect_type(result_simple$reason, "character")
  expect_type(result_with_stay$reason, "character")
})

test_that("get_pts_and_tox reason field contains valid stopping reasons", {
  result <- get_pts_and_tox(
    n_trials = 100,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  valid_reasons <- c(
    "trial_completed",
    "n_earlystop_reached",
    "n_earlystop_with_stay",
    "lowest_dose_too_toxic",
    "lowest_dose_eliminated",
    "max_cohorts_reached"
  )

  expect_true(all(result$reason %in% valid_reasons))
})

test_that("get_pts_and_tox handles different cohort sizes", {
  cohort_sizes <- c(1, 3, 6)

  for (size in cohort_sizes) {
    result <- get_pts_and_tox(
      n_trials = 10,
      target = 0.30,
      p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
      n_cohort = 10,
      cohort_size = size,
      seed = 123
    )

    expect_type(result, "list")
    expect_equal(dim(result$n_pts), c(10, 5))
  }
})

test_that("get_pts_and_tox handles variable cohort sizes", {
  result <- get_pts_and_tox(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = c(3, 3, 3, 6, 6, 6, 6, 6, 6, 6),
    seed = 123
  )

  expect_type(result, "list")
  expect_equal(dim(result$n_pts), c(10, 5))
})

test_that("get_pts_and_tox decision_table is correctly formatted", {
  result <- get_pts_and_tox(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  expect_true(is.matrix(result$decision_table))

  # Check that all decisions are valid
  valid_decisions <- c("E", "S", "D", "DE", NA_character_)
  unique_decisions <- unique(as.vector(result$decision_table))
  expect_true(all(unique_decisions %in% valid_decisions))
})

test_that("get_pts_and_tox stop_bound is correctly formatted", {
  result <- get_pts_and_tox(
    n_trials = 10,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = 10,
    cohort_size = 3,
    extrasafe = TRUE,
    seed = 123
  )

  expect_true(is.matrix(result$stop_bound))
})

test_that("get_pts_and_tox handles edge case with all toxic doses", {
  result <- get_pts_and_tox(
    n_trials = 20,
    target = 0.30,
    p_true = c(0.60, 0.70, 0.80, 0.90, 0.95),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  # Should have early terminations
  expect_true(any(result$reason != "trial_completed"))
})

test_that("get_pts_and_tox handles edge case with all safe doses", {
  result <- get_pts_and_tox(
    n_trials = 20,
    target = 0.30,
    p_true = c(0.01, 0.03, 0.05, 0.08, 0.10),
    n_cohort = 10,
    cohort_size = 3,
    seed = 123
  )

  expect_type(result, "list")
  # Most patients should be at higher doses
  total_pts_per_dose <- colSums(result$n_pts)
  expect_true(total_pts_per_dose[5] > total_pts_per_dose[1])
})

test_that("get_pts_and_tox total patients respects cohort limits", {
  n_cohort <- 10
  cohort_size <- 3
  max_pts <- n_cohort * cohort_size

  result <- get_pts_and_tox(
    n_trials = 50,
    target = 0.30,
    p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
    n_cohort = n_cohort,
    cohort_size = cohort_size,
    seed = 123
  )

  # Total patients per trial should not exceed max
  total_pts <- rowSums(result$n_pts)
  expect_true(all(total_pts <= max_pts))
})
