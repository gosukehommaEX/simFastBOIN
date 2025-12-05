# tests/testthat/test-get_boin_boundary.R

# Test for get_boin_boundary function
test_that("get_boin_boundary calculates correct boundaries", {
  # Test with default parameters
  result <- get_boin_boundary(target = 0.30)

  expect_type(result, "list")
  expect_named(result, c("lambda_e", "lambda_d"))
  expect_true(result$lambda_e < result$lambda_d)
  expect_true(result$lambda_e > 0 && result$lambda_e < 1)
  expect_true(result$lambda_d > 0 && result$lambda_d < 1)
})

test_that("get_boin_boundary handles custom p_saf and p_tox", {
  result <- get_boin_boundary(
    target = 0.25,
    p_saf = 0.12,
    p_tox = 0.40
  )

  expect_type(result, "list")
  expect_true(result$lambda_e < 0.25)
  expect_true(result$lambda_d > 0.25)
})

test_that("get_boin_boundary works for different target values", {
  targets <- c(0.20, 0.25, 0.30, 0.35)

  for (target in targets) {
    result <- get_boin_boundary(target = target)
    expect_true(result$lambda_e < target)
    expect_true(result$lambda_d > target)
  }
})
