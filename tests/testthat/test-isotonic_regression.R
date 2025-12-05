#' Test for isotonic_regression Function
#'
#' @description
#'   Test suite for the isotonic_regression function which applies isotonic
#'   regression to estimate dose-toxicity curves.
#'
#' @details
#'   Tests include:
#'   - Valid estimates with monotonicity constraint
#'   - Edge cases (all zeros, all DLTs)
#'   - Single dose handling
#'
#' @importFrom testthat test_that expect_type expect_length expect_true expect_equal

# Test for isotonic_regression function
test_that("isotonic_regression produces valid estimates", {
  # Convert to matrix format (1 trial x 5 doses)
  n_pts_mat <- matrix(c(3, 6, 9, 12, 15), nrow = 1)
  n_tox_mat <- matrix(c(0, 1, 3, 7, 12), nrow = 1)

  result <- isotonic_regression(n_pts_mat, n_tox_mat)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 5)
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
  expect_true(all(diff(result[1, ]) >= 0, na.rm = TRUE))  # Non-decreasing
})

test_that("isotonic_regression handles edge cases", {
  # All zeros - convert to matrix format
  n_pts_mat <- matrix(c(3, 3, 3), nrow = 1)
  n_tox_mat <- matrix(c(0, 0, 0), nrow = 1)
  result <- isotonic_regression(n_pts_mat, n_tox_mat)
  expect_true(all(result < 0.1, na.rm = TRUE))  # Near zero with pseudocounts

  # All DLTs - convert to matrix format
  n_pts_mat <- matrix(c(3, 3, 3), nrow = 1)
  n_tox_mat <- matrix(c(3, 3, 3), nrow = 1)
  result <- isotonic_regression(n_pts_mat, n_tox_mat)
  expect_true(all(result > 0.9, na.rm = TRUE))  # Near 1 with pseudocounts
})

test_that("isotonic_regression handles single dose", {
  # Convert to matrix format (1 trial x 1 dose)
  n_pts_mat <- matrix(6, nrow = 1, ncol = 1)
  n_tox_mat <- matrix(2, nrow = 1, ncol = 1)

  result <- isotonic_regression(n_pts_mat, n_tox_mat)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 1)
  expect_true(result[1, 1] >= 0 && result[1, 1] <= 1)
})

test_that("isotonic_regression handles matrix input", {
  n_pts_mat <- matrix(c(3, 6, 9, 3, 6, 9), nrow = 2, byrow = TRUE)
  n_tox_mat <- matrix(c(0, 1, 3, 0, 1, 2), nrow = 2, byrow = TRUE)

  result <- isotonic_regression(n_pts_mat, n_tox_mat)

  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(n_pts_mat))
  expect_true(all(result >= 0 & result <= 1))
})
