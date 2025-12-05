# tests/testthat/test-isotonic_regression.R

# Test for isotonic_regression function
test_that("isotonic_regression produces valid estimates", {
  n_pts <- c(3, 6, 9, 12, 15)
  n_tox <- c(0, 1, 3, 7, 12)

  result <- isotonic_regression(n_pts, n_tox)

  expect_type(result, "double")
  expect_length(result, length(n_pts))
  expect_true(all(result >= 0 & result <= 1))
  expect_true(all(diff(result) >= 0))  # Non-decreasing
})

test_that("isotonic_regression handles edge cases", {
  # All zeros
  result <- isotonic_regression(c(3, 3, 3), c(0, 0, 0))
  expect_true(all(result == 0))

  # All DLTs
  result <- isotonic_regression(c(3, 3, 3), c(3, 3, 3))
  expect_true(all(result == 1))
})

test_that("isotonic_regression handles single dose", {
  result <- isotonic_regression(n_pts = 6, n_tox = 2)
  expect_equal(result, 2 / 6)
})
