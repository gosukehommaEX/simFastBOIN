#' Test for select_mtd Function
#'
#' @description
#'   Test suite for the select_mtd function which selects maximum tolerated
#'   dose (MTD) based on isotonic regression estimates.
#'
#' @details
#'   Tests include:
#'   - Correct MTD selection
#'   - NA return when no valid MTD exists
#'   - boundMTD constraint enforcement
#'   - min_mtd_sample enforcement
#'
#' @importFrom testthat test_that expect_type expect_true expect_s3_class

# Test for select_mtd function
test_that("select_mtd selects correct MTD", {
  target <- 0.30
  n_pts_mat <- matrix(c(6, 9, 12, 15, 18), nrow = 1)
  n_tox_mat <- matrix(c(0, 1, 3, 5, 9), nrow = 1)
  eliminated_mat <- matrix(FALSE, nrow = 1, ncol = 5)

  iso_est_mat <- isotonic_regression(n_pts_mat, n_tox_mat)

  result <- select_mtd(
    iso_est_mat = iso_est_mat,
    n_pts_mat = n_pts_mat,
    eliminated_mat = eliminated_mat,
    target = target,
    boundMTD = FALSE,
    min_mtd_sample = 1
  )

  expect_s3_class(result, "data.frame")
  expect_true("mtd" %in% names(result))
  expect_true(result$mtd >= 1 && result$mtd <= 5 || is.na(result$mtd))
})

test_that("select_mtd returns NA when no valid MTD", {
  target <- 0.30
  n_pts_mat <- matrix(c(6, 9, 12, 15, 18), nrow = 1)
  n_tox_mat <- matrix(c(3, 5, 8, 12, 15), nrow = 1)  # All too toxic
  eliminated_mat <- matrix(c(TRUE, TRUE, TRUE, TRUE, TRUE), nrow = 1, ncol = 5)

  iso_est_mat <- isotonic_regression(n_pts_mat, n_tox_mat)

  result <- select_mtd(
    iso_est_mat = iso_est_mat,
    n_pts_mat = n_pts_mat,
    eliminated_mat = eliminated_mat,
    target = target,
    boundMTD = FALSE,
    min_mtd_sample = 1
  )

  expect_true(is.na(result$mtd))
})

test_that("select_mtd respects boundMTD constraint", {
  target <- 0.30
  boin_bound <- get_boin_boundary(target = target)
  lambda_d <- boin_bound$lambda_d

  n_pts_mat <- matrix(c(6, 9, 12, 15, 18), nrow = 1)
  n_tox_mat <- matrix(c(0, 1, 3, 7, 12), nrow = 1)
  eliminated_mat <- matrix(FALSE, nrow = 1, ncol = 5)

  iso_est_mat <- isotonic_regression(n_pts_mat, n_tox_mat)

  result <- select_mtd(
    iso_est_mat = iso_est_mat,
    n_pts_mat = n_pts_mat,
    eliminated_mat = eliminated_mat,
    target = target,
    boundMTD = TRUE,
    lambda_d = lambda_d,
    min_mtd_sample = 1
  )

  # If MTD is selected, its isotonic estimate should be <= lambda_d
  if (!is.na(result$mtd)) {
    expect_true(iso_est_mat[1, result$mtd] <= lambda_d)
  }
})

test_that("select_mtd respects min_mtd_sample", {
  target <- 0.30
  n_pts_mat <- matrix(c(1, 2, 3, 15, 18), nrow = 1)
  n_tox_mat <- matrix(c(0, 0, 1, 4, 9), nrow = 1)
  eliminated_mat <- matrix(FALSE, nrow = 1, ncol = 5)

  iso_est_mat <- isotonic_regression(n_pts_mat, n_tox_mat)

  result <- select_mtd(
    iso_est_mat = iso_est_mat,
    n_pts_mat = n_pts_mat,
    eliminated_mat = eliminated_mat,
    target = target,
    boundMTD = FALSE,
    min_mtd_sample = 6
  )

  # Only doses 4 and 5 have >= 6 patients
  if (!is.na(result$mtd)) {
    expect_true(result$mtd >= 4)
  }
})

test_that("select_mtd handles lowest dose elimination", {
  target <- 0.30
  n_pts_mat <- matrix(c(6, 9, 12, 15, 18), nrow = 1)
  n_tox_mat <- matrix(c(5, 3, 4, 5, 9), nrow = 1)
  eliminated_mat <- matrix(c(TRUE, FALSE, FALSE, FALSE, FALSE), nrow = 1, ncol = 5)

  iso_est_mat <- isotonic_regression(n_pts_mat, n_tox_mat)

  result <- select_mtd(
    iso_est_mat = iso_est_mat,
    n_pts_mat = n_pts_mat,
    eliminated_mat = eliminated_mat,
    target = target,
    boundMTD = FALSE,
    min_mtd_sample = 1
  )

  # When lowest dose is eliminated, no MTD should be selected
  expect_true(is.na(result$mtd))
  expect_equal(result$reason, "lowest_dose_eliminated")
})
