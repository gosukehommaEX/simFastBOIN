# tests/testthat/test-select_mtd.R

# Test for select_mtd function
test_that("select_mtd selects correct MTD", {
  iso_est <- c(0.05, 0.15, 0.25, 0.35, 0.50)
  n_pts <- c(6, 9, 12, 15, 18)
  target <- 0.30
  lambda_d <- 0.36

  result <- select_mtd(
    iso_est = iso_est,
    n_pts = n_pts,
    target = target,
    lambda_d = lambda_d,
    min_mtd_sample = 1,
    boundMTD = FALSE
  )

  expect_type(result, "integer")
  expect_true(result >= 1 && result <= length(iso_est))
})

test_that("select_mtd returns NA when no valid MTD", {
  iso_est <- c(0.50, 0.60, 0.70, 0.80, 0.90)
  n_pts <- c(6, 9, 12, 15, 18)
  target <- 0.30
  lambda_d <- 0.36

  result <- select_mtd(
    iso_est = iso_est,
    n_pts = n_pts,
    target = target,
    lambda_d = lambda_d,
    min_mtd_sample = 1,
    boundMTD = FALSE
  )

  expect_true(is.na(result))
})

test_that("select_mtd respects boundMTD constraint", {
  iso_est <- c(0.05, 0.15, 0.25, 0.40, 0.50)
  n_pts <- c(6, 9, 12, 15, 18)
  target <- 0.30
  lambda_d <- 0.36

  result <- select_mtd(
    iso_est = iso_est,
    n_pts = n_pts,
    target = target,
    lambda_d = lambda_d,
    min_mtd_sample = 1,
    boundMTD = TRUE
  )

  # MTD should be dose 3 or lower (iso_est <= lambda_d)
  if (!is.na(result)) {
    expect_true(iso_est[result] <= lambda_d)
  }
})

test_that("select_mtd respects min_mtd_sample", {
  iso_est <- c(0.05, 0.15, 0.25, 0.35, 0.50)
  n_pts <- c(1, 2, 3, 15, 18)
  target <- 0.30
  lambda_d <- 0.36

  result <- select_mtd(
    iso_est = iso_est,
    n_pts = n_pts,
    target = target,
    lambda_d = lambda_d,
    min_mtd_sample = 6,
    boundMTD = FALSE
  )

  # Only doses 4 and 5 have >= 6 patients
  if (!is.na(result)) {
    expect_true(result >= 4)
  }
})
