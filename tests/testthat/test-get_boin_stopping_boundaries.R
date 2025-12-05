# tests/testthat/test-get_boin_stopping_boundaries.R

# Test for get_boin_stopping_boundaries function
test_that("get_boin_stopping_boundaries calculates correct boundaries", {
  result <- get_boin_stopping_boundaries(
    target = 0.30,
    max_sample_size = 30,
    cutoff_stop = 0.95
  )

  expect_type(result, "character")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 31)  # 0 to 30 DLTs
  expect_equal(ncol(result), 30)  # 1 to 30 patients
})

test_that("get_boin_stopping_boundaries handles different cutoffs", {
  cutoffs <- c(0.90, 0.95, 0.99)

  for (cutoff in cutoffs) {
    result <- get_boin_stopping_boundaries(
      target = 0.30,
      max_sample_size = 30,
      cutoff_stop = cutoff
    )
    expect_true(is.matrix(result))
    expect_equal(nrow(result), 31)
    expect_equal(ncol(result), 30)
  }
})

test_that("get_boin_stopping_boundaries decisions are valid", {
  result <- get_boin_stopping_boundaries(
    target = 0.30,
    max_sample_size = 30,
    cutoff_stop = 0.95
  )

  valid_decisions <- c("STOP", "GO", NA_character_)
  unique_decisions <- unique(as.vector(result))

  expect_true(all(unique_decisions %in% valid_decisions))
})
