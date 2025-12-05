#' Test for get_boin_decision Function
#'
#' @description
#'   Test suite for the get_boin_decision function which generates decision
#'   tables for dose escalation and de-escalation.
#'
#' @details
#'   Tests include:
#'   - Correct decision table structure generation
#'   - Different sample sizes
#'   - Valid decision values (E, S, D, DE, NA)
#'
#' @importFrom testthat test_that expect_type expect_true expect_equal

# Test for get_boin_decision function
test_that("get_boin_decision generates correct decision table", {
  boin_bound <- get_boin_boundary(target = 0.30)
  result <- get_boin_decision(
    target = 0.30,
    lambda_e = boin_bound$lambda_e,
    lambda_d = boin_bound$lambda_d,
    max_sample_size = 18,
    cutoff_eli = 0.95
  )

  expect_type(result, "character")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 19)  # 0 to 18 DLTs
  expect_equal(ncol(result), 18)  # 1 to 18 patients
})

test_that("get_boin_decision handles different sample sizes", {
  boin_bound <- get_boin_boundary(target = 0.30)

  for (size in c(9, 18, 30)) {
    result <- get_boin_decision(
      target = 0.30,
      lambda_e = boin_bound$lambda_e,
      lambda_d = boin_bound$lambda_d,
      max_sample_size = size,
      cutoff_eli = 0.95
    )
    expect_equal(nrow(result), size + 1)
    expect_equal(ncol(result), size)
  }
})

test_that("get_boin_decision decisions are valid", {
  boin_bound <- get_boin_boundary(target = 0.30)
  result <- get_boin_decision(
    target = 0.30,
    lambda_e = boin_bound$lambda_e,
    lambda_d = boin_bound$lambda_d,
    max_sample_size = 18,
    cutoff_eli = 0.95
  )

  valid_decisions <- c("E", "S", "D", "DE", NA_character_)
  unique_decisions <- unique(as.vector(result))

  expect_true(all(unique_decisions %in% valid_decisions))
})
