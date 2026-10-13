test_that("print.comb_boin_oc_multi shows one row per scenario", {
  oc <- sim_comb_boin_multi(
    target = 0.30,
    scenarios = list(
      "MTD in the middle" = matrix(c(0.10, 0.20, 0.30,
                                     0.20, 0.30, 0.45), 2, byrow = TRUE),
      "All toxic"         = matrix(c(0.40, 0.50, 0.60,
                                     0.50, 0.60, 0.70), 2, byrow = TRUE)
    ),
    n_cohort = 8, cohort_size = 3, n_trials = 100, seed = 1
  )
  out <- capture.output(print(oc))

  expect_identical(
    out[1L], "BOIN drug-combination operating characteristics across 2 scenarios"
  )
  expect_true(any(grepl("MTD in the middle", out, fixed = TRUE)))
  expect_true(any(grepl("All toxic", out, fixed = TRUE)))
  expect_true(any(grepl("MTD selected (%)", out, fixed = TRUE)))

  result <- quiet_print(oc)
  expect_false(result$visible)
  expect_identical(result$value, oc)
})

test_that("kable output is available", {
  skip_if_not_installed("knitr")
  expect_true(have_package("knitr"))

  oc <- sim_comb_boin_multi(
    target = 0.30,
    scenarios = list(a = matrix(c(0.10, 0.20, 0.20, 0.30), 2),
                     b = matrix(c(0.30, 0.45, 0.45, 0.60), 2)),
    n_cohort = 6, cohort_size = 3, n_trials = 20, seed = 2
  )
  out <- capture.output(print(oc, kable = TRUE))

  expect_true(any(grepl("|", out, fixed = TRUE)))
  expect_false(any(grepl("across 2 scenarios", out, fixed = TRUE)))
})
