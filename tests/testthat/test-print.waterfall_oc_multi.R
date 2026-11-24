waterfall_multi_print_example <- function() {
  sim_waterfall_multi(
    target = 0.30,
    scenarios = list("Scenario 1" = zhang2016_p_true(1),
                     "Scenario 4" = zhang2016_p_true(4)),
    cohort_size = 3, n_trials = 100, seed = 123
  )
}

test_that("print.waterfall_oc_multi shows the design and the summary table", {
  oc <- waterfall_multi_print_example()
  out <- capture.output(print(oc))

  expect_identical(
    out[1L],
    "Waterfall drug-combination operating characteristics across 2 scenarios"
  )
  expect_true(any(grepl("2 of drug A, 3 of drug B", out, fixed = TRUE)))
  expect_true(any(grepl("trials each     : 100", out, fixed = TRUE)))
  expect_true(any(grepl("Contour selected (%)", out, fixed = TRUE)))
  expect_true(any(grepl("A1B3, A2B2", out, fixed = TRUE)))
  expect_true(any(grepl("Scenario 4", out, fixed = TRUE)))
})

test_that("print.waterfall_oc_multi returns its input invisibly", {
  oc <- waterfall_multi_print_example()
  result <- quiet_print(oc)

  expect_false(result$visible)
  expect_identical(result$value, oc)
})

test_that("kable output is available", {
  skip_if_not_installed("knitr")
  expect_true(have_package("knitr"))

  oc <- waterfall_multi_print_example()
  out <- capture.output(print(oc, kable = TRUE))

  expect_true(any(grepl("|", out, fixed = TRUE)))
  expect_false(any(grepl("Waterfall drug-combination operating",
                         out, fixed = TRUE)))
})
