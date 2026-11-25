uboin_multi_print_example <- function() {
  sim_uboin_multi(tox_upper = 0.30, resp_lower = 0.20,
                  scenarios = list("Scenario 1" = zhou2019_scenario(1),
                                   "Scenario 4" = zhou2019_scenario(4)),
                  n_cohort = 18, cohort_size = 3, n_trials = 100, seed = 123)
}

test_that("print.uboin_oc_multi shows the design and the summary table", {
  oc <- uboin_multi_print_example()
  out <- capture.output(print(oc))
  expect_identical(out[1L],
                   "U-BOIN operating characteristics across 2 scenarios")
  expect_true(any(grepl("trials each          : 100", out, fixed = TRUE)))
  expect_true(any(grepl("OBD selected (%)", out, fixed = TRUE)))
  expect_true(any(grepl("Scenario 4", out, fixed = TRUE)))
})

test_that("print.uboin_oc_multi returns its input invisibly", {
  oc <- uboin_multi_print_example()
  result <- quiet_print(oc)
  expect_false(result$visible)
  expect_identical(result$value, oc)
})

test_that("kable output is available", {
  skip_if_not_installed("knitr")
  expect_true(have_package("knitr"))
  oc <- uboin_multi_print_example()
  out <- capture.output(print(oc, kable = TRUE))
  expect_true(any(grepl("|", out, fixed = TRUE)))
  expect_false(any(grepl("U-BOIN operating characteristics", out,
                         fixed = TRUE)))
})
