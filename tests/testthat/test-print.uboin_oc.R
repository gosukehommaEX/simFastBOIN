uboin_print_example <- function(n_trials = 200, seed = 123) {
  sc <- zhou2019_scenario(1)
  sim_uboin(tox_upper = 0.30, resp_lower = 0.20, p_true = sc$p_true,
            p_resp = sc$p_resp, n_cohort = 18, cohort_size = 3,
            resp_assoc = 0.2, n_trials = n_trials, seed = seed)
}

test_that("print.uboin_oc shows the design, the table and the endings", {
  oc <- uboin_print_example()
  out <- capture.output(print(oc))

  expect_identical(out[1L], "U-BOIN operating characteristics")
  expect_true(any(grepl("DLT upper limit      : 0.3", out, fixed = TRUE)))
  expect_true(any(grepl("stage I target       : 0.25", out, fixed = TRUE)))
  expect_true(any(grepl("18 of size 3", out, fixed = TRUE)))
  expect_true(any(grepl("12 / 54 patients at a dose", out, fixed = TRUE)))
  expect_true(any(grepl("utilities            : 0, 30, 50, 100", out,
                        fixed = TRUE)))
  expect_true(any(grepl("True utility", out, fixed = TRUE)))
  expect_true(any(grepl("OBD selected (%)", out, fixed = TRUE)))
  expect_true(any(grepl("Reason for ending (%):", out, fixed = TRUE)))
})

test_that("percent shows the patients as percentages", {
  oc <- uboin_print_example(n_trials = 50, seed = 2)
  out <- capture.output(print(oc, percent = TRUE))
  expect_true(any(grepl("Patients treated (%)", out, fixed = TRUE)))
  expect_error(print(oc, percent = NA), "'percent' must be TRUE or FALSE")
})

test_that("print.uboin_oc returns its input invisibly", {
  oc <- uboin_print_example(n_trials = 20, seed = 2)
  result <- quiet_print(oc)
  expect_false(result$visible)
  expect_identical(result$value, oc)
})

test_that("kable output is available", {
  skip_if_not_installed("knitr")
  expect_true(have_package("knitr"))
  oc <- uboin_print_example(n_trials = 20, seed = 2)
  out <- capture.output(print(oc, kable = TRUE))
  expect_true(any(grepl("|", out, fixed = TRUE)))
  expect_false(any(grepl("U-BOIN operating characteristics", out,
                         fixed = TRUE)))
})
