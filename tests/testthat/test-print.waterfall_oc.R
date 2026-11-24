waterfall_print_example <- function(n_trials = 200, seed = 123) {
  sim_waterfall(
    target = 0.30, p_true = zhang2016_p_true(1), cohort_size = 3,
    n_trials = n_trials, seed = seed
  )
}

test_that("print.waterfall_oc shows the design, the table and the summaries", {
  oc <- waterfall_print_example()
  out <- capture.output(print(oc))

  expect_identical(out[1L],
                   "Waterfall drug-combination operating characteristics")
  expect_true(any(grepl("2 of drug A, 3 of drug B", out, fixed = TRUE)))
  expect_true(any(grepl("at most 6/3 per subtrial, of size 3", out,
                        fixed = TRUE)))
  expect_true(any(grepl("max sample size : 27", out, fixed = TRUE)))
  expect_true(any(grepl("MTD selected (%)", out, fixed = TRUE)))
  expect_true(any(grepl("Patients with DLT", out, fixed = TRUE)))
  contour_line <- out[grepl("True MTD contour", out, fixed = TRUE)]
  expect_length(contour_line, 1L)
  expect_true(grepl(": A1B3, A2B2", contour_line, fixed = TRUE))
  selected_line <- out[grepl("Trials selecting exactly the contour", out,
                             fixed = TRUE)]
  expect_length(selected_line, 1L)
  expect_true(grepl(paste0(": ", format(round(oc$pct_contour_selected, 1)),
                           "%"), selected_line, fixed = TRUE))
})

test_that("the printed table holds the selection percentages", {
  oc <- waterfall_print_example()
  tab <- comb_oc_table(oc, digits = 1)

  expect_equal(nrow(tab), 4L * 2L)
  selected <- tab[3:4, c("B1", "B2", "B3")]
  expect_equal(unname(as.matrix(selected)),
               unname(round(oc$sel_percent, 1)))
})

test_that("print.waterfall_oc returns its input invisibly", {
  oc <- waterfall_print_example(n_trials = 20, seed = 2)
  result <- quiet_print(oc)

  expect_false(result$visible)
  expect_identical(result$value, oc)
})

test_that("kable output is available", {
  skip_if_not_installed("knitr")
  expect_true(have_package("knitr"))

  oc <- waterfall_print_example(n_trials = 20, seed = 2)
  out <- capture.output(print(oc, kable = TRUE))

  expect_true(any(grepl("|", out, fixed = TRUE)))
  expect_false(any(grepl("Waterfall drug-combination operating",
                         out, fixed = TRUE)))
})
