comb_print_example <- function(n_trials = 200, seed = 123) {
  sim_comb_boin(
    target = 0.25,
    p_true = matrix(c(0.02, 0.04, 0.08, 0.14,
                      0.08, 0.25, 0.42, 0.48,
                      0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE),
    n_cohort = 16, cohort_size = 3, n_trials = n_trials, seed = seed
  )
}

test_that("print.comb_boin_oc shows the design, the table and the summaries", {
  oc <- comb_print_example()
  out <- capture.output(print(oc))

  expect_identical(out[1L], "BOIN drug-combination operating characteristics")
  expect_true(any(grepl("3 of drug A, 4 of drug B", out, fixed = TRUE)))
  expect_true(any(grepl("MTD selected (%)", out, fixed = TRUE)))
  expect_true(any(grepl("Patients with DLT", out, fixed = TRUE)))
  expect_true(any(grepl("A3B1, A2B2", out, fixed = TRUE)))
  selected_line <- out[grepl("Selected as the MTD", out, fixed = TRUE)]
  expect_length(selected_line, 1L)
  expect_true(grepl(paste0(": ", format(round(oc$pct_true_mtd_selected, 1)),
                           "%"), selected_line, fixed = TRUE))
  # Five combinations have a true DLT rate above the target 0.25.
  expect_true(any(grepl("Combinations above a true DLT rate of 0.25: 5 of 12",
                        out, fixed = TRUE)))
})

test_that("the printed table holds the selection percentages", {
  oc <- comb_print_example()
  tab <- comb_oc_table(oc, digits = 1)

  expect_equal(nrow(tab), 4L * 3L)
  expect_equal(unique(tab$Item[tab$Item != ""]),
               c("True DLT rate (%)", "MTD selected (%)", "Patients treated",
                 "Patients with DLT"))
  selected <- tab[4:6, c("B1", "B2", "B3", "B4")]
  expect_equal(unname(as.matrix(selected)),
               unname(round(oc$sel_percent, 1)))
})

test_that("print.comb_boin_oc returns its input invisibly", {
  oc <- comb_print_example(n_trials = 20, seed = 2)
  result <- quiet_print(oc)

  expect_false(result$visible)
  expect_identical(result$value, oc)
})

test_that("kable output is available", {
  skip_if_not_installed("knitr")
  expect_true(have_package("knitr"))

  oc <- comb_print_example(n_trials = 20, seed = 2)
  out <- capture.output(print(oc, kable = TRUE))

  expect_true(any(grepl("|", out, fixed = TRUE)))
  expect_false(any(grepl("BOIN drug-combination operating characteristics",
                         out, fixed = TRUE)))
})
