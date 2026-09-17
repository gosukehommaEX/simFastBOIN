example_oc <- function() {
  sim_boin(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, n_trials = 100, seed = 1
  )
}

test_that("print.boin_boundary shows the boundaries", {
  bd <- boin_boundary(target = 0.30, max_n = 18, extrasafe = TRUE)

  expect_output(print(bd), "BOIN decision boundaries")
  expect_output(print(bd), "lambda_e")
  expect_output(print(bd), "Escalate if")
  expect_output(print(bd), "Stop at lowest dose")
  expect_output(print(bd, cohort_size = 3), "Deescalate if")

  returned <- quiet_print(bd)
  expect_false(returned$visible)
  expect_identical(returned$value, bd)
})

test_that("print.boin_boundary restricts the table to cohort ends", {
  bd <- boin_boundary(target = 0.30, max_n = 18)

  full <- capture.output(print(bd))
  reduced <- capture.output(print(bd, cohort_size = 3))

  expect_lt(sum(nchar(reduced)), sum(nchar(full)))
})

test_that("print.boin_boundary rejects a bad cohort size without printing", {
  bd <- boin_boundary(target = 0.30, max_n = 18)

  expect_error(print(bd, cohort_size = 0), "at least 1")
  # The argument is checked before anything is written, so nothing leaks out.
  expect_silent(try(print(bd, cohort_size = 0), silent = TRUE))
})

test_that("print.boin_decision_table labels the two dimensions", {
  decisions <- boin_decision_table(target = 0.30, max_n = 12)

  expect_output(print(decisions), "DLTs")
  expect_output(print(decisions), "Patients")
  expect_output(print(decisions), "E = escalate")
  # Missing combinations are blanked rather than printed as NA.
  expect_false(any(grepl("NA", capture.output(print(decisions)))))

  returned <- quiet_print(decisions)
  expect_false(returned$visible)
  expect_identical(returned$value, decisions)
})

test_that("print.boin_decision_table restricts the table to cohort ends", {
  decisions <- boin_decision_table(target = 0.30, max_n = 18)

  full <- capture.output(print(decisions))
  reduced <- capture.output(print(decisions, cohort_size = 3))

  expect_lt(sum(nchar(reduced)), sum(nchar(full)))
  expect_error(print(decisions, cohort_size = 0), "at least 1")
  expect_silent(try(print(decisions, cohort_size = 0), silent = TRUE))
})

test_that("print.boin_trials summarises rather than dumps the matrices", {
  trials <- boin_simulate(
    target = 0.30, p_true = c(0.10, 0.25, 0.40),
    n_cohort = 8, cohort_size = 3, n_trials = 50, seed = 1
  )

  expect_output(print(trials), "Simulated BOIN trials")
  expect_output(print(trials), "Average per trial")
  expect_output(print(trials), "Stopping reason")
  expect_lt(length(capture.output(print(trials))), 30L)

  returned <- quiet_print(trials)
  expect_false(returned$visible)
  expect_identical(returned$value, trials)
})

test_that("print.boin_oc shows the summary table", {
  oc <- example_oc()

  expect_output(print(oc), "BOIN operating characteristics")
  expect_output(print(oc), "True DLT rate")
  expect_output(print(oc), "MTD selected")
  expect_output(print(oc), "DL1")
  expect_output(print(oc), "Total / No MTD")

  returned <- quiet_print(oc)
  expect_false(returned$visible)
  expect_identical(returned$value, oc)
})

test_that("the percent option relabels the patient rows", {
  oc <- example_oc()

  expect_output(print(oc, percent = TRUE), "Patients treated \\(%\\)")
  expect_output(print(oc, percent = FALSE), "Patients treated")
})

test_that("print.boin_oc can produce a kable table", {
  skip_if_not_installed("knitr")
  expect_true(have_package("knitr"))

  oc <- example_oc()
  expect_output(print(oc, kable = TRUE), "\\|")
  expect_output(print(oc, kable = TRUE, kable_format = "pipe"), "MTD selected")
})

test_that("print.boin_oc_multi shows every scenario", {
  oc <- sim_boin_multi(
    target = 0.30,
    scenarios = list(Conservative = c(0.05, 0.15, 0.30),
                     Aggressive = c(0.30, 0.45, 0.60)),
    n_cohort = 8, cohort_size = 3, n_trials = 100, seed = 1
  )

  out <- capture.output(print(oc))
  expect_true(any(grepl("across 2 scenarios", out)))
  expect_true(any(grepl("Conservative", out)))
  expect_true(any(grepl("Aggressive", out)))

  returned <- quiet_print(oc)
  expect_false(returned$visible)
  expect_identical(returned$value, oc)
})
