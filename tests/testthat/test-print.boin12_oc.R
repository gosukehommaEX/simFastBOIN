test_that("the operating characteristics print as a table", {
  oc <- sim_boin12(target = 0.35, target_resp = 0.25,
                   p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                   p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
                   n_cohort = 12, cohort_size = 3, n_trials = 50, seed = 1)
  expect_output(print(oc), "BOIN12 operating characteristics")
  expect_output(print(oc), "OBD selected")
  expect_output(print(oc, percent = TRUE), "Patients treated \\(%\\)")
  out <- capture.output(res <- print(oc))
  expect_identical(res, oc)
  skip_if_not_installed("knitr")
  expect_output(print(oc, kable = TRUE), "\\|")
})
