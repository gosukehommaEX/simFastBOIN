test_that("print.boin_et_oc shows the table and the settings", {
  oc <- sim_boin_et(target = 0.30, target_resp = 0.60,
                    p_true = c(0.05, 0.15, 0.25, 0.30),
                    p_resp = c(0.05, 0.35, 0.55, 0.65), n_cohort = 8,
                    cohort_size = 3, n_trials = 200, seed = 1)
  out <- capture.output(res <- print(oc))

  expect_identical(res, oc)
  expect_true(any(grepl("BOIN-ET operating characteristics", out)))
  expect_true(any(grepl("OBD selected \\(%\\)", out)))
  expect_true(any(grepl("Patients with response", out)))
  expect_true(any(grepl("Total / No OBD", out)))
  expect_true(any(grepl("Trials stopped early", out)))

  out <- capture.output(print(oc, percent = TRUE))
  expect_true(any(grepl("Patients treated \\(%\\)", out)))
  expect_error(print(oc, percent = NA), "percent")

  skip_if_not_installed("knitr")
  out <- capture.output(print(oc, kable = TRUE))
  expect_true(any(grepl("^\\|", out)))
})
