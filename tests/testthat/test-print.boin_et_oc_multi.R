test_that("print.boin_et_oc_multi shows six rows per scenario", {
  scenarios <- list(
    list(p_true = c(0.05, 0.15, 0.30), p_resp = c(0.1, 0.4, 0.6)),
    list(p_true = c(0.10, 0.30, 0.50), p_resp = c(0.5, 0.5, 0.5))
  )
  res <- sim_boin_et_multi(target = 0.30, target_resp = 0.60,
                           scenarios = scenarios, n_cohort = 6,
                           cohort_size = 3, n_trials = 100, seed = 1)
  out <- capture.output(r <- print(res))

  expect_identical(r, res)
  expect_true(any(grepl("across 2 scenarios", out)))
  expect_identical(sum(grepl("True response rate", out)), 2L)
  expect_identical(sum(grepl("OBD selected", out)), 2L)

  out <- capture.output(print(res, percent = TRUE))
  expect_true(any(grepl("Patients with DLT \\(%\\)", out)))
  expect_error(print(res, percent = "yes"), "percent")

  skip_if_not_installed("knitr")
  out <- capture.output(print(res, kable = TRUE))
  expect_true(any(grepl("^\\|", out)))
})
