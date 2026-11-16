test_that("the scenarios print in one table", {
  scenarios <- list(
    "Plateau" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                     p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
    "Increasing" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                        p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
  )
  oc <- sim_boin12_multi(target = 0.35, target_resp = 0.25,
                         scenarios = scenarios, n_cohort = 12,
                         cohort_size = 3, n_trials = 50, seed = 1)
  expect_output(print(oc), "across 2 scenarios")
  expect_output(print(oc), "Increasing")
  expect_output(print(oc, percent = TRUE), "\\(%\\)")
  out <- capture.output(res <- print(oc))
  expect_identical(res, oc)
  expect_error(capture.output(print(oc, percent = NA)), "percent")
})
