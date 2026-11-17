test_that("print.comb_boin12_oc shows the table and the summary", {
  oc <- sim_comb_boin12(
    target = 0.33, target_resp = 0.2,
    p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), 2),
    p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), 2),
    n_cohort = 14, cohort_size = 3, n_trials = 50, est_method = "observed",
    seed = 1)
  out <- capture.output(res <- print(oc))
  expect_identical(res, oc)
  expect_true(any(grepl("Comb-BOIN12 operating characteristics", out)))
  expect_true(any(grepl("ties / no candidate  : higher / nearest", out,
                        fixed = TRUE)))
  expect_true(any(grepl("OBD selection        : observed response, rectangle",
                        out, fixed = TRUE)))
  expect_true(any(grepl("True OBD (highest true utility)      : A1B3", out,
                        fixed = TRUE)))
  expect_true(any(grepl("OBD selected (%)", out, fixed = TRUE)))
  pct <- capture.output(print(oc, percent = TRUE))
  expect_true(any(grepl("Patients treated (%)", pct, fixed = TRUE)))
  skip_if_not_installed("knitr")
  md <- capture.output(print(oc, kable = TRUE))
  expect_true(any(grepl("|Item", md, fixed = TRUE)))
})
