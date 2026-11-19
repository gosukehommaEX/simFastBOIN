test_that("print.boin_etc_oc shows the design and the table", {
  p_true <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
  p_resp <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
  oc <- sim_boin_etc(target = 0.30, target_resp = 0.60, p_true = p_true,
                     p_resp = p_resp, n_cohort = 8, cohort_size = 3,
                     n_trials = 50, seed = 1)
  expect_output(print(oc), "BOIN-ETC1 operating characteristics")
  expect_output(print(oc), "start of a row       : diagonal")
  expect_output(print(oc), "lambda1/lambda2/eta1 : 0.14 / 0.35 / 0.48")
  expect_output(print(oc, percent = TRUE), "Patients treated \\(%\\)")
  printed <- capture.output(out <- print(oc))
  expect_true(length(printed) > 10L)
  expect_identical(out, oc)
  skip_if_not_installed("knitr")
  expect_output(print(oc, kable = TRUE), "\\|")
})
