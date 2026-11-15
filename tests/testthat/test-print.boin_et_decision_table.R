test_that("print.boin_et_decision_table shows one matrix per sample size", {
  tab <- boin_et_decision_table(target = 0.30, target_resp = 0.60, n = c(3, 12))
  out <- capture.output(res <- print(tab))

  expect_identical(res, tab)
  expect_true(any(grepl("^3 patients", out)))
  expect_true(any(grepl("^12 patients", out)))
  expect_true(any(grepl("E/S/D", out)))
  # 12 patients without a response: the dose is no longer admissible.
  expect_true(any(grepl("\\*", out)))
})
