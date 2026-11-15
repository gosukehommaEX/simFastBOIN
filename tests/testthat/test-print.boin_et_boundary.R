test_that("print.boin_et_boundary shows the three boundaries", {
  bd <- boin_et_boundary(target = 0.30, target_resp = 0.60)
  out <- capture.output(res <- print(bd))

  expect_identical(res, bd)
  expect_true(any(grepl("lambda1 \\(escalate\\)   : 0.14", out)))
  expect_true(any(grepl("lambda2 \\(de-escalate\\): 0.35", out)))
  expect_true(any(grepl("eta1 \\(response\\)      : 0.48", out)))
})
