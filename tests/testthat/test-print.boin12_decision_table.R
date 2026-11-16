test_that("the table prints with the scores and the eliminated outcomes", {
  tab <- boin12_decision_table(target = 0.35, target_resp = 0.25,
                               n = c(0, 3))
  expect_output(print(tab), "desirability score table")
  expect_output(print(tab), "E")
  expect_output(print(tab, desirability = TRUE), "0.295")
  out <- capture.output(res <- print(tab))
  expect_identical(res, tab)
  expect_error(capture.output(print(tab, desirability = NA)), "desirability")
})
