test_that("the comparison prints one block per scenario", {
  p <- c(0.10, 0.25, 0.40, 0.55)
  truth <- oc_scenarios(list(A = p, B = p / 2), target = 0.30)
  cmp <- oc_compare("3+3" = list(A = oc_3p3(p), B = oc_3p3(p / 2)),
                    truth = truth)
  expect_output(print(cmp), "Comparison of 1 design\\(s\\)")
  expect_output(print(cmp), "above 0.3")
  expect_output(print(cmp), "3\\+3 \\(MTD\\)")
  out <- capture.output(print(cmp))
  expect_identical(sum(out %in% c("A", "B")), 2L)
  expect_false(any(grepl("Duration", out)))
  expect_output(print(cmp, digits = 3),
                formatC(cmp$pct_correct[1], format = "f", digits = 3),
                fixed = TRUE)
  r <- quiet_print(cmp)
  expect_identical(r$value, cmp)
  expect_false(r$visible)
  expect_error(print.oc_compare(data.frame(a = 1)), "oc_compare")
})
