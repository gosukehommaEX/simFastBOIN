test_that("the scenarios print with the true doses marked", {
  truth <- oc_scenarios(
    list("Scenario 1" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                             p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))),
    target = 0.35, target_resp = 0.25
  )
  expect_output(print(truth), "True scenarios: 1 scenario\\(s\\), 5 doses")
  expect_output(print(truth), "response probability >= 0.25")
  expect_output(print(truth), "utility +50.0 +56.2")
  expect_output(print(truth), "true OBD")
  res <- quiet_print(truth)
  expect_identical(res$value, truth)
  expect_false(res$visible)

  m <- matrix(c(0.05, 0.10, 0.20, 0.10, 0.20, 0.30), nrow = 3)
  comb <- oc_scenarios(list(A = m), target = 0.30)
  expect_output(print(comb), "3 x 2 combinations")
  expect_output(print(comb), "0.30\\*")
  expect_output(print(comb, digits = 3), "0.300\\*")
  expect_error(print.oc_scenarios(data.frame(a = 1)), "oc_scenarios")
  expect_error(print(truth, digits = -1), "at least 0")
})
