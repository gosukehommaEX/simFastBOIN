test_that("the scenarios are drawn as curves or tiles", {
  skip_if_not_installed("ggplot2")
  expect_true(have_package("ggplot2"))
  truth <- oc_scenarios(
    list(A = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                  p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
         B = list(p_true = c(0.25, 0.55, 0.67, 0.87, 0.95),
                  p_resp = c(0.01, 0.55, 0.55, 0.55, 0.55))),
    target = 0.35, target_resp = 0.25
  )
  p <- plot(truth)
  expect_s3_class(p, "ggplot")
  # Three curves of five doses in two scenarios.
  expect_identical(nrow(p$data), 30L)
  expect_identical(levels(p$data$series),
                   c("DLT probability", "Response probability",
                     "Utility / 100"))
  # One MTD in each scenario and one OBD (none in B).
  marks <- p$layers[[3]]$data
  expect_identical(marks$label, c("MTD", "MTD", "OBD"))

  only_tox <- plot(oc_scenarios(list(c(0.1, 0.3, 0.5)), target = 0.3))
  expect_identical(levels(only_tox$data$series), "DLT probability")

  m <- matrix(c(0.05, 0.10, 0.20, 0.10, 0.20, 0.30), nrow = 3)
  tiles <- plot(oc_scenarios(list(A = m), target = 0.30))
  expect_s3_class(tiles, "ggplot")
  expect_identical(nrow(tiles$data), 6L)
  expect_identical(nrow(tiles$layers[[3]]$data), 1L)

  expect_error(plot(truth, text_size = 0), "positive number")
  expect_error(plot.oc_scenarios(data.frame(a = 1)), "oc_scenarios")
})
