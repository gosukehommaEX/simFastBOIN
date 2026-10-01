test_that("plot.decision_table_3p3 draws one cell per state", {
  skip_if_not_installed("ggplot2")
  expect_true(have_package("ggplot2"))

  prev <- decision_table_3p3()
  p <- plot(prev)

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), nrow(prev))
  expect_equal(levels(p$data$decision), c("E", "S", "STOP"))
  expect_equal(levels(p$data$panel), "Dose escalation")
  expect_identical(as.character(p$data$decision), prev$decision)

  expand <- decision_table_3p3(mtd_rule = "expand")
  p <- plot(expand)

  expect_equal(nrow(p$data), nrow(expand))
  expect_equal(levels(p$data$decision), c("E", "S", "STOP", "MTD", "D"))
  expect_equal(levels(p$data$panel), c("Dose escalation", "Search for the MTD"))
  expect_identical(as.character(p$data$decision), expand$decision)
  expect_identical(as.integer(as.character(p$data$n_pts)), expand$n)
})

test_that("invalid arguments are rejected", {
  skip_if_not_installed("ggplot2")

  expand <- decision_table_3p3(mtd_rule = "expand")
  mine <- c(E = "#4DAF4A", S = "#377EB8", STOP = "#E41A1C", MTD = "#984EA3",
            D = "#FF7F00")

  expect_s3_class(plot(expand, colors = mine), "ggplot")
  expect_s3_class(plot(decision_table_3p3(), colors = mine[c("E", "S", "STOP")]),
                  "ggplot")
  expect_error(plot(expand, colors = mine[c("E", "S", "STOP")]), "named")
  expect_error(plot(expand, text_size = 0), "positive number")
  expect_error(plot.decision_table_3p3(data.frame(a = 1)), "decision_table_3p3")
})
