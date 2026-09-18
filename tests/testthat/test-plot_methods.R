test_that("plot.boin_decision_table builds a ggplot with one cell per decision", {
  skip_if_not_installed("ggplot2")
  expect_true(have_package("ggplot2"))

  decisions <- boin_decision_table(target = 0.30, max_n = 12)
  p <- plot(decisions)

  expect_s3_class(p, "ggplot")
  # One row per non-missing entry of the table.
  expect_equal(nrow(p$data), sum(!is.na(decisions)))
  expect_equal(levels(p$data$decision), c("E", "S", "D", "DE"))
  expect_true(all(p$data$n_tox <= p$data$n_pts))
})

test_that("the plotted decisions agree with the table", {
  skip_if_not_installed("ggplot2")

  decisions <- boin_decision_table(target = 0.30, max_n = 9)
  cells <- plot(decisions)$data

  for (i in seq_len(nrow(cells))) {
    expect_identical(
      as.character(cells$decision[i]),
      decisions[as.character(cells$n_tox[i]), as.character(cells$n_pts[i])]
    )
  }
})

test_that("a custom palette is accepted and a bad one is rejected", {
  skip_if_not_installed("ggplot2")

  decisions <- boin_decision_table(target = 0.30, max_n = 9)
  mine <- c(E = "#4DAF4A", S = "#377EB8", D = "#FF7F00", DE = "#E41A1C")

  expect_s3_class(plot(decisions, colors = mine), "ggplot")
  expect_error(plot(decisions, colors = c(E = "red", S = "blue")), "named")
  expect_error(plot(decisions, colors = 1:4), "named")
  expect_error(plot(decisions, text_size = 0), "positive number")
  expect_error(plot(decisions, text_size = c(1, 2)), "positive number")
})

test_that("the figure can be rendered without error", {
  skip_if_not_installed("ggplot2")

  decisions <- boin_decision_table(target = 0.30, max_n = 9)
  built <- ggplot2::ggplot_build(plot(decisions, text_size = 2))

  expect_s3_class(built, "ggplot_built")
  expect_gt(nrow(built$data[[1]]), 0L)
})
