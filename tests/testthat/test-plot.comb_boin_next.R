zero_matrix <- function() matrix(0, nrow = 3, ncol = 4)

test_that("plot.comb_boin_next draws one cell per combination with its status", {
  skip_if_not_installed("ggplot2")
  expect_true(have_package("ggplot2"))

  n_pts <- rbind(c(6, 6, 0, 0), c(3, 3, 0, 0), c(0, 0, 0, 0))
  n_tox <- rbind(c(0, 1, 0, 0), c(1, 2, 0, 0), c(0, 0, 0, 0))
  nxt <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                        current = c(2, 2))
  p <- plot(nxt)

  expect_s3_class(p, "ggplot")
  cells <- p$data
  expect_equal(nrow(cells), 12L)
  expect_equal(levels(cells$status),
               c("Next", "Current", "Eliminated", "Treated", "Not treated"))
  at <- function(a, b) cells$dose_a == a & cells$dose_b == b
  expect_identical(as.character(cells$status[at(1, 2)]), "Next")
  expect_identical(as.character(cells$status[at(2, 2)]), "Current")
  expect_identical(as.character(cells$status[at(1, 1)]), "Treated")
  expect_identical(as.character(cells$status[at(3, 4)]), "Not treated")
  expect_identical(cells$label[at(1, 2)], "1/6\nnext\nPr = 0.22")
  expect_identical(cells$label[at(2, 1)], "1/3\nPr = 0.17")
  expect_identical(cells$label[at(2, 2)], "2/3\ncurrent")
  expect_identical(cells$label[at(3, 1)], "-")
  expect_identical(p$labels$title, "De-escalate to A1B2")
})

test_that("ties, stays and stops are drawn", {
  skip_if_not_installed("ggplot2")

  n_pts <- zero_matrix()
  n_pts[1, 1] <- 3
  tie <- plot(comb_boin_next(target = 0.25, n_pts = n_pts,
                             n_tox = zero_matrix(), current = c(1, 1)))
  expect_equal(sum(tie$data$status == "Next"), 2L)
  expect_identical(tie$labels$title, "Escalate to A2B1 or A1B2")

  n_pts <- rbind(c(3, 3, 0, 0), c(3, 9, 0, 0), c(0, 0, 0, 0))
  n_tox <- rbind(c(0, 0, 0, 0), c(0, 2, 0, 0), c(0, 0, 0, 0))
  stay <- plot(comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                              current = c(2, 2)))
  expect_equal(sum(stay$data$status == "Next"), 0L)
  expect_identical(stay$labels$title, "Stay at A2B2")
  expect_true(any(stay$data$label == "2/9\ncurrent, stay"))

  n_tox <- zero_matrix()
  n_tox[1, 1] <- 3
  n_pts <- zero_matrix()
  n_pts[1, 1] <- 3
  stopped <- plot(comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                                 current = c(1, 1)))
  expect_identical(stopped$labels$title, "Stop the trial")
  expect_equal(sum(stopped$data$status == "Eliminated"), 11L)
  expect_true(any(stopped$data$label == "3/3\ncurrent, eliminated"))
})

test_that("a custom palette is accepted and invalid arguments are rejected", {
  skip_if_not_installed("ggplot2")

  n_pts <- zero_matrix()
  n_pts[1, 1] <- 3
  nxt <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = zero_matrix(),
                        current = c(1, 1))
  mine <- c(Next = "#4DAF4A", Current = "#377EB8", Eliminated = "#E41A1C",
            Treated = "grey80", "Not treated" = "white")

  expect_s3_class(plot(nxt, colors = mine), "ggplot")
  expect_error(plot(nxt, colors = mine[1:3]), "named")
  expect_error(plot(nxt, text_size = 0), "positive number")
})

test_that("the figure can be rendered without error", {
  skip_if_not_installed("ggplot2")

  n_pts <- rbind(c(6, 6, 0, 0), c(3, 3, 0, 0), c(0, 0, 0, 0))
  n_tox <- rbind(c(0, 1, 0, 0), c(1, 2, 0, 0), c(0, 0, 0, 0))
  nxt <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                        current = c(2, 2))
  built <- ggplot2::ggplot_build(plot(nxt, text_size = 2))

  expect_s3_class(built, "ggplot_built")
  expect_equal(nrow(built$data[[1]]), 12L)
})
