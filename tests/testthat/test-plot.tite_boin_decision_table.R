test_that("plot.tite_boin_decision_table draws one cell per state", {
  skip_if_not_installed("ggplot2")
  expect_true(have_package("ggplot2"))

  tab <- tite_boin_decision_table(target = 0.2, max_n = 15)
  p <- plot(tab, n = c(9, 12, 15))

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), sum(tab$n %in% c(9, 12, 15)))
  expect_equal(levels(p$data$fill_group),
               c("E", "S", "D", "DE", "SUS", "Depends"))
  expect_equal(levels(p$data$panel),
               paste(c(9, 12, 15), "patients treated"))
  expect_true(all(p$data$n_tox + p$data$n_pending <=
                    as.integer(sub(" .*$", "", as.character(p$data$panel)))))
})

test_that("the labels carry the decisions and the boundaries", {
  skip_if_not_installed("ggplot2")

  tab <- tite_boin_decision_table(target = 0.2, max_n = 9)
  cells <- plot(tab, n = 9)$data
  nine <- tab[tab$n == 9, ]

  expect_identical(sub("\n.*$", "", cells$label), nine$decision)
  k <- which(cells$n_tox == 1 & cells$n_pending == 4)
  expect_identical(cells$label[k], "E/S\n>=2.15")
  k <- which(cells$n_tox == 2 & cells$n_pending == 1)
  expect_identical(cells$label[k], "S/D\n<=0.52")
  expect_identical(as.character(cells$fill_group[k]), "Depends")

  ess <- tite_boin_decision_table(target = 0.3, max_n = 6, method = "ess")
  cells <- plot(ess, n = 6)$data
  k <- which(cells$n_tox == 1 & cells$n_pending == 5)
  expect_identical(cells$label[k], "SUS/S/D\n>=4.23\n<=2.79")
})

test_that("escalations that require a minimum follow-up are marked", {
  skip_if_not_installed("ggplot2")

  tab <- tite_boin_decision_table(target = 0.25, max_n = 6,
                                  max_pending_ratio = 0.49,
                                  min_follow_up = 0.25)
  p <- plot(tab, n = 6)
  cells <- p$data
  k <- which(cells$n_tox == 1 & cells$n_pending == 1)
  expect_identical(cells$label[k], "E/S*\n>=0.22")
  k <- which(cells$n_tox == 0 & cells$n_pending == 2)
  expect_identical(cells$label[k], "E*")
  k <- which(cells$n_tox == 0 & cells$n_pending == 0)
  expect_identical(cells$label[k], "E")
  expect_match(p$labels$caption, "at least 25% of the assessment window",
               fixed = TRUE)

  plain <- plot(tite_boin_decision_table(target = 0.25, max_n = 6), n = 6)
  expect_false(any(grepl("*", plain$data$label, fixed = TRUE)))
  expect_null(plain$labels$caption)
})

test_that("cohort_size selects the panels", {
  skip_if_not_installed("ggplot2")

  tab <- tite_boin_decision_table(target = 0.3, max_n = 12)
  p <- plot(tab, cohort_size = 3)

  expect_equal(levels(p$data$panel), paste(c(3, 6, 9, 12), "patients treated"))
  expect_error(plot(tab, n = 4, cohort_size = 3), "no number of patients")
})

test_that("invalid arguments are rejected", {
  skip_if_not_installed("ggplot2")

  tab <- tite_boin_decision_table(target = 0.3, max_n = 9)
  mine <- c(E = "#4DAF4A", S = "#377EB8", D = "#FF7F00", DE = "#E41A1C",
            SUS = "#999999", Depends = "#DDDDDD")

  expect_s3_class(plot(tab, n = 9, colors = mine), "ggplot")
  expect_error(plot(tab, colors = mine[1:4]), "named")
  expect_error(plot(tab, n = 10), "not in the table")
  expect_error(plot(tab, n = 2.5), "whole numbers")
  expect_error(plot(tab, text_size = 0), "positive number")
  expect_error(plot(tab, digits = -1), "digits")
  expect_error(plot(tab[, c("n", "decision")]), "tite_boin_decision_table")
})

test_that("the figure can be rendered without error", {
  skip_if_not_installed("ggplot2")

  tab <- tite_boin_decision_table(target = 0.3, max_n = 9)
  built <- ggplot2::ggplot_build(plot(tab, cohort_size = 3, text_size = 1.5))

  expect_s3_class(built, "ggplot_built")
  expect_gt(nrow(built$data[[1]]), 0L)
})
