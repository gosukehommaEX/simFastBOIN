test_that("the comparison is drawn as points by scenario", {
  skip_if_not_installed("ggplot2")
  expect_true(have_package("ggplot2"))
  p <- c(0.10, 0.25, 0.40, 0.55)
  truth <- oc_scenarios(list(A = p, B = p / 2), target = 0.30)
  boin <- sim_boin_multi(target = 0.30, scenarios = list(A = p, B = p / 2),
                         n_cohort = 6, cohort_size = 3, n_trials = 50,
                         seed = 1)
  cmp <- oc_compare(BOIN = boin,
                    "3+3" = list(A = oc_3p3(p), B = oc_3p3(p / 2)),
                    truth = truth)
  g <- plot(cmp)
  expect_s3_class(g, "ggplot")
  # Four quantities for two designs in two scenarios.
  expect_identical(nrow(g$data), 16L)
  expect_identical(levels(g$data$metric)[1], "Correct selection (%)")
  # Lines join the two scenarios of every design.
  expect_identical(length(g$layers), 2L)
  one <- plot(cmp[cmp$scenario == "A", ], metrics = "pct_correct")
  expect_identical(length(one$layers), 1L)
  expect_error(plot(cmp, metrics = "x"), "must name columns")
  expect_error(plot.oc_compare(data.frame(a = 1)), "oc_compare")
})
