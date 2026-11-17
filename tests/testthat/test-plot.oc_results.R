test_that("the results are drawn as bars or tiles", {
  skip_if_not_installed("ggplot2")
  expect_true(have_package("ggplot2"))
  p <- c(0.10, 0.25, 0.40, 0.55)
  boin <- sim_boin(target = 0.30, p_true = p, n_cohort = 6, cohort_size = 3,
                   n_trials = 50, seed = 1)
  res <- oc_results(BOIN = boin, "3+3" = oc_3p3(p))
  g <- plot(res)
  expect_s3_class(g, "ggplot")
  expect_identical(nrow(g$data), 8L)
  expect_identical(g$data$value[g$data$design == "BOIN"],
                   unname(boin$sel_percent))
  expect_identical(levels(g$data$dose), c("1", "2", "3", "4"))
  g2 <- plot(res, metric = "n_pts")
  expect_identical(g2$data$value[g2$data$design == "3+3"],
                   oc_3p3(p)$n_pts_dose)
  expect_error(plot(res, metric = "n_resp"), "no design records")
  expect_error(plot(res, metric = "total_n_pts"), "must be one of")

  m <- matrix(c(0.05, 0.10, 0.20, 0.10, 0.20, 0.30), nrow = 3)
  comb <- sim_comb_boin(target = 0.30, p_true = m, n_cohort = 4,
                        cohort_size = 3, n_trials = 20, seed = 1)
  tiles <- plot(oc_results(comb))
  expect_s3_class(tiles, "ggplot")
  expect_identical(nrow(tiles$data), 6L)
  expect_error(plot.oc_results(data.frame(a = 1)), "oc_results")
})
