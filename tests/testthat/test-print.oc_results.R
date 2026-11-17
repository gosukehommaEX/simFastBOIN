test_that("the results print one table per scenario", {
  p <- c(0.10, 0.25, 0.40, 0.55)
  boin <- sim_boin(target = 0.30, p_true = p, n_cohort = 6, cohort_size = 3,
                   n_trials = 50, seed = 1)
  res <- oc_results(BOIN = boin, "3+3" = oc_3p3(p))
  expect_output(print(res), "2 design\\(s\\) in 1 scenario\\(s\\)")
  expect_output(print(res), "Selection \\(%\\) by dose")
  expect_output(print(res), "Patients by dose")
  expect_output(print(res, metric = "n_tox"), "DLTs by dose")
  expect_output(print(res), "No selection \\(%\\)")
  out <- capture.output(print(res, digits = 2))
  expect_true(any(grepl(formatC(boin$sel_percent[2], format = "f",
                                digits = 2), out, fixed = TRUE)))
  r <- quiet_print(res)
  expect_identical(r$value, res)
  expect_false(r$visible)

  m <- matrix(c(0.05, 0.10, 0.20, 0.10, 0.20, 0.30), nrow = 3)
  comb <- sim_comb_boin(target = 0.30, p_true = m, n_cohort = 4,
                        cohort_size = 3, n_trials = 20, seed = 1)
  expect_output(print(oc_results(comb)),
                "Selection \\(%\\) by combination, BOIN combination")
  expect_error(print(res, metric = "duration"), "per-dose quantities")
  expect_error(print.oc_results(data.frame(a = 1)), "oc_results")
})
