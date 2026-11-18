scenarios <- list(
  "Scenario 16" = list(
    p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2),
    p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)),
  "Scenario 17" = list(
    p_true = matrix(c(0.05, 0.10, 0.15, 0.20, 0.40, 0.50), nrow = 2),
    p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.55), nrow = 2))
)

test_that("sim_tite_comb_boin12_multi matches sim_tite_comb_boin12", {
  args <- list(target = 0.33, target_resp = 0.2, n_cohort = 14,
               cohort_size = 3, window = 1, resp_window = 2,
               accrual_rate = 2, n_trials = 100, est_method = "observed",
               seed = 1)
  multi <- do.call(sim_tite_comb_boin12_multi,
                   c(args, list(scenarios = scenarios)))
  expect_s3_class(multi, "comb_boin12_oc_multi")
  expect_equal(multi$scenario_names, names(scenarios))
  single <- do.call(sim_tite_comb_boin12,
                    c(args, list(p_true = scenarios[[2]]$p_true,
                                 p_resp = scenarios[[2]]$p_resp)))
  expect_identical(multi$results[[2]]$sel_percent, single$sel_percent)
  expect_equal(multi$results[[2]]$duration, single$duration)
  tab <- multi$summary_table
  expect_equal(tab$Duration[2], round(single$duration, 1))
  expect_equal(multi$results[[1]]$design, "TITE-Comb-BOIN12")
  out <- utils::capture.output(res <- print(multi))
  expect_identical(res, multi)
  expect_true(any(grepl("Duration", out, fixed = TRUE)))
})
