test_that("simulated trials print their summary", {
  trials <- boin12_simulate(target = 0.35, target_resp = 0.25,
                            p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                            p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
                            n_cohort = 12, cohort_size = 3, n_trials = 50,
                            seed = 1)
  expect_output(print(trials), "Simulated BOIN12 trials")
  expect_output(print(trials), "Stopping reason")
  out <- capture.output(res <- print(trials))
  expect_identical(res, trials)
})
