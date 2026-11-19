test_that("print.boin_etc_trials shows the averages and the stopping reasons", {
  p_true <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
  p_resp <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
  trials <- boin_etc_simulate(target = 0.30, target_resp = 0.60,
                              p_true = p_true, p_resp = p_resp, n_cohort = 8,
                              cohort_size = 3, n_trials = 50,
                              subtrials = "etc3", seed = 1)
  expect_output(print(trials), "Simulated BOIN-ETC3 trials")
  expect_output(print(trials), "Patients treated, average per trial")
  expect_output(print(trials), "Stopping reason")
})
