test_that("print.comb_boin12_trials shows the averages and the reasons", {
  trials <- comb_boin12_simulate(
    target = 0.33, target_resp = 0.2,
    p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), 2),
    p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), 2),
    n_cohort = 14, cohort_size = 3, n_trials = 50, seed = 1)
  out <- capture.output(res <- print(trials))
  expect_identical(res, trials)
  expect_true(any(grepl("Simulated Comb-BOIN12 trials", out)))
  expect_true(any(grepl("2 of drug A, 3 of drug B", out)))
  expect_true(any(grepl("Patients treated, average per trial", out)))
  expect_true(any(grepl("Stopping reason", out)))
})
