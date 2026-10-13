test_that("print.comb_boin_trials shows a compact summary", {
  trials <- comb_boin_simulate(
    target = 0.30,
    p_true = matrix(c(0.05, 0.15, 0.30,
                      0.15, 0.30, 0.45), nrow = 2, byrow = TRUE),
    n_cohort = 10, cohort_size = 3, n_trials = 100, seed = 1
  )
  out <- capture.output(print(trials))

  expect_identical(out[1L], "Simulated BOIN drug-combination trials")
  expect_true(any(grepl("2 of drug A, 3 of drug B", out, fixed = TRUE)))
  expect_true(any(grepl("Average patients treated per trial:", out,
                        fixed = TRUE)))
  expect_true(any(grepl("Stopping reason (%):", out, fixed = TRUE)))

  result <- quiet_print(trials)
  expect_false(result$visible)
  expect_identical(result$value, trials)
})
