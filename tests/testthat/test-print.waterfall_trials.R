test_that("print.waterfall_trials shows a compact summary", {
  trials <- waterfall_simulate(
    target = 0.30, p_true = zhang2016_p_true(1), cohort_size = 3,
    n_trials = 100, seed = 1
  )
  out <- capture.output(print(trials))

  expect_identical(out[1L], "Simulated waterfall drug-combination trials")
  expect_true(any(grepl("2 of drug A, 3 of drug B", out, fixed = TRUE)))
  expect_true(any(grepl("at most 6/3 per subtrial, of size 3", out,
                        fixed = TRUE)))
  expect_true(any(grepl("Average patients treated per trial:", out,
                        fixed = TRUE)))
  expect_true(any(grepl("Reason for ending (%):", out, fixed = TRUE)))

  result <- quiet_print(trials)
  expect_false(result$visible)
  expect_identical(result$value, trials)
})
