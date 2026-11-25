test_that("print.uboin_trials shows a compact summary", {
  sc <- zhou2019_scenario(1)
  trials <- uboin_simulate(tox_upper = 0.30, resp_lower = 0.20,
                           p_true = sc$p_true, p_resp = sc$p_resp,
                           n_cohort = 18, cohort_size = 3, n_trials = 100,
                           seed = 1)
  out <- capture.output(print(trials))

  expect_identical(out[1L], "Simulated U-BOIN trials")
  expect_true(any(grepl("doses                : 5", out, fixed = TRUE)))
  expect_true(any(grepl("Averages per trial:", out, fixed = TRUE)))
  expect_true(any(grepl("Responses observed", out, fixed = TRUE)))
  expect_true(any(grepl("Reason for ending (%):", out, fixed = TRUE)))

  result <- quiet_print(trials)
  expect_false(result$visible)
  expect_identical(result$value, trials)
})
