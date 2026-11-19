test_that("print.boin_etc_oc_multi shows one row per scenario", {
  scenarios <- list(
    A = list(p_true = matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), 2),
             p_resp = matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), 2)),
    B = list(p_true = matrix(c(0.05, 0.15, 0.15, 0.35, 0.30, 0.55), 2),
             p_resp = matrix(c(0.25, 0.60, 0.60, 0.65, 0.65, 0.70), 2))
  )
  res <- sim_boin_etc_multi(target = 0.30, target_resp = 0.60,
                            scenarios = scenarios, n_cohort = 8,
                            cohort_size = 3, n_trials = 50,
                            subtrials = "etc2")
  expect_output(print(res), "BOIN-ETC2 operating characteristics across 2")
  expect_output(print(res), "rule 4               : paper")
  out <- capture.output(print(res))
  expect_true(any(grepl("^ +A ", out)))
  expect_true(any(grepl("^ +B ", out)))
})
