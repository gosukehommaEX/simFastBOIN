test_that("print.backfill_oc_multi shows one table for every scenario", {
  scenarios <- list(
    list(name = "Rising", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
         p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58)),
    list(name = "Plateau", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
         p_resp = c(0.30, 0.32, 0.35, 0.36, 0.36))
  )
  oc <- sim_bf_boin_multi(
    target = 0.25, scenarios = scenarios, n_cohort = 10, cohort_size = 3,
    window = 1, accrual_rate = 3, n_trials = 100, seed = 1
  )
  out <- capture.output(print(oc))

  expect_match(out[1L], "BF-BOIN operating characteristics across 2 scenarios",
               fixed = TRUE)
  expect_equal(sum(grepl("Backfilled patients", out, fixed = TRUE)), 2L)
  expect_equal(sum(grepl("Trial duration (mean)", out, fixed = TRUE)), 2L)

  pct <- capture.output(print(oc, percent = TRUE))
  expect_equal(sum(grepl("Backfilled patients (%)", pct, fixed = TRUE)), 2L)

  result <- quiet_print(oc)
  expect_false(result$visible)
  expect_identical(result$value, oc)
  expect_error(print(oc, percent = "yes"), "percent")
})
