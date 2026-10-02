test_that("print.backfill_trials shows the design and the backfilling", {
  args <- list(target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
               p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58), n_cohort = 10,
               cohort_size = 3, accrual_rate = 3, n_trials = 100, seed = 1)

  bf <- do.call(bf_boin_simulate, c(args, list(window = 1)))
  out <- capture.output(print(bf))
  expect_match(out[1L], "BF-BOIN", fixed = TRUE)
  expect_true(any(grepl("of which backfilled", out, fixed = TRUE)))
  expect_true(any(grepl("Trials with backfill", out, fixed = TRUE)))
  expect_false(any(grepl("turned away", out, fixed = TRUE)))

  be <- do.call(be_boin_simulate, c(args, list(window = 3, no_slot = "leave")))
  out <- capture.output(print(be))
  expect_match(out[1L], "BE-BOIN", fixed = TRUE)
  expect_true(any(grepl("Patients turned away", out, fixed = TRUE)))
})

test_that("print.backfill_trials returns its input invisibly", {
  trials <- bf_boin_simulate(
    target = 0.25, p_true = c(0.10, 0.25, 0.40), p_resp = c(0.2, 0.3, 0.4),
    n_cohort = 5, cohort_size = 3, window = 1, accrual_rate = 3,
    n_trials = 50, seed = 2
  )
  result <- quiet_print(trials)

  expect_false(result$visible)
  expect_identical(result$value, trials)
})
