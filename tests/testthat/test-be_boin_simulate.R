test_that("without any response the trials are those of tite_boin_simulate", {
  # No dose is ever opened for backfilling, so the engine must reproduce
  # TITE-BOIN with the same suspension rules trial by trial, timing included.
  p_true <- c(0.04, 0.12, 0.25, 0.43, 0.63)
  outputs <- c("n_pts", "n_tox", "eliminated", "cohorts_used", "stop_reason",
               "duration", "n_suspensions", "time_suspended")
  settings <- list(
    list(max_pending_ratio = 0.49, min_follow_up = 0.25, extrasafe = FALSE,
         accrual = "exponential", n_earlystop = 9),
    list(max_pending_ratio = 0.5, min_follow_up = 0, extrasafe = TRUE,
         accrual = "exponential", n_earlystop = 18),
    list(max_pending_ratio = 0.49, min_follow_up = 0.25, extrasafe = FALSE,
         accrual = "uniform", n_earlystop = 12)
  )

  for (set in settings) {
    common <- list(target = 0.25, p_true = p_true, n_cohort = 10,
                   cohort_size = 3, window = 3, accrual_rate = 2,
                   accrual = set$accrual,
                   max_pending_ratio = set$max_pending_ratio,
                   min_follow_up = set$min_follow_up,
                   extrasafe = set$extrasafe, n_earlystop = set$n_earlystop,
                   n_trials = 500, seed = 31)
    be <- do.call(be_boin_simulate, c(common, list(p_resp = rep(0, 5))))
    tite <- do.call(tite_boin_simulate, c(common, list(method = "imputation")))
    for (component in outputs) {
      expect_identical(unname(be[[component]]), unname(tite[[component]]),
                       info = component)
    }
    expect_true(all(be$n_bf == 0L))
  }
})

test_that("backfilling follows its rules", {
  trials <- be_boin_simulate(
    target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
    p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58), n_cohort = 10, cohort_size = 3,
    window = 3, accrual_rate = 2, n_earlystop = 9, n_trials = 1000, seed = 2
  )

  expect_gt(mean(rowSums(trials$n_bf) > 0), 0.5)
  expect_true(all(trials$n_bf <= trials$n_pts))
  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(trials$n_bf[, 5] == 0L))
  expect_true(all(trials$n_bf <= 12L))
  expect_true(all(rowSums(trials$n_pts - trials$n_bf) <= 30L))
})

test_that("patients can be turned away instead of waiting", {
  args <- list(target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
               p_resp = rep(0, 5), n_cohort = 10, cohort_size = 3,
               window = 3, accrual_rate = 2, n_trials = 300, seed = 5)

  wait <- do.call(be_boin_simulate, args)
  leave <- do.call(be_boin_simulate, c(args, list(no_slot = "leave")))
  expect_true(all(wait$n_turned_away == 0L))
  expect_gt(mean(leave$n_turned_away), 0)
})

test_that("the defaults are those of the design", {
  trials <- be_boin_simulate(
    target = 0.25, p_true = c(0.10, 0.25, 0.40), p_resp = c(0.2, 0.3, 0.4),
    n_cohort = 5, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 100, seed = 6
  )

  expect_s3_class(trials, "be_boin_trials")
  expect_s3_class(trials, "backfill_trials")
  expect_identical(trials$settings$design, "be_boin")
  expect_identical(trials$settings$conflict_dose, "lowest")
  expect_identical(trials$settings$backfill_dose, "highest")
  expect_equal(trials$settings$max_pending_ratio, 0.49)
  expect_equal(trials$settings$min_follow_up, 0.25)
  expect_identical(trials$settings$min_completed, 0L)
  expect_null(trials$settings$stay_on_1_of_3)
})

test_that("invalid arguments are rejected", {
  call_with <- function(...) {
    args <- list(target = 0.25, p_true = c(0.10, 0.25, 0.40),
                 p_resp = c(0.2, 0.3, 0.4), n_cohort = 5, cohort_size = 3,
                 window = 3, accrual_rate = 2, n_trials = 10)
    new <- list(...)
    args[names(new)] <- new
    do.call(be_boin_simulate, args)
  }

  expect_error(call_with(max_pending_ratio = 0), "max_pending_ratio")
  expect_error(call_with(min_follow_up = 2), "min_follow_up")
  expect_error(call_with(min_completed = -1), "min_completed")
  expect_error(call_with(conflict_dose = "middle"), "should be one of")
  expect_error(call_with(p_resp = c(0.2, 0.3)), "p_resp")
})
