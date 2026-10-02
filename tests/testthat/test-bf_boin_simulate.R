test_that("without any response the trials are those of boin_simulate", {
  # No dose is ever opened for backfilling, so the engine must reproduce the
  # BOIN design trial by trial with the same seed.
  p_true <- c(0.04, 0.12, 0.25, 0.43, 0.63)
  outputs <- c("n_pts", "n_tox", "eliminated", "cohorts_used", "stop_reason")
  settings <- list(
    list(stay_on_1_of_3 = FALSE, extrasafe = FALSE, n_earlystop = 9,
         start_dose = 1),
    list(stay_on_1_of_3 = TRUE, extrasafe = FALSE, n_earlystop = 9,
         start_dose = 1),
    list(stay_on_1_of_3 = FALSE, extrasafe = TRUE, n_earlystop = 12,
         start_dose = 2)
  )

  for (set in settings) {
    bf <- bf_boin_simulate(
      target = 0.25, p_true = p_true, p_resp = rep(0, 5), n_cohort = 10,
      cohort_size = 3, window = 1, accrual_rate = 3, n_trials = 500,
      start_dose = set$start_dose, n_earlystop = set$n_earlystop,
      extrasafe = set$extrasafe, stay_on_1_of_3 = set$stay_on_1_of_3,
      seed = 21
    )
    boin <- boin_simulate(
      target = 0.25, p_true = p_true, n_cohort = 10, cohort_size = 3,
      n_trials = 500, start_dose = set$start_dose,
      n_earlystop = set$n_earlystop, extrasafe = set$extrasafe,
      stay_on_1_of_3 = set$stay_on_1_of_3, seed = 21
    )
    for (component in outputs) {
      expect_identical(unname(bf[[component]]), unname(boin[[component]]),
                       info = component)
    }
    expect_true(all(bf$n_bf == 0L))
    expect_true(all(bf$n_resp == 0L))
  }
})

test_that("backfilling follows its rules", {
  # Scenario 3 of Zhao et al. (2024).
  trials <- bf_boin_simulate(
    target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
    p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58), n_cohort = 10, cohort_size = 3,
    window = 1, accrual_rate = 3, n_earlystop = 9, stay_on_1_of_3 = TRUE,
    n_trials = 1000, seed = 2
  )

  expect_gt(mean(rowSums(trials$n_bf) > 0), 0.5)
  expect_true(all(trials$n_bf <= trials$n_pts))
  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(trials$n_resp <= trials$n_pts))
  # Backfilling is below the current dose, so never at the highest dose, and
  # stops once a dose has n_cap patients.
  expect_true(all(trials$n_bf[, 5] == 0L))
  expect_true(all(trials$n_bf <= 12L))
  # The escalation itself keeps its size.
  expect_true(all(rowSums(trials$n_pts - trials$n_bf) <= 30L))
  expect_true(all(trials$n_turned_away == 0L))
})

test_that("a response is needed at or below a backfilled dose", {
  args <- list(target = 0.25, p_true = c(0.02, 0.05, 0.08, 0.11, 0.25),
               n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
               n_trials = 300, seed = 3)

  # Responses only at the highest dose: nothing can be backfilled.
  top <- do.call(bf_boin_simulate, c(args, list(p_resp = c(0, 0, 0, 0, 0.5))))
  expect_true(all(top$n_bf == 0L))
  expect_gt(sum(top$n_resp), 0L)

  # Responses only at dose 2: doses 2 to 4 can be backfilled, dose 1 cannot.
  second <- do.call(bf_boin_simulate, c(args, list(p_resp = c(0, 0.5, 0, 0, 0))))
  expect_true(all(second$n_bf[, 1] == 0L))
  expect_gt(sum(second$n_bf[, 2:4]), 0L)
})

test_that("n_cap and backfill_dose change where patients are backfilled", {
  args <- list(target = 0.25, p_true = c(0.02, 0.05, 0.08, 0.11, 0.25),
               p_resp = rep(0.5, 5), n_cohort = 10, cohort_size = 3,
               window = 1, accrual_rate = 3, n_trials = 300, seed = 4)

  capped <- do.call(bf_boin_simulate, c(args, list(n_cap = 4)))
  expect_true(all(capped$n_bf <= 4L))

  highest <- do.call(bf_boin_simulate, args)
  lowest <- do.call(bf_boin_simulate, c(args, list(backfill_dose = "lowest")))
  expect_gt(mean(lowest$n_bf[, 1]), mean(highest$n_bf[, 1]))
  expect_identical(lowest$settings$backfill_dose, "lowest")
})

test_that("patients can be turned away instead of waiting", {
  args <- list(target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
               p_resp = rep(0, 5), n_cohort = 10, cohort_size = 3,
               window = 1, accrual_rate = 3, n_trials = 300, seed = 5)

  wait <- do.call(bf_boin_simulate, args)
  leave <- do.call(bf_boin_simulate, c(args, list(no_slot = "leave")))
  expect_true(all(wait$n_turned_away == 0L))
  expect_gt(mean(leave$n_turned_away), 0)
  # Turning patients away changes only the timing, not the decisions.
  expect_identical(leave$n_pts, wait$n_pts)
  expect_gt(mean(leave$duration), mean(wait$duration))
})

test_that("the result has the expected shape and settings", {
  trials <- bf_boin_simulate(
    target = 0.25, p_true = c(0.10, 0.25, 0.40), p_resp = c(0.2, 0.3, 0.4),
    n_cohort = 5, cohort_size = 3, window = 1, accrual_rate = 3,
    n_trials = 100, seed = 6
  )

  expect_s3_class(trials, "bf_boin_trials")
  expect_s3_class(trials, "backfill_trials")
  for (component in c("n_pts", "n_tox", "n_bf", "n_resp", "eliminated")) {
    expect_identical(dim(trials[[component]]), c(100L, 3L), info = component)
    expect_identical(colnames(trials[[component]]), c("DL1", "DL2", "DL3"))
  }
  expect_length(trials$duration, 100L)
  expect_identical(trials$settings$design, "bf_boin")
  expect_identical(trials$settings$conflict_dose, "highest")
  expect_identical(trials$settings$no_slot, "wait")
  expect_equal(trials$settings$resp_window, 1)
  expect_false(trials$settings$stay_on_1_of_3)
  expect_null(trials$settings$min_follow_up)
})

test_that("the same seed reproduces the trials and the stream is restored", {
  args <- list(target = 0.25, p_true = c(0.10, 0.25, 0.40),
               p_resp = c(0.2, 0.3, 0.4), n_cohort = 5, cohort_size = 3,
               window = 1, accrual_rate = 3, n_trials = 100)

  first <- do.call(bf_boin_simulate, c(args, list(seed = 7)))
  second <- do.call(bf_boin_simulate, c(args, list(seed = 7)))
  expect_identical(first$n_bf, second$n_bf)
  expect_identical(first$duration, second$duration)

  set.seed(99)
  before <- .Random.seed
  invisible(do.call(bf_boin_simulate, c(args, list(seed = 1))))
  expect_identical(.Random.seed, before)
})

test_that("invalid arguments are rejected", {
  call_with <- function(...) {
    args <- list(target = 0.25, p_true = c(0.10, 0.25, 0.40),
                 p_resp = c(0.2, 0.3, 0.4), n_cohort = 5, cohort_size = 3,
                 window = 1, accrual_rate = 3, n_trials = 10)
    new <- list(...)
    args[names(new)] <- new
    do.call(bf_boin_simulate, args)
  }

  expect_error(call_with(p_resp = c(0.2, 0.3)), "p_resp")
  expect_error(call_with(p_resp = c(0.2, 0.3, 1.2)), "p_resp")
  expect_error(call_with(p_resp = c(0.2, 0.3, 1)), "below 1")
  expect_s3_class(call_with(p_resp = c(0.2, 0.3, 1), dlt_time = "uniform"),
                  "bf_boin_trials")
  expect_error(call_with(n_cap = 0), "n_cap")
  expect_error(call_with(resp_window = 0), "resp_window")
  expect_error(call_with(resp_late_fraction = 1), "resp_late_fraction")
  expect_error(call_with(resp_cor = 1), "resp_cor")
  expect_error(call_with(backfill_dose = "middle"), "should be one of")
  expect_error(call_with(no_slot = "queue"), "should be one of")
  expect_error(call_with(stay_on_1_of_3 = NA), "stay_on_1_of_3")
})

test_that("the trials agree with a separate implementation", {
  # Reference values from a separate implementation in Python, written from
  # the documented rules and run on replicas of R's Mersenne-Twister and of
  # xoshiro256**. It agreed with this function on every one of 2,700 trials
  # in nine settings of both designs; these are 300 of them.
  trials <- bf_boin_simulate(
    target = 0.25, p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
    p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58), n_cohort = 10, cohort_size = 3,
    window = 1, accrual_rate = 3, n_earlystop = 9, stay_on_1_of_3 = TRUE,
    n_trials = 300, seed = 11
  )

  expect_identical(unname(colSums(trials$n_pts)), c(1812, 3313, 2978, 1082, 123))
  expect_identical(unname(colSums(trials$n_tox)), c(70, 396, 744, 478, 82))
  expect_identical(unname(colSums(trials$n_bf)), c(531, 1072, 617, 89, 0))
  expect_identical(unname(colSums(trials$n_resp)), c(182, 691, 936, 474, 76))
  expect_identical(sum(trials$n_suspensions), 2104L)
  expect_identical(sum(trials$stop_reason == "n_earlystop"), 216L)
  expect_equal(sum(trials$duration), 4062.0022202492, tolerance = 1e-12)
  expect_identical(unname(trials$n_pts[2, ]), c(3L, 3L, 14L, 6L, 0L))
  expect_identical(unname(trials$n_bf[2, ]), c(0L, 0L, 5L, 0L, 0L))
  expect_equal(trials$duration[1:3], c(5.4004614590, 14.3393459536, 7.5675044043),
               tolerance = 1e-9)
})
