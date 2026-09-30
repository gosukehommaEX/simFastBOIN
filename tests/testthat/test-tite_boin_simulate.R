test_that("the trials agree with an independent implementation", {
  # Reference values from a separate Python implementation of the design,
  # run on a replica of R's Mersenne-Twister stream after set.seed() and on the
  # xoshiro256** arrival stream seeded with the same seed.
  as_mat <- function(x, n_row) matrix(as.integer(x), nrow = n_row, byrow = TRUE)

  a <- tite_boin_simulate(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 5, seed = 2026
  )
  expect_identical(unname(a$n_pts), as_mat(c(
    12, 18, 0, 0, 0, 3, 12, 12, 3, 0, 3, 9, 18, 0, 0, 3, 6, 15, 6, 0,
    6, 18, 0, 0, 0), 5))
  expect_identical(unname(a$n_tox), as_mat(c(
    1, 6, 0, 0, 0, 0, 2, 2, 2, 0, 0, 1, 6, 0, 0, 0, 1, 6, 3, 0,
    0, 5, 0, 0, 0), 5))
  expect_false(any(a$eliminated))
  expect_identical(a$cohorts_used, c(10L, 10L, 10L, 10L, 8L))
  expect_identical(a$stop_reason, c(rep("max_sample_size", 4), "n_earlystop"))
  expect_equal(a$duration, c(22.56294001594462, 22.867098284522033,
                             26.52971612305278, 30.239684588806284,
                             20.305728181527975), tolerance = 1e-10)
  expect_identical(a$n_suspensions, c(2L, 4L, 4L, 4L, 3L))
  expect_equal(a$time_suspended, c(2.2623615510612023, 7.716774943525113,
                                   6.980862588110228, 8.67892822454493,
                                   5.7854150823373764), tolerance = 1e-10)

  b <- tite_boin_simulate(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    method = "ess", accrual = "uniform", dlt_time = "uniform",
    extrasafe = TRUE, n_trials = 5, seed = 7
  )
  expect_identical(unname(b$n_pts), as_mat(c(
    6, 15, 9, 0, 0, 3, 3, 9, 12, 3, 3, 3, 15, 9, 0, 9, 3, 12, 6, 0,
    3, 6, 9, 9, 3), 5))
  expect_identical(unname(b$n_tox), as_mat(c(
    0, 3, 3, 0, 0, 0, 0, 1, 6, 1, 0, 0, 3, 5, 0, 1, 0, 2, 4, 0,
    0, 1, 1, 3, 3), 5))
  expect_identical(unname(b$eliminated), matrix(as.logical(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 1)), nrow = 5, byrow = TRUE))
  expect_identical(b$cohorts_used, rep(10L, 5))
  expect_identical(b$stop_reason, rep("max_sample_size", 5))
  expect_equal(b$duration, c(18.841098900311575, 27.631077603952534,
                             23.30592069787039, 20.976379341180554,
                             22.701111913120936), tolerance = 1e-10)
  expect_identical(b$n_suspensions, c(1L, 5L, 4L, 3L, 4L))
  expect_equal(b$time_suspended, c(1.8816213086497962, 8.407485412178046,
                                   5.188505216809293, 2.8403812716033263,
                                   6.835534050265028), tolerance = 1e-10)
})

test_that("without pending patients the trials are those of boin_simulate", {
  # Arrivals one unit apart and a window of half a unit: every patient has
  # completed the assessment before the next decision.
  settings <- list(
    list(target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
         n_cohort = 10, cohort_size = 3),
    list(target = 0.30, p_true = c(0.45, 0.55, 0.65), n_cohort = 10,
         cohort_size = 3, extrasafe = TRUE, offset = 0.10, n_earlystop = 9,
         n_earlystop_rule = "simple"),
    list(target = 0.30, p_true = c(0.02, 0.05, 0.08, 0.12), n_cohort = 7,
         cohort_size = c(1, 1, 3, 3, 3, 3, 3), start_dose = 2, n_earlystop = 9),
    list(target = 0.20, p_true = c(0.10, 0.20, 0.35, 0.50), n_cohort = 12,
         cohort_size = 3, extrasafe = TRUE),
    list(target = 0.30, p_true = c(0.25, 0.40, 0.55), n_cohort = 15,
         cohort_size = 2, p_saf = 0.20, p_tox = 0.40)
  )

  for (i in seq_along(settings)) {
    set <- settings[[i]]
    boin <- do.call(boin_simulate, c(set, list(n_trials = 500, seed = 100 + i)))
    for (method in c("imputation", "ess")) {
      tite <- do.call(tite_boin_simulate, c(set, list(
        window = 0.5, accrual_rate = 1, accrual = "fixed", method = method,
        n_trials = 500, seed = 100 + i
      )))
      info <- paste("setting", i, method)
      expect_identical(tite$n_pts, boin$n_pts, info = info)
      expect_identical(tite$n_tox, boin$n_tox, info = info)
      expect_identical(tite$eliminated, boin$eliminated, info = info)
      expect_identical(tite$cohorts_used, boin$cohorts_used, info = info)
      expect_identical(tite$stop_reason, boin$stop_reason, info = info)
      expect_true(all(tite$n_suspensions == 0L), info = info)
    }
  }
})

test_that("tite_boin_simulate returns a well formed object", {
  trials <- tite_boin_simulate(
    target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 300, seed = 1
  )

  expect_s3_class(trials, "tite_boin_trials")
  expect_equal(dim(trials$n_pts), c(300L, 5L))
  expect_equal(dim(trials$n_tox), c(300L, 5L))
  expect_equal(dim(trials$eliminated), c(300L, 5L))
  expect_equal(colnames(trials$n_pts), paste0("DL", 1:5))
  expect_type(trials$n_pts, "integer")
  expect_type(trials$eliminated, "logical")
  expect_length(trials$duration, 300L)
  expect_length(trials$n_suspensions, 300L)
  expect_length(trials$time_suspended, 300L)

  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(rowSums(trials$n_pts) <= 30))
  expect_true(all(trials$cohorts_used >= 1 & trials$cohorts_used <= 10))
  expect_true(all(trials$duration > 0))
  expect_true(all(trials$time_suspended >= 0))
  expect_true(all(trials$time_suspended[trials$n_suspensions == 0L] == 0))
  expect_true(all(trials$time_suspended < trials$duration))
  expect_true(all(trials$stop_reason %in%
                    c("lowest_dose_eliminated", "lowest_dose_too_toxic",
                      "n_earlystop", "max_sample_size", "max_cohorts")))
  expect_identical(trials$settings$method, "imputation")
  expect_equal(trials$settings$max_pending_ratio, 0.5)
  expect_identical(trials$settings$min_completed, 0L)
})

test_that("the suspension rules act as specified", {
  args <- list(target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
               n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
               n_trials = 300, seed = 5)

  # Six patients arrive per window, so accrual is often suspended by default.
  default <- do.call(tite_boin_simulate, args)
  expect_gt(mean(default$n_suspensions > 0), 0.5)

  # Without either rule accrual is never suspended.
  no_rule <- do.call(tite_boin_simulate,
                     c(args, list(max_pending_ratio = 1, min_completed = 0)))
  expect_true(all(no_rule$n_suspensions == 0L))
  expect_true(all(no_rule$time_suspended == 0))
  expect_lt(mean(no_rule$duration), mean(default$duration))

  # The effective sample size method suspends only for escalation.
  ess <- do.call(tite_boin_simulate, c(args, list(method = "ess")))
  expect_identical(ess$settings$min_completed, 2L)
  expect_gt(sum(ess$n_suspensions), 0L)
})

test_that("the same seed reproduces the same trials", {
  args <- list(target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
               n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
               n_trials = 100)

  first <- do.call(tite_boin_simulate, c(args, list(seed = 7)))
  second <- do.call(tite_boin_simulate, c(args, list(seed = 7)))
  other <- do.call(tite_boin_simulate, c(args, list(seed = 8)))

  expect_identical(first$n_pts, second$n_pts)
  expect_identical(first$duration, second$duration)
  expect_false(identical(first$duration, other$duration))
})

test_that("the state of the random number generator is restored", {
  set.seed(99)
  before <- .Random.seed
  invisible(tite_boin_simulate(
    target = 0.30, p_true = c(0.10, 0.25, 0.40), n_cohort = 5,
    cohort_size = 3, window = 3, accrual_rate = 2, n_trials = 50, seed = 1
  ))
  expect_identical(.Random.seed, before)
})

test_that("seed = NULL draws from the current stream", {
  args <- list(target = 0.30, p_true = c(0.10, 0.25, 0.40), n_cohort = 5,
               cohort_size = 3, window = 3, accrual_rate = 2, n_trials = 50,
               seed = NULL)

  set.seed(11)
  first <- do.call(tite_boin_simulate, args)
  set.seed(11)
  second <- do.call(tite_boin_simulate, args)

  expect_identical(first$n_pts, second$n_pts)
  expect_identical(first$duration, second$duration)
})

test_that("invalid arguments are rejected", {
  call_with <- function(...) {
    args <- list(target = 0.30, p_true = c(0.10, 0.25, 0.40), n_cohort = 5,
                 cohort_size = 3, window = 3, accrual_rate = 2, n_trials = 10)
    new <- list(...)
    args[names(new)] <- new
    do.call(tite_boin_simulate, args)
  }

  expect_error(call_with(window = 0), "window")
  expect_error(call_with(window = c(1, 2)), "window")
  expect_error(call_with(accrual_rate = -1), "accrual_rate")
  expect_error(call_with(late_fraction = 1), "late_fraction")
  expect_error(call_with(p_true = c(0.5, 1)), "below 1")
  expect_s3_class(call_with(p_true = c(0.5, 1), dlt_time = "uniform"),
                  "tite_boin_trials")
  expect_error(call_with(seed = 1.5), "seed")
  expect_error(call_with(start_dose = 4), "start_dose")
  expect_error(call_with(method = "crm"), "should be one of")
  expect_error(call_with(accrual = "poisson"), "should be one of")
  expect_error(call_with(max_pending_ratio = 0), "max_pending_ratio")
  expect_warning(call_with(n_earlystop = 6), "n_earlystop")
})
