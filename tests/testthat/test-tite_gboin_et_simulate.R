test_that("tite_gboin_et_simulate returns the trials of the rules", {
  trials <- tite_takeda_trials(n_trials = 4, seed = 3)

  # Expected trials from the C++ core run outside of R on R's uniform stream
  # for this seed, and from an implementation of the rules in Python
  n_pts <- rbind(c(3, 3, 3, 3, 3, 21), c(3, 3, 3, 3, 3, 21),
                 c(3, 3, 3, 3, 3, 21), c(3, 3, 3, 6, 9, 12))
  ets_total <- rbind(c(0, 0.5, 1, 1, 1, 10), c(0, 0, 0, 1, 1, 6.5),
                     c(0, 0.5, 0, 1, 0.5, 8.5), c(0, 0, 1, 2, 3.5, 2.5))
  ees_total <- rbind(c(0, 3.25, 0.75, 1, 0.5, 39.5),
                     c(0.5, 0.25, 0.5, 1.25, 1.25, 21.5),
                     c(1.25, 3.25, 1.25, 0.25, 0.75, 37.25),
                     c(0.75, 0.25, 0.25, 2.5, 7, 12.75))
  grades_trial_1 <- rbind(c(3, 0, 0, 0), c(2, 1, 0, 0), c(2, 0, 1, 0),
                          c(2, 0, 1, 0), c(2, 0, 1, 0), c(9, 5, 6, 1))
  responses_trial_1 <- rbind(c(3, 0, 0, 0), c(1, 1, 0, 1), c(0, 3, 0, 0),
                             c(2, 0, 1, 0), c(1, 2, 0, 0), c(0, 2, 9, 10))
  duration <- c(530.630772, 469.797048, 583.925200, 525.159290)

  expect_s3_class(trials, "boin_et_trials")
  expect_identical(trials$design, "TITE-gBOIN-ET")
  expect_equal(unname(trials$n_pts), n_pts)
  expect_equal(unname(trials$n_tox) * 1.5, ets_total, tolerance = 1e-12)
  expect_equal(unname(trials$n_resp) * 3, ees_total, tolerance = 1e-12)
  expect_equal(unname(trials$n_grade[1, , ]), grades_trial_1)
  expect_equal(unname(trials$n_resp_grade[1, , ]), responses_trial_1)
  expect_lt(max(abs(trials$duration - duration)), 1e-5)
  expect_identical(trials$stop_reason, rep("max_cohorts", 4))
  expect_identical(trials$early_stop, rep(FALSE, 4))
  expect_identical(trials$cohorts_used, rep(12L, 4))
  expect_identical(dimnames(trials$n_grade)[[2]], paste0("DL", 1:6))
  expect_identical(dimnames(trials$n_grade)[[3]], paste0("G", 1:4))
  expect_identical(dimnames(trials$n_resp_grade)[[3]], paste0("R", 1:4))
})

test_that("the settings hold the true scores, the boundaries and the timing", {
  trials <- tite_takeda_trials(n_trials = 10)

  # Mean normalized scores of the doses, by hand: sum of p times score,
  # divided by the largest score
  expect_equal(trials$settings$p_true,
               c(0.035, 0.08, 0.14, 0.23, 0.295, 0.425) / 1.5,
               tolerance = 1e-12)
  expect_equal(trials$settings$p_resp,
               c(0.155, 0.24, 0.35, 0.55, 0.90, 1.7375) / 3,
               tolerance = 1e-12)
  ref <- boin_et_boundary(0.47 / 1.5, 1.75 / 3)
  expect_identical(trials$boundary$lambda1, ref$lambda1)
  expect_identical(trials$boundary$lambda2, ref$lambda2)
  expect_identical(trials$boundary$eta1, ref$eta1)
  # The boundaries printed by Takeda et al. (2023): 0.14, 0.37 and 0.46
  expect_equal(round(c(ref$lambda1, ref$lambda2, ref$eta1), 2),
               c(0.14, 0.37, 0.46))
  expect_equal(trials$settings$score_norm, c(0, 1 / 3, 2 / 3, 1))
  expect_equal(trials$settings$score_resp_norm, c(0, 0.25 / 3, 1 / 3, 1))
  expect_identical(trials$settings$window, 30)
  expect_identical(trials$settings$resp_window, 60)
  expect_identical(trials$settings$accrual, "uniform")
  expect_identical(trials$settings$event_time, "weibull")
  expect_identical(trials$settings$min_completed_share, 0.5)
  expect_identical(trials$settings$p_eli, 0.47 / 1.5)
  expect_identical(trials$settings$p_eli_resp, ref$p_saf_resp)
})

test_that("the category counts agree with the patients and the scores", {
  trials <- tite_takeda_trials(4, n_trials = 200, seed = 21)

  expect_equal(apply(trials$n_grade, c(1, 2), sum), trials$n_pts,
               ignore_attr = TRUE)
  expect_equal(apply(trials$n_resp_grade, c(1, 2), sum), trials$n_pts,
               ignore_attr = TRUE)
  tox_total <- apply(trials$n_grade, c(1, 2), function(k) {
    sum(k * c(0, 1 / 3, 2 / 3, 1))
  })
  resp_total <- apply(trials$n_resp_grade, c(1, 2), function(k) {
    sum(k * c(0, 0.25 / 3, 1 / 3, 1))
  })
  expect_equal(tox_total, trials$n_tox, ignore_attr = TRUE,
               tolerance = 1e-12)
  expect_equal(resp_total, trials$n_resp, ignore_attr = TRUE,
               tolerance = 1e-12)
  expect_identical(rowSums(trials$n_pts), 3 * trials$cohorts_used)
  expect_true(all(trials$duration > 0))
  # A full trial ends when the longer window has passed after the last
  # arrival.
  full <- trials$cohorts_used == 12L
  expect_true(all(trials$duration[full] >= 60))
  expect_identical(trials$early_stop,
                   trials$stop_reason %in% c("no_admissible_dose",
                                             "no_lower_admissible_dose"))
})

test_that("the simulated trials have the expected averages", {
  # Scenario 4. Expected values as in the first test.
  trials <- tite_takeda_trials(4, n_trials = 1000, seed = 2023)

  expect_equal(unname(colMeans(trials$n_pts)),
               c(4.332, 8.187, 15.702, 5.679, 1.716, 0.339),
               tolerance = 1e-10)
  expect_equal(sum(trials$stop_reason == "no_admissible_dose"), 2)
  expect_equal(sum(trials$stop_reason == "max_cohorts"), 998)
  expect_lt(abs(mean(trials$duration) - 480.258290), 1e-5)

  # Without waiting the cohorts follow each other faster and reach higher
  # doses.
  fast <- tite_takeda_trials(4, n_trials = 1000, seed = 2023,
                             min_completed_share = 0)
  expect_equal(unname(colMeans(fast$n_pts)),
               c(4.326, 7.698, 14.130, 6.078, 2.853, 0.897),
               tolerance = 1e-10)
  expect_lt(abs(mean(fast$duration) - 410.136485), 1e-5)
})

test_that("the trial stops at n_earlystop and when no dose is admissible", {
  # Scenario 7: the lowest dose is the optimal dose. Expected values as in
  # the first test.
  trials <- tite_takeda_trials(7, n_trials = 1000, n_earlystop = 12,
                               seed = 6)

  expect_equal(sum(trials$stop_reason == "no_admissible_dose"), 17)
  expect_equal(sum(trials$stop_reason == "n_earlystop"), 982)
  expect_equal(sum(trials$stop_reason == "max_cohorts"), 1)
  expect_equal(unname(colMeans(trials$n_pts)),
               c(10.719, 3.498, 0.672, 0.048, 0.003, 0),
               tolerance = 1e-10)
  expect_true(all(trials$n_pts <= 12))
  expect_lt(abs(mean(trials$duration) - 179.744192), 1e-5)
})

test_that("two categories with scores 0 and 1 give the TITE-BOIN-ET trials", {
  p <- c(0.05, 0.15, 0.30, 0.45, 0.60)
  q <- c(0.10, 0.30, 0.55, 0.60, 0.62)
  settings <- list(
    list(cs = 3, es = NULL, cor = 0.2, compat = FALSE, start = 1,
         accrual = "uniform", event_time = "weibull", share = 0.5),
    list(cs = 3, es = 9, cor = -0.3, compat = TRUE, start = 2,
         accrual = "exponential", event_time = "uniform", share = 0),
    list(cs = c(1, 2, 3), es = NULL, cor = 0, compat = FALSE, start = 1,
         accrual = "fixed", event_time = "weibull", share = 1)
  )
  for (s in settings) {
    ref <- tite_boin_et_simulate(
      target = 0.30, target_resp = 0.60, p_true = p, p_resp = q,
      n_cohort = 10, cohort_size = s$cs, window = 30, accrual_rate = 0.1,
      resp_window = 45, accrual = s$accrual, event_time = s$event_time,
      resp_cor = s$cor, min_completed_share = s$share, n_trials = 300,
      start_dose = s$start, n_earlystop = s$es,
      boinet_compatible = s$compat, seed = 99
    )
    new <- tite_gboin_et_simulate(
      target = 0.30, target_resp = 0.60, p_grade = rbind(1 - p, p),
      p_resp_grade = rbind(1 - q, q), score = c(0, 1), score_resp = c(0, 1),
      n_cohort = 10, cohort_size = s$cs, window = 30, accrual_rate = 0.1,
      resp_window = 45, accrual = s$accrual, event_time = s$event_time,
      resp_cor = s$cor, min_completed_share = s$share, n_trials = 300,
      start_dose = s$start, n_earlystop = s$es,
      boinet_compatible = s$compat, seed = 99
    )
    expect_identical(new$n_pts, ref$n_pts)
    expect_equal(new$n_tox, ref$n_tox, ignore_attr = TRUE)
    expect_equal(new$n_resp, ref$n_resp, ignore_attr = TRUE)
    expect_identical(new$stop_reason, ref$stop_reason)
    expect_identical(new$cohorts_used, ref$cohorts_used)
    expect_identical(new$early_stop, ref$early_stop)
    expect_equal(new$duration, ref$duration, tolerance = 1e-12)
  }
})

test_that("arrivals longer apart than the windows give the gBOIN-ET trials", {
  # With fixed arrivals 100 days apart every outcome is known at every
  # decision, and no random number is drawn for the arrivals.
  for (s in list(list(k = 1, compat = FALSE, cor = 0.2),
                 list(k = 4, compat = TRUE, cor = 0))) {
    args <- tite_takeda_args(s$k, n_trials = 300, resp_cor = s$cor,
                             boinet_compatible = s$compat, seed = 8)
    new <- do.call(tite_gboin_et_simulate,
                   utils::modifyList(args, list(accrual = "fixed",
                                                accrual_rate = 0.01)))
    drop <- c("window", "resp_window", "accrual_rate")
    ref <- do.call(gboin_et_simulate, args[setdiff(names(args), drop)])
    expect_identical(new$n_pts, ref$n_pts)
    expect_identical(new$n_grade, ref$n_grade)
    expect_identical(new$n_resp_grade, ref$n_resp_grade)
    expect_identical(new$n_tox, ref$n_tox)
    expect_identical(new$n_resp, ref$n_resp)
    expect_identical(new$stop_reason, ref$stop_reason)
  }
})

test_that("fixed arrivals without events give the expected timeline", {
  # Arrivals every 10 days, a 30-day window and no events. Cohort 1 arrives
  # at 0, 10 and 20; the decision waits from day 30 until two of the three
  # patients have finished (day 40) and escalates. Cohort 2 arrives at 40,
  # 50 and 60 and the trial ends 30 days after the last arrival.
  none <- rbind(c(1, 1), c(0, 0))
  args <- list(target = 0.30, target_resp = 0.60, p_grade = none,
               p_resp_grade = none, score = c(0, 1), score_resp = c(0, 1),
               n_cohort = 2, cohort_size = 3, window = 30,
               accrual_rate = 0.1, accrual = "fixed", n_trials = 5, seed = 1)
  tr <- do.call(tite_gboin_et_simulate, args)
  expect_true(all(tr$n_pts[, 1] == 3L & tr$n_pts[, 2] == 3L))
  expect_equal(tr$duration, rep(90, 5))

  # Without waiting cohort 2 arrives at 30, 40 and 50, and the trial ends at
  # 80; waiting for all three patients puts the decision at day 50 and the
  # end at 100.
  tr <- do.call(tite_gboin_et_simulate, c(args, min_completed_share = 0))
  expect_equal(tr$duration, rep(80, 5))
  tr <- do.call(tite_gboin_et_simulate, c(args, min_completed_share = 1))
  expect_equal(tr$duration, rep(100, 5))
})

test_that("only the events seen at a decision are counted", {
  # One dose, cohorts of three, fixed arrivals every 10 days, a 30-day
  # window, every patient with a toxicity event of score 1 at a uniform time
  # within the window, and no efficacy event. The first decision waits until
  # day 40, when the first two patients have finished. If the event of the
  # third patient, who arrived on day 20, has been seen by then, three
  # events in three patients exclude the dose and the trial ends on day 40;
  # otherwise two events in three patients do not, a second cohort is
  # treated, and at least four events in six patients end the trial on day
  # 70. Expected number of trials of each kind from the C++ core run outside
  # of R on R's uniform stream for this seed.
  args <- list(target = 0.30, target_resp = 0.60,
               p_grade = rbind(0, 1), p_resp_grade = rbind(1, 0),
               score = c(0, 1), score_resp = c(0, 1), n_cohort = 3,
               cohort_size = 3, window = 30, accrual_rate = 0.1,
               accrual = "fixed", event_time = "uniform", n_trials = 50,
               seed = 2)
  tr <- do.call(tite_gboin_et_simulate, args)
  short <- tr$n_pts[, 1] == 3L
  expect_equal(sum(short), 31)
  expect_true(all(tr$n_pts[!short, 1] == 6L))
  expect_true(all(tr$n_tox[short, 1] == 3))
  expect_true(all(tr$n_tox[!short, 1] >= 4))
  expect_equal(tr$duration, ifelse(short, 40, 70))
  expect_identical(tr$stop_reason, rep("no_admissible_dose", 50))
  expect_true(all(tr$n_grade[, 1, 1] + tr$n_grade[, 1, 2] == tr$n_pts[, 1]))
})

test_that("the boundaries can be given directly", {
  ref <- boin_et_boundary(0.47 / 1.5, 1.75 / 3)
  default <- tite_takeda_trials(n_trials = 100, seed = 11)
  same <- tite_takeda_trials(n_trials = 100, seed = 11,
                             lambda1 = ref$lambda1, lambda2 = ref$lambda2,
                             eta1 = ref$eta1)
  expect_identical(same$n_grade, default$n_grade)
  expect_identical(same$duration, default$duration)
  expect_true(is.na(same$boundary$prob_incorrect))
  expect_error(tite_takeda_trials(lambda1 = 0.14), "all of")
})

test_that("the trials are reproducible and leave the session seed alone", {
  set.seed(42)
  before <- .Random.seed
  a <- tite_takeda_trials(n_trials = 50, seed = 7)
  expect_identical(.Random.seed, before)
  b <- tite_takeda_trials(n_trials = 50, seed = 7)
  expect_identical(a$n_grade, b$n_grade)
  expect_identical(a$duration, b$duration)
  c2 <- tite_takeda_trials(n_trials = 50, seed = 8)
  expect_false(identical(a$duration, c2$duration))
})

test_that("invalid input is rejected", {
  expect_error(tite_takeda_trials(window = 0), "window")
  expect_error(tite_takeda_trials(accrual_rate = -1), "accrual_rate")
  expect_error(tite_takeda_trials(late_fraction = 1), "late_fraction")
  expect_error(tite_takeda_trials(accrual = "other"))
  expect_error(tite_takeda_trials(min_completed_share = 1.5),
               "min_completed_share")
  expect_error(tite_takeda_trials(score = c(0, 1)), "score")
  expect_error(tite_takeda_trials(start_dose = 7), "start_dose")
  expect_error(tite_takeda_trials(boinet_compatible = NA),
               "boinet_compatible")

  # A first category of probability zero needs uniform event times.
  sure <- rbind(c(0, 0.5), c(1, 0.5))
  args <- list(target = 0.30, target_resp = 0.60, p_grade = sure,
               p_resp_grade = rbind(c(0.5, 0.5), c(0.5, 0.5)),
               score = c(0, 1), score_resp = c(0, 1), n_cohort = 2,
               cohort_size = 3, window = 30, accrual_rate = 0.1,
               n_trials = 5)
  expect_error(do.call(tite_gboin_et_simulate, args), "first rows")
  expect_s3_class(do.call(tite_gboin_et_simulate,
                          c(args, event_time = "uniform")),
                  "boin_et_trials")
})
