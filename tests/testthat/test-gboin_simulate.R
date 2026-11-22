mu_trials <- function(s = 1, ...) {
  gboin_simulate(
    target = ets_target(), p_grade = mu_scenario(s), score = ets_score(),
    n_cohort = 10, cohort_size = 3, ...
  )
}

test_that("gboin_simulate returns the trials of the rules", {
  trials <- mu_trials(n_trials = 4, n_earlystop = 100, seed = 3)

  # Expected trials from the C++ core run outside of R on R's uniform stream
  # for this seed, and from an implementation of the rules in Python
  n_pts <- rbind(c(3, 3, 9, 9, 6, 0), c(3, 3, 18, 6, 0, 0),
                 c(3, 3, 3, 18, 3, 0), c(3, 9, 18, 0, 0, 0))
  ets_total <- rbind(c(0.5, 0.0, 2.0, 4.5, 4.5, 0.0),
                     c(0.0, 0.0, 5.0, 6.0, 0.0, 0.0),
                     c(0.5, 0.5, 0.5, 4.5, 3.5, 0.0),
                     c(1.0, 2.0, 7.5, 0.0, 0.0, 0.0))
  grades_trial_1 <- rbind(c(2, 1, 0, 0), c(3, 0, 0, 0), c(6, 2, 1, 0),
                          c(6, 0, 0, 3), c(3, 0, 0, 3), c(0, 0, 0, 0))

  expect_s3_class(trials, "gboin_trials")
  expect_s3_class(trials, "boin_trials")
  expect_equal(unname(trials$n_pts), n_pts)
  expect_equal(unname(trials$n_tox) * 1.5, ets_total, tolerance = 1e-12)
  expect_equal(unname(trials$n_grade[1, , ]), grades_trial_1)
  expect_identical(trials$stop_reason, rep("max_sample_size", 4))
  expect_identical(trials$design, "gBOIN")
  expect_identical(dimnames(trials$n_grade)[[2]], paste0("DL", 1:6))
  expect_identical(dimnames(trials$n_grade)[[3]], paste0("G", 1:4))
})

test_that("the settings hold the true scores and the BOIN boundaries", {
  trials <- mu_trials(n_trials = 10)

  # Mean normalized scores of the doses, by hand: sum of p times score / 1.5
  expect_equal(trials$settings$p_true,
               c(0.115, 0.19, 0.335, 0.475, 0.755, 1.05) / 1.5,
               tolerance = 1e-12)
  ref <- boin_lambda(ets_target())
  expect_identical(trials$lambda_e, ref$lambda_e)
  expect_identical(trials$lambda_d, ref$lambda_d)
  expect_equal(trials$lambda_e, 0.2470995617, tolerance = 1e-9)
  expect_equal(trials$lambda_d, 0.3745942485, tolerance = 1e-9)
  expect_equal(trials$settings$score_norm, c(0, 1 / 3, 2 / 3, 1))
})

test_that("the grade counts agree with the patients and the scores", {
  trials <- mu_trials(n_trials = 200, seed = 21)

  expect_equal(apply(trials$n_grade, c(1, 2), sum), trials$n_pts,
               ignore_attr = TRUE)
  score_total <- apply(trials$n_grade, c(1, 2), function(k) {
    sum(k * c(0, 1 / 3, 2 / 3, 1))
  })
  expect_equal(score_total, trials$n_tox, ignore_attr = TRUE,
               tolerance = 1e-12)
  expect_true(all(trials$n_tox <= trials$n_pts + 1e-12))
  expect_true(all(rowSums(trials$n_pts) <= 30))
})

test_that("scores that differ by a factor give the same trials", {
  # The totals of the scores 0, 0.5, 1, 1.5 and of the scores 0, 1, 2, 3 are
  # exact, and each total is divided once by the largest score.
  half <- mu_trials(n_trials = 200, seed = 21)
  whole <- gboin_simulate(
    target = ets_target(), p_grade = mu_scenario(1), score = c(0, 1, 2, 3),
    n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 21
  )
  expect_identical(whole$n_pts, half$n_pts)
  expect_identical(whole$n_grade, half$n_grade)
  expect_identical(whole$n_tox, half$n_tox)
  expect_identical(whole$stop_reason, half$stop_reason)
})

test_that("two categories with scores 0 and 1 give the BOIN trials", {
  settings <- list(
    list(p = c(0.08, 0.15, 0.31, 0.45, 0.55), target = 0.30, cs = 3,
         es = 18, rule = "with_stay", xs = FALSE, start = 1),
    list(p = c(0.45, 0.55, 0.65, 0.75, 0.85), target = 0.30, cs = 3,
         es = 12, rule = "simple", xs = TRUE, start = 2),
    list(p = c(0.02, 0.05, 0.10, 0.18, 0.30), target = 0.10,
         cs = c(1, 2, 3), es = 9, rule = "with_stay", xs = FALSE,
         start = 1)
  )
  for (s in settings) {
    ref <- boin_simulate(
      target = s$target, p_true = s$p, n_cohort = 10, cohort_size = s$cs,
      n_trials = 300, start_dose = s$start, n_earlystop = s$es,
      extrasafe = s$xs, n_earlystop_rule = s$rule, seed = 99
    )
    for (rule in c("posterior", "count")) {
      new <- gboin_simulate(
        target = s$target, p_grade = rbind(1 - s$p, s$p), score = c(0, 1),
        n_cohort = 10, cohort_size = s$cs, n_trials = 300,
        start_dose = s$start, n_earlystop = s$es, extrasafe = s$xs,
        elim_rule = rule, n_earlystop_rule = s$rule, seed = 99
      )
      expect_identical(new$n_pts, ref$n_pts)
      expect_equal(new$n_tox, ref$n_tox, ignore_attr = TRUE)
      expect_identical(new$eliminated, ref$eliminated)
      expect_identical(new$stop_reason, ref$stop_reason)
      expect_identical(new$cohorts_used, ref$cohorts_used)
    }
  }
})

test_that("elimination by count stops fewer trials than by posterior", {
  # Scenario 6: the lowest dose is close to the target. Expected numbers of
  # stopped trials as in the first test.
  posterior <- mu_trials(6, n_trials = 1000, n_earlystop = 100, seed = 6)
  count <- mu_trials(6, n_trials = 1000, n_earlystop = 100, seed = 6,
                     elim_rule = "count")

  expect_equal(sum(posterior$stop_reason == "lowest_dose_eliminated"), 63)
  expect_equal(sum(count$stop_reason == "lowest_dose_eliminated"), 29)
  expect_identical(posterior$settings$elim_rule, "posterior")
  expect_identical(count$settings$elim_rule, "count")
})

test_that("gboin_simulate is reproducible and leaves the session seed alone", {
  set.seed(42)
  before <- runif(1)
  set.seed(42)
  a <- mu_trials(n_trials = 50, seed = 8)
  after <- runif(1)
  b <- mu_trials(n_trials = 50, seed = 8)

  expect_identical(a$n_grade, b$n_grade)
  expect_identical(before, after)
})

test_that("row names of p_grade name the toxicity categories", {
  p <- mu_scenario(1)
  rownames(p) <- c("Grade 0-1", "Grade 2", "Grade 3", "Grade 4")
  trials <- gboin_simulate(ets_target(), p, ets_score(), 5, 3, n_trials = 5)
  expect_identical(dimnames(trials$n_grade)[[3]], rownames(p))

  # Row names that are only partly given are not used
  rownames(p) <- c("", "Grade 2", "Grade 3", "Grade 4")
  trials <- gboin_simulate(ets_target(), p, ets_score(), 5, 3, n_trials = 5)
  expect_identical(dimnames(trials$n_grade)[[3]], paste0("G", 1:4))
})

test_that("gboin_simulate rejects invalid grades and scores", {
  p <- mu_scenario(1)
  bad <- p
  bad[1, 1] <- 0.5
  expect_error(gboin_simulate(ets_target(), bad, ets_score(), 5, 3),
               "sum to 1")
  expect_error(gboin_simulate(ets_target(), p[1, ], ets_score(), 5, 3),
               "p_grade")
  expect_error(gboin_simulate(ets_target(), p, c(0, 1), 5, 3), "score")
  expect_error(gboin_simulate(ets_target(), p, c(0, -1, 1, 2), 5, 3),
               "score")
  expect_error(gboin_simulate(ets_target(), p, c(0, 0, 0, 0), 5, 3),
               "score")
  expect_error(mu_trials(start_dose = 7), "start_dose")
  expect_error(mu_trials(elim_rule = "other"), "should be one of")
  expect_error(gboin_simulate(0.47, p, ets_score(), 5, 3, p_saf = 0.46),
               "p_saf")
  expect_warning(mu_trials(n_trials = 10, n_earlystop = 6), "n_earlystop")
})
