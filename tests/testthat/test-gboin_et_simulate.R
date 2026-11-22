test_that("gboin_et_simulate returns the trials of the rules", {
  trials <- takeda_trials(n_trials = 4, seed = 3)

  # Expected trials from the C++ core run outside of R on R's uniform stream
  # for this seed, and from an implementation of the rules in Python
  n_pts <- rbind(c(3, 3, 9, 15, 0, 0), c(3, 3, 12, 12, 0, 0),
                 c(3, 3, 21, 3, 0, 0), c(3, 9, 12, 6, 0, 0))
  ets_total <- rbind(c(0.5, 1.5, 3.5, 7.5, 0, 0), c(0.0, 0.5, 1.5, 6.5, 0, 0),
                     c(0.5, 0.5, 4.5, 4.0, 0, 0), c(0.0, 1.5, 6.0, 3.5, 0, 0))
  ees_total <- rbind(c(0, 1.5, 3.5, 13.0, 0, 0), c(0, 2.0, 5.5, 14.5, 0, 0),
                     c(0, 0.0, 8.5, 4.0, 0, 0), c(0, 3.0, 4.0, 5.5, 0, 0))
  grades_trial_1 <- rbind(c(2, 1, 0, 0), c(1, 1, 1, 0), c(5, 2, 1, 1),
                          c(8, 3, 0, 4), c(0, 0, 0, 0), c(0, 0, 0, 0))
  responses_trial_1 <- rbind(c(3, 0, 0), c(2, 0, 1), c(6, 2, 1), c(6, 1, 8),
                             c(0, 0, 0), c(0, 0, 0))

  expect_s3_class(trials, "boin_et_trials")
  expect_identical(trials$design, "gBOIN-ET")
  expect_equal(unname(trials$n_pts), n_pts)
  expect_equal(unname(trials$n_tox) * 1.5, ets_total, tolerance = 1e-12)
  expect_equal(unname(trials$n_resp) * 1.5, ees_total, tolerance = 1e-12)
  expect_equal(unname(trials$n_grade[1, , ]), grades_trial_1)
  expect_equal(unname(trials$n_resp_grade[1, , ]), responses_trial_1)
  expect_identical(trials$stop_reason, rep("max_cohorts", 4))
  expect_identical(trials$early_stop, rep(FALSE, 4))
  expect_identical(trials$cohorts_used, rep(10L, 4))
  expect_identical(dimnames(trials$n_grade)[[2]], paste0("DL", 1:6))
  expect_identical(dimnames(trials$n_grade)[[3]], paste0("G", 1:4))
  expect_identical(dimnames(trials$n_resp_grade)[[3]], paste0("R", 1:3))
})

test_that("the settings hold the true scores and the BOIN-ET boundaries", {
  trials <- takeda_trials(n_trials = 10)

  # Mean normalized scores of the doses, by hand: sum of p times score / 1.5
  expect_equal(trials$settings$p_true,
               c(0.115, 0.19, 0.335, 0.475, 0.755, 1.05) / 1.5,
               tolerance = 1e-12)
  expect_equal(trials$settings$p_resp,
               c(0.125, 0.35, 0.575, 0.95, 1.025, 1.10) / 1.5,
               tolerance = 1e-12)
  ref <- boin_et_boundary(0.47 / 1.5, 0.95 / 1.5)
  expect_identical(trials$boundary$lambda1, ref$lambda1)
  expect_identical(trials$boundary$lambda2, ref$lambda2)
  expect_identical(trials$boundary$eta1, ref$eta1)
  # The values of the program of Takeda et al. (2022): 0.14133, 0.37333, 0.50
  expect_equal(c(ref$lambda1, ref$lambda2, ref$eta1),
               c(0.47 / 1.5 * 0.1 + 0.11, 0.47 / 1.5 + 0.06, 0.50),
               tolerance = 1e-12)
  expect_equal(trials$settings$score_norm, c(0, 1 / 3, 2 / 3, 1))
  expect_equal(trials$settings$score_resp_norm, c(0, 2 / 3, 1))
  expect_identical(trials$settings$p_eli, 0.47 / 1.5)
  expect_identical(trials$settings$p_eli_resp, ref$p_saf_resp)
})

test_that("the category counts agree with the patients and the scores", {
  trials <- takeda_trials(n_trials = 200, seed = 21)

  expect_equal(apply(trials$n_grade, c(1, 2), sum), trials$n_pts,
               ignore_attr = TRUE)
  expect_equal(apply(trials$n_resp_grade, c(1, 2), sum), trials$n_pts,
               ignore_attr = TRUE)
  tox_total <- apply(trials$n_grade, c(1, 2), function(k) {
    sum(k * c(0, 1 / 3, 2 / 3, 1))
  })
  resp_total <- apply(trials$n_resp_grade, c(1, 2), function(k) {
    sum(k * c(0, 2 / 3, 1))
  })
  expect_equal(tox_total, trials$n_tox, ignore_attr = TRUE,
               tolerance = 1e-12)
  expect_equal(resp_total, trials$n_resp, ignore_attr = TRUE,
               tolerance = 1e-12)
  expect_true(all(trials$n_tox <= trials$n_pts + 1e-12))
  expect_true(all(trials$n_resp <= trials$n_pts + 1e-12))
  expect_true(all(rowSums(trials$n_pts) <= 30))
  expect_identical(trials$early_stop,
                   trials$stop_reason %in% c("no_admissible_dose",
                                             "no_lower_admissible_dose"))
})

test_that("two categories with scores 0 and 1 give the BOIN-ET trials", {
  p <- c(0.05, 0.15, 0.30, 0.45, 0.60)
  q <- c(0.10, 0.30, 0.55, 0.60, 0.62)
  settings <- list(
    list(cs = 3, es = NULL, cor = 0.2, compat = FALSE, start = 1),
    list(cs = 3, es = 9, cor = -0.3, compat = TRUE, start = 2),
    list(cs = c(1, 2, 3), es = NULL, cor = 0, compat = FALSE, start = 1)
  )
  for (s in settings) {
    ref <- boin_et_simulate(
      target = 0.30, target_resp = 0.60, p_true = p, p_resp = q,
      n_cohort = 10, cohort_size = s$cs, n_trials = 300,
      start_dose = s$start, n_earlystop = s$es, resp_cor = s$cor,
      boinet_compatible = s$compat, seed = 99
    )
    new <- gboin_et_simulate(
      target = 0.30, target_resp = 0.60, p_grade = rbind(1 - p, p),
      p_resp_grade = rbind(1 - q, q), score = c(0, 1), score_resp = c(0, 1),
      n_cohort = 10, cohort_size = s$cs, n_trials = 300,
      start_dose = s$start, n_earlystop = s$es, resp_cor = s$cor,
      boinet_compatible = s$compat, seed = 99
    )
    expect_identical(new$n_pts, ref$n_pts)
    expect_equal(new$n_tox, ref$n_tox, ignore_attr = TRUE)
    expect_equal(new$n_resp, ref$n_resp, ignore_attr = TRUE)
    expect_identical(new$stop_reason, ref$stop_reason)
    expect_identical(new$cohorts_used, ref$cohorts_used)
    expect_identical(new$early_stop, ref$early_stop)
  }
})

test_that("the boundaries can be given directly", {
  ref <- boin_et_boundary(0.47 / 1.5, 0.95 / 1.5)
  default <- takeda_trials(n_trials = 100, seed = 11)
  same <- takeda_trials(n_trials = 100, seed = 11, lambda1 = ref$lambda1,
                        lambda2 = ref$lambda2, eta1 = ref$eta1)
  expect_identical(same$n_grade, default$n_grade)
  expect_identical(same$n_resp_grade, default$n_resp_grade)
  expect_true(is.na(same$boundary$prob_incorrect))

  # Expected numbers of patients as in the first test: the rounded
  # boundaries of the article give other trials
  rounded <- takeda_trials(n_trials = 1000, seed = 2022, lambda1 = 0.14,
                           lambda2 = 0.37, eta1 = 0.50)
  expect_identical(rounded$boundary$eta1, 0.50)
  expect_equal(unname(colMeans(rounded$n_pts)),
               c(3.396, 5.511, 8.940, 10.491, 1.530, 0.114),
               tolerance = 1e-10)
  expect_equal(sum(rounded$early_stop), 3)

  expect_error(takeda_trials(lambda1 = 0.14), "all of")
  expect_error(takeda_trials(lambda1 = 0.4, lambda2 = 0.3, eta1 = 0.5),
               "lambda1")
})

test_that("the trial stops at n_earlystop and when no dose is admissible", {
  # Scenario 6: the lowest dose is close to the target toxicity score.
  # Expected numbers of trials as in the first test.
  trials <- takeda_trials(6, n_trials = 1000, n_earlystop = 12, seed = 6)

  expect_equal(sum(trials$stop_reason == "no_admissible_dose"), 51)
  expect_equal(sum(trials$stop_reason == "n_earlystop"), 948)
  expect_equal(sum(trials$stop_reason == "max_cohorts"), 1)
  expect_equal(unname(colMeans(trials$n_pts)),
               c(10.965, 2.205, 0.213, 0.003, 0, 0), tolerance = 1e-10)
  expect_true(all(trials$n_pts <= 12))
})

test_that("the trials are reproducible and leave the session seed alone", {
  set.seed(42)
  before <- runif(1)
  set.seed(42)
  a <- takeda_trials(n_trials = 50, seed = 8)
  after <- runif(1)
  b <- takeda_trials(n_trials = 50, seed = 8)

  expect_identical(a$n_grade, b$n_grade)
  expect_identical(a$n_resp_grade, b$n_resp_grade)
  expect_identical(before, after)
})

test_that("row names name the categories", {
  sc <- takeda_scenario(1)
  rownames(sc$p_grade) <- c("Grade 0-1", "Grade 2", "Grade 3", "Grade 4")
  rownames(sc$p_resp_grade) <- c("SD or PD", "PR", "CR")
  trials <- gboin_et_simulate(0.47 / 1.5, 0.95 / 1.5, sc$p_grade,
                              sc$p_resp_grade, c(0, 0.5, 1, 1.5),
                              c(0, 1, 1.5), 5, 3, n_trials = 5)
  expect_identical(dimnames(trials$n_grade)[[3]], rownames(sc$p_grade))
  expect_identical(dimnames(trials$n_resp_grade)[[3]],
                   rownames(sc$p_resp_grade))

  # Row names that are only partly given are not used
  sc <- takeda_scenario(1)
  expect_identical(rownames(sc$p_resp_grade), c("", "pr", "cr"))
  trials <- gboin_et_simulate(0.47 / 1.5, 0.95 / 1.5, sc$p_grade,
                              sc$p_resp_grade, c(0, 0.5, 1, 1.5),
                              c(0, 1, 1.5), 5, 3, n_trials = 5)
  expect_identical(dimnames(trials$n_resp_grade)[[3]], paste0("R", 1:3))
})

test_that("gboin_et_simulate rejects invalid categories and scores", {
  sc <- takeda_scenario(1)
  run <- function(p_grade = sc$p_grade, p_resp_grade = sc$p_resp_grade,
                  score = c(0, 0.5, 1, 1.5), score_resp = c(0, 1, 1.5),
                  ...) {
    gboin_et_simulate(0.47 / 1.5, 0.95 / 1.5, p_grade, p_resp_grade, score,
                      score_resp, 5, 3, n_trials = 5, ...)
  }
  bad <- sc$p_resp_grade
  bad[1, 1] <- 0.5
  expect_error(run(p_resp_grade = bad), "sum to 1")
  expect_error(run(p_resp_grade = sc$p_resp_grade[, 1:5]), "one column")
  expect_error(run(p_resp_grade = sc$p_resp_grade[1, ]), "p_resp_grade")
  expect_error(run(score_resp = c(0, 1)), "score_resp")
  expect_error(run(score_resp = c(0, 0, 0)), "score_resp")
  expect_error(run(score = c(0, 1)), "score")
  expect_error(run(p_grade = sc$p_grade[1, ]), "p_grade")
  expect_error(run(start_dose = 7), "start_dose")
  expect_error(run(resp_cor = 1.5), "resp_cor")
  expect_error(run(boinet_compatible = NA), "boinet_compatible")
})
