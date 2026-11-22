mu_oc <- function(s = 1, ...) {
  sim_gboin(
    target = ets_target(), p_grade = mu_scenario(s), score = ets_score(),
    n_cohort = 10, cohort_size = 3, ...
  )
}

test_that("sim_gboin gives the operating characteristics of the rules", {
  # Expected values from the C++ core run outside of R on R's uniform stream
  # for this seed, and from an implementation of the rules in Python
  oc <- mu_oc(n_trials = 1000, n_earlystop = 100, seed = 2019)

  expect_s3_class(oc, "gboin_oc")
  expect_s3_class(oc, "boin_oc")
  expect_equal(unname(oc$sel_percent), c(0.0, 5.0, 32.2, 53.6, 8.7, 0.4),
               tolerance = 1e-10)
  expect_equal(oc$percent_no_mtd, 0.1, tolerance = 1e-10)
  expect_equal(unname(oc$n_pts_dose),
               c(3.570, 5.751, 9.432, 8.187, 2.697, 0.336),
               tolerance = 1e-10)
  expect_lt(max(abs(unname(oc$n_tox_dose) -
                      c(0.290667, 0.762667, 2.150000, 2.618333, 1.381333,
                        0.230667))), 1e-6)
  expect_equal(oc$total_n_pts, 29.973, tolerance = 1e-10)
  expect_lt(abs(oc$total_n_tox - 7.433667), 1e-6)

  early <- mu_oc(n_trials = 1000, seed = 2019)
  expect_equal(unname(early$sel_percent), c(0.0, 4.6, 33.6, 51.5, 10.2, 0.1),
               tolerance = 1e-10)
  expect_equal(unname(early$n_pts_dose),
               c(3.609, 5.667, 9.027, 8.025, 2.832, 0.399),
               tolerance = 1e-10)
  expect_equal(early$total_n_pts, 29.559, tolerance = 1e-10)
  expect_equal(as.numeric(early$stop_reason_percent["n_earlystop"]), 10.8,
               tolerance = 1e-10)
})

test_that("the two elimination rules give their own results", {
  # Scenario 6 of Mu et al. (2019). Expected values as above.
  posterior <- mu_oc(6, n_trials = 1000, n_earlystop = 100, seed = 6)
  count <- mu_oc(6, n_trials = 1000, n_earlystop = 100, seed = 6,
                 elim_rule = "count")

  expect_equal(unname(posterior$sel_percent),
               c(68.0, 23.4, 1.9, 0.0, 0.0, 0.0), tolerance = 1e-10)
  expect_equal(posterior$percent_no_mtd, 6.7, tolerance = 1e-10)
  expect_equal(unname(posterior$n_pts_dose),
               c(19.740, 7.203, 1.557, 0.171, 0.003, 0.000),
               tolerance = 1e-10)
  expect_equal(unname(count$sel_percent),
               c(70.2, 24.2, 2.3, 0.0, 0.0, 0.0), tolerance = 1e-10)
  expect_equal(count$percent_no_mtd, 3.3, tolerance = 1e-10)
  expect_equal(unname(count$n_pts_dose),
               c(19.917, 7.818, 1.569, 0.186, 0.003, 0.000),
               tolerance = 1e-10)
})

test_that("two categories with scores 0 and 1 give the result of sim_boin", {
  p <- c(0.08, 0.15, 0.31, 0.45, 0.55)
  new <- sim_gboin(0.30, rbind(1 - p, p), c(0, 1), 10, 3, n_trials = 400,
                   seed = 31, keep_trials = TRUE)
  ref <- sim_boin(0.30, p, 10, 3, n_trials = 400, seed = 31,
                  keep_trials = TRUE)

  expect_identical(new$sel_percent, ref$sel_percent)
  expect_identical(new$n_pts_dose, ref$n_pts_dose)
  expect_equal(new$n_tox_dose, ref$n_tox_dose)
  expect_identical(new$percent_no_mtd, ref$percent_no_mtd)
  expect_identical(new$trials$mtd, ref$trials$mtd)
  expect_equal(new$p_true, ref$p_true)
  expect_identical(new$lambda_e, ref$lambda_e)
  expect_identical(new$overdose$doses, ref$overdose$doses)
})

test_that("sim_gboin returns a coherent summary", {
  oc <- mu_oc(n_trials = 300, seed = 4, keep_trials = TRUE)

  expect_equal(sum(oc$sel_percent) + oc$percent_no_mtd, 100,
               tolerance = 1e-10)
  expect_equal(sum(oc$n_pts_dose), oc$total_n_pts, tolerance = 1e-10)
  expect_equal(dim(oc$n_grade_dose), c(4L, 6L))
  expect_equal(unname(colSums(oc$n_grade_dose)), unname(oc$n_pts_dose),
               tolerance = 1e-10)
  expect_equal(unname(colSums(oc$n_grade_dose * c(0, 1 / 3, 2 / 3, 1))),
               unname(oc$n_tox_dose), tolerance = 1e-10)
  expect_identical(rownames(oc$n_grade_dose), paste0("G", 1:4))
  expect_identical(colnames(oc$n_grade_dose), paste0("DL", 1:6))
  expect_identical(oc$design, "gBOIN")
  expect_identical(oc$score, ets_score())
  expect_equal(dim(oc$trials$n_grade), c(300L, 6L, 4L))
  # True scores above the target: doses 4 (0.475 / 1.5 = 0.317) to 6
  expect_identical(oc$overdose$doses, 4:6)
  expect_null(mu_oc(n_trials = 50)$trials)
  reported <- capture_messages(mu_oc(n_trials = 50, verbose = TRUE))
  expect_true(any(grepl("Simulating 50 trials", reported)))
})

test_that("the result prints with the labels of a score", {
  oc <- mu_oc(n_trials = 200, seed = 4)

  expect_output(print(oc), "gBOIN operating characteristics")
  expect_output(print(oc), "target score    : ")
  expect_output(print(oc), "True score \\(%\\)")
  expect_output(print(oc), "Toxicity score")
  expect_output(print(oc), "Doses above a true score of ")
  expect_output(print(oc, percent = TRUE), "Toxicity score \\(%\\)")
  out <- capture.output(print(oc))
  expect_false(any(grepl("DLT", out)))

  trials <- gboin_simulate(ets_target(), mu_scenario(1), ets_score(), 10,
                           3, n_trials = 20)
  expect_output(print(trials), "Simulated gBOIN trials")
  expect_output(print(trials), "Toxicity score")

  res <- oc_results(oc)
  expect_identical(levels(res$design), "gBOIN")
})
