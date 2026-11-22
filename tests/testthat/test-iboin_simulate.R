zhou_design <- function(...) {
  iboin_simulate(
    target = 0.30, p_true = c(0.08, 0.15, 0.31, 0.45, 0.55),
    skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54), prior_ess = 3,
    n_cohort = 10, cohort_size = 3, ...
  )
}

test_that("iboin_simulate returns the trials of the rules", {
  trials <- zhou_design(n_trials = 6, n_earlystop = 100, seed = 7)

  # Expected trials from the C++ core run outside of R on R's uniform stream
  # for this seed, and from an implementation of the rules in Python
  n_pts <- rbind(
    c(3, 21, 6, 0, 0), c(3, 3, 12, 12, 0), c(3, 3, 18, 6, 0),
    c(3, 3, 18, 6, 0), c(3, 3, 12, 6, 6), c(3, 3, 6, 12, 6)
  )
  n_tox <- rbind(
    c(0, 3, 4, 0, 0), c(0, 0, 1, 6, 0), c(0, 0, 4, 4, 0),
    c(1, 0, 3, 4, 0), c(0, 0, 3, 2, 2), c(0, 0, 0, 4, 5)
  )

  expect_s3_class(trials, "iboin_trials")
  expect_s3_class(trials, "boin_trials")
  expect_equal(unname(trials$n_pts), n_pts)
  expect_equal(unname(trials$n_tox), n_tox)
  expect_identical(trials$stop_reason, rep("max_sample_size", 6))
  expect_identical(trials$design, "iBOIN")
  expect_s3_class(trials$boundary, "iboin_boundary")
  expect_identical(colnames(trials$n_pts), paste0("DL", 1:5))
})

test_that("without prior information the trials are those of boin_simulate", {
  settings <- list(
    list(p = c(0.08, 0.15, 0.31, 0.45, 0.55), cs = 3, es = 18,
         rule = "with_stay", xs = FALSE, start = 1),
    list(p = c(0.45, 0.55, 0.65, 0.75, 0.85), cs = 3, es = 12,
         rule = "simple", xs = TRUE, start = 2),
    list(p = c(0.05, 0.08, 0.10, 0.14, 0.30), cs = c(1, 2, 3), es = 9,
         rule = "with_stay", xs = FALSE, start = 1)
  )
  for (s in settings) {
    ref <- boin_simulate(
      target = 0.30, p_true = s$p, n_cohort = 10, cohort_size = s$cs,
      n_trials = 300, start_dose = s$start, n_earlystop = s$es,
      extrasafe = s$xs, n_earlystop_rule = s$rule, seed = 99
    )
    new <- iboin_simulate(
      target = 0.30, p_true = s$p, skeleton = c(0.1, 0.2, 0.3, 0.4, 0.5),
      prior_ess = 0, n_cohort = 10, cohort_size = s$cs, n_trials = 300,
      start_dose = s$start, n_earlystop = s$es, extrasafe = s$xs,
      n_earlystop_rule = s$rule, seed = 99
    )
    expect_identical(new$n_pts, ref$n_pts)
    expect_identical(new$n_tox, ref$n_tox)
    expect_identical(new$eliminated, ref$eliminated)
    expect_identical(new$stop_reason, ref$stop_reason)
    expect_identical(new$cohorts_used, ref$cohorts_used)
  }
})

test_that("the prior information changes the trials", {
  informed <- zhou_design(n_trials = 300, seed = 5)
  plain <- boin_simulate(
    target = 0.30, p_true = c(0.08, 0.15, 0.31, 0.45, 0.55),
    n_cohort = 10, cohort_size = 3, n_trials = 300, seed = 5
  )
  expect_false(identical(informed$n_pts, plain$n_pts))
})

test_that("iboin_simulate respects the design limits", {
  trials <- zhou_design(n_trials = 300, n_earlystop = 9, seed = 3)

  expect_true(all(rowSums(trials$n_pts) <= 30))
  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(trials$stop_reason %in%
                    c("lowest_dose_eliminated", "lowest_dose_too_toxic",
                      "n_earlystop", "max_sample_size", "max_cohorts")))
  stopped <- trials$stop_reason == "n_earlystop"
  expect_true(any(stopped))
  expect_true(all(apply(trials$n_pts[stopped, , drop = FALSE], 1, max) >= 9))
  expect_identical(trials$settings$prior_ess, rep(3L, 5))
})

test_that("iboin_simulate is reproducible and leaves the session seed alone", {
  set.seed(42)
  before <- runif(1)
  set.seed(42)
  a <- zhou_design(n_trials = 50, seed = 8)
  after <- runif(1)
  b <- zhou_design(n_trials = 50, seed = 8)

  expect_identical(a$n_pts, b$n_pts)
  expect_identical(before, after)
})

test_that("iboin_simulate rejects inconsistent arguments", {
  expect_error(
    iboin_simulate(0.30, c(0.1, 0.3, 0.5), c(0.1, 0.3), 3, 5, 3),
    "skeleton"
  )
  expect_error(zhou_design(start_dose = 6), "start_dose")
  expect_error(zhou_design(n_earlystop_rule = "other"), "should be one of")
  expect_warning(zhou_design(n_trials = 10, n_earlystop = 6), "n_earlystop")
})
