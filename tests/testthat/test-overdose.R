overdose_example <- function(...) {
  sim_boin(
    target = 0.30, p_true = c(0.10, 0.20, 0.30, 0.42, 0.55),
    n_cohort = 20, cohort_size = 3, n_trials = 500, keep_trials = TRUE,
    seed = 31, ...
  )
}

test_that("the overdose summary counts the doses above the cutoff", {
  oc <- overdose_example(overdose_cutoff = 0.33)

  expect_equal(oc$overdose$cutoff, 0.33)
  expect_equal(oc$overdose$doses, c(4L, 5L), ignore_attr = TRUE)
  expect_true(all(c("pct_patients", "pct_patients_by_trial", "avg_n_patients",
                    "pct_trials_any", "pct_trials_over_60",
                    "pct_trials_over_80") %in% names(oc$overdose)))
})

test_that("the percentages agree with the trial level data", {
  oc <- overdose_example(overdose_cutoff = 0.33)

  above <- c(0.10, 0.20, 0.30, 0.42, 0.55) > 0.33
  n_above <- rowSums(oc$trials$n_pts[, above, drop = FALSE])
  n_total <- rowSums(oc$trials$n_pts)

  expect_equal(oc$overdose$pct_patients, sum(n_above) / sum(n_total) * 100)
  expect_equal(oc$overdose$pct_patients_by_trial, mean(n_above / n_total) * 100)
  expect_equal(oc$overdose$avg_n_patients, mean(n_above))
  expect_equal(oc$overdose$pct_trials_any, mean(n_above > 0) * 100)
  expect_equal(oc$overdose$pct_trials_over_60,
               mean(n_above > 0.6 * oc$settings$max_total_pts) * 100)
})

test_that("the cutoff defaults to the target", {
  expect_equal(overdose_example()$overdose$cutoff, 0.30)
  expect_equal(overdose_example()$overdose$doses, c(4L, 5L), ignore_attr = TRUE)
})

test_that("no dose above the cutoff gives zero exposure", {
  oc <- overdose_example(overdose_cutoff = 0.60)

  expect_length(oc$overdose$doses, 0L)
  expect_equal(oc$overdose$pct_patients, 0)
  expect_equal(oc$overdose$pct_patients_by_trial, 0)
  expect_equal(oc$overdose$avg_n_patients, 0)
  expect_equal(oc$overdose$pct_trials_any, 0)
  expect_equal(oc$overdose$pct_trials_over_60, 0)
})

test_that("every dose above the cutoff gives full exposure", {
  oc <- overdose_example(overdose_cutoff = 0.05)

  expect_length(oc$overdose$doses, 5L)
  expect_equal(oc$overdose$pct_patients, 100)
  expect_equal(oc$overdose$pct_patients_by_trial, 100)
  expect_equal(oc$overdose$avg_n_patients, mean(rowSums(oc$trials$n_pts)))
})

test_that("exposure falls as the cutoff rises", {
  exposure <- vapply(
    c(0.15, 0.25, 0.33, 0.45, 0.50),
    function(cut) overdose_example(overdose_cutoff = cut)$overdose$pct_patients,
    numeric(1)
  )
  expect_true(all(diff(exposure) <= 0))
  expect_true(all(exposure >= 0 & exposure <= 100))
})

test_that("a tighter de-escalation boundary lowers exposure above it", {
  # The comparison the design team asks for: keep the default boundary, or
  # tighten it to 0.33, and see how many patients are dosed above 0.33.
  p_tox <- boin_p_tox(target = 0.30, lambda_d = 0.33)
  args <- list(target = 0.30, p_true = c(0.10, 0.20, 0.30, 0.42, 0.55),
               n_cohort = 20, cohort_size = 3, n_trials = 1000,
               overdose_cutoff = 0.33, seed = 32)

  default <- do.call(sim_boin, args)
  tighter <- do.call(sim_boin, c(args, list(p_tox = p_tox)))

  expect_lt(tighter$overdose$pct_patients, default$overdose$pct_patients)
})

test_that("overdose_cutoff is validated", {
  expect_error(overdose_example(overdose_cutoff = 0), "strictly between")
  expect_error(overdose_example(overdose_cutoff = 1), "strictly between")
  expect_error(overdose_example(overdose_cutoff = c(0.2, 0.3)),
               "single finite number")
})
