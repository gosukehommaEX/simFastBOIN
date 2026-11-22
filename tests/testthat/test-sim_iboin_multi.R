iboin_scenarios <- function() {
  list(
    "MTD at dose 3" = c(0.08, 0.15, 0.31, 0.45, 0.55),
    "MTD at dose 4" = c(0.09, 0.12, 0.15, 0.30, 0.45)
  )
}

test_that("sim_iboin_multi runs every scenario with the same design", {
  skeleton <- c(0.10, 0.19, 0.30, 0.42, 0.54)
  res <- sim_iboin_multi(
    target = 0.30, scenarios = iboin_scenarios(), skeleton = skeleton,
    prior_ess = 3, n_cohort = 10, cohort_size = 3, n_trials = 200,
    seed = 17
  )

  expect_s3_class(res, "boin_oc_multi")
  expect_length(res$results, 2L)
  expect_identical(res$scenario_names, c("MTD at dose 3", "MTD at dose 4"))
  expect_identical(res$n_doses, 5L)
  expect_equal(nrow(res$summary_table), 8L)

  for (i in 1:2) {
    single <- sim_iboin(
      target = 0.30, p_true = iboin_scenarios()[[i]], skeleton = skeleton,
      prior_ess = 3, n_cohort = 10, cohort_size = 3, n_trials = 200,
      seed = 17
    )
    expect_identical(res$results[[i]]$sel_percent, single$sel_percent)
    expect_identical(res$results[[i]]$n_pts_dose, single$n_pts_dose)
    expect_s3_class(res$results[[i]], "iboin_oc")
  }
})

test_that("sim_iboin_multi passes the options on", {
  skeleton <- c(0.10, 0.19, 0.30, 0.42, 0.54)
  res <- sim_iboin_multi(
    target = 0.30, scenarios = iboin_scenarios(), skeleton = skeleton,
    prior_ess = 3, n_cohort = 10, cohort_size = 3, n_trials = 100,
    n_earlystop = 100, robust = TRUE, keep_trials = TRUE, seed = 2
  )

  expect_identical(res$results[[1]]$prior_ess, c(3L, 3L, 3L, 0L, 0L))
  expect_false(is.null(res$results[[2]]$trials))
  # n_earlystop = 100 switches the early stop off; a trial can still end
  # before 30 patients when the lowest dose is eliminated
  stops <- res$results[[1]]$trials$stop_reason
  expect_false(any(stops == "n_earlystop"))
  expect_true(all(rowSums(res$results[[1]]$trials$n_pts)[
    stops == "max_sample_size"] == 30))
  reported <- capture_messages(
    sim_iboin_multi(0.30, iboin_scenarios(), skeleton, 3, 10, 3,
                    n_trials = 20, verbose = TRUE)
  )
  expect_true(any(grepl("Scenario 1 of 2", reported)))
  expect_true(any(grepl("Scenario 2 of 2", reported)))
})

test_that("the result prints as iBOIN and enters the comparison tables", {
  res <- sim_iboin_multi(
    target = 0.30, scenarios = iboin_scenarios(),
    skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54), prior_ess = 3,
    n_cohort = 10, cohort_size = 3, n_trials = 100, seed = 2
  )

  expect_output(print(res), "iBOIN operating characteristics across 2")
  expect_output(print(res), "MTD at dose 4")
  expect_identical(levels(oc_results(res)$design), "iBOIN")
  expect_error(
    sim_iboin_multi(0.30, list(c(0.1, 0.2), c(0.1, 0.2, 0.3)),
                    c(0.1, 0.2), 3, 5, 3, n_trials = 10),
    "same number of doses"
  )
})
