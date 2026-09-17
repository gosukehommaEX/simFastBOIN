scenarios_example <- function() {
  list(
    list(name = "MTD at dose 3", p_true = c(0.05, 0.15, 0.30, 0.45, 0.60)),
    list(name = "MTD at dose 1", p_true = c(0.30, 0.45, 0.55, 0.65, 0.75)),
    list(name = "All safe",      p_true = c(0.02, 0.04, 0.06, 0.08, 0.10))
  )
}

test_that("sim_boin_multi runs every scenario", {
  oc <- sim_boin_multi(
    target = 0.30, scenarios = scenarios_example(),
    n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 123
  )

  expect_s3_class(oc, "boin_oc_multi")
  expect_length(oc$results, 3L)
  expect_equal(oc$scenario_names,
               c("MTD at dose 3", "MTD at dose 1", "All safe"))
  expect_equal(oc$n_doses, 5L)
  expect_true(all(vapply(oc$results, inherits, logical(1), "boin_oc")))
})

test_that("the summary table has one block of four rows per scenario", {
  oc <- sim_boin_multi(
    target = 0.30, scenarios = scenarios_example(),
    n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 123
  )

  expect_s3_class(oc$summary_table, "data.frame")
  expect_equal(nrow(oc$summary_table), 12L)
  expect_equal(ncol(oc$summary_table), 8L)
  expect_equal(oc$summary_table$Scenario[c(1, 5, 9)], oc$scenario_names)
  expect_equal(unique(oc$summary_table$Item),
               c("True DLT rate (%)", "MTD selected (%)",
                 "Patients treated", "Patients with DLT"))
})

test_that("a scenario matches the same scenario run on its own", {
  scenarios <- scenarios_example()
  multi <- sim_boin_multi(
    target = 0.30, scenarios = scenarios,
    n_cohort = 10, cohort_size = 3, n_trials = 300, seed = 42
  )

  for (i in seq_along(scenarios)) {
    single <- sim_boin(
      target = 0.30, p_true = scenarios[[i]]$p_true,
      n_cohort = 10, cohort_size = 3, n_trials = 300, seed = 42
    )
    expect_equal(multi$results[[i]]$sel_percent, single$sel_percent)
    expect_equal(multi$results[[i]]$n_pts_dose, single$n_pts_dose)
  }
})

test_that("a named list of probability vectors is accepted", {
  oc <- sim_boin_multi(
    target = 0.30,
    scenarios = list(low = c(0.05, 0.15, 0.30), high = c(0.30, 0.45, 0.60)),
    n_cohort = 8, cohort_size = 3, n_trials = 100, seed = 1
  )
  expect_equal(oc$scenario_names, c("low", "high"))
})

test_that("design arguments reach every scenario", {
  oc <- sim_boin_multi(
    target = 0.30, scenarios = scenarios_example(),
    n_cohort = 12, cohort_size = 2, n_trials = 100,
    extrasafe = TRUE, titration = TRUE, bound_mtd = TRUE, seed = 1
  )

  for (result in oc$results) {
    expect_equal(result$settings$max_total_pts, 24)
    expect_true(result$settings$extrasafe)
    expect_true(result$settings$titration)
  }
})

test_that("sim_boin_multi validates the scenario list", {
  expect_error(sim_boin_multi(0.30, list(), 10, 3), "non-empty list")
  expect_error(
    sim_boin_multi(0.30, list(list(name = "a", p_true = c(0.1, 0.2)),
                              list(name = "b", p_true = c(0.1, 0.2, 0.3))),
                   10, 3),
    "same number of doses"
  )
  expect_error(
    sim_boin_multi(0.30, list(list(name = "a", p_true = c(0.1, 0.2)),
                              list(name = "a", p_true = c(0.1, 0.2))),
                   10, 3),
    "unique"
  )
  expect_error(sim_boin_multi(0.30, list(list(name = "a")), 10, 3), "p_true")
})
