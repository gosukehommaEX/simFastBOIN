comb_scenarios_example <- function() {
  list(
    list(name = "Two MTD combinations",
         p_true = matrix(c(0.05, 0.10, 0.20,
                           0.10, 0.20, 0.30,
                           0.20, 0.30, 0.45), nrow = 3, byrow = TRUE)),
    list(name = "MTD at the top",
         p_true = matrix(c(0.02, 0.04, 0.06,
                           0.04, 0.06, 0.10,
                           0.06, 0.10, 0.30), nrow = 3, byrow = TRUE)),
    list(name = "All toxic",
         p_true = matrix(c(0.40, 0.50, 0.60,
                           0.50, 0.60, 0.70,
                           0.60, 0.70, 0.80), nrow = 3, byrow = TRUE))
  )
}

test_that("sim_comb_boin_multi runs every scenario", {
  oc <- sim_comb_boin_multi(
    target = 0.30, scenarios = comb_scenarios_example(),
    n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 123
  )

  expect_s3_class(oc, "comb_boin_oc_multi")
  expect_length(oc$results, 3L)
  expect_equal(oc$scenario_names,
               c("Two MTD combinations", "MTD at the top", "All toxic"))
  expect_equal(oc$n_rows, 3L)
  expect_equal(oc$n_cols, 3L)
  expect_true(all(vapply(oc$results, inherits, logical(1), "comb_boin_oc")))
})

test_that("the summary table has one row per scenario", {
  oc <- sim_comb_boin_multi(
    target = 0.30, scenarios = comb_scenarios_example(),
    n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 123
  )
  tab <- oc$summary_table

  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 3L)
  expect_equal(tab$Scenario, oc$scenario_names)
  expect_equal(tab[["True MTD"]], c("A3B2, A2B3", "A3B3", "A1B1"))
  expect_equal(tab[["MTD selected (%)"]],
               vapply(oc$results, function(r) round(r$pct_true_mtd_selected, 1),
                      numeric(1), USE.NAMES = FALSE))
  expect_equal(tab[["Patients"]],
               vapply(oc$results, function(r) round(r$total_n_pts, 1),
                      numeric(1), USE.NAMES = FALSE))
})

test_that("a scenario matches the same scenario run on its own", {
  scenarios <- comb_scenarios_example()
  multi <- sim_comb_boin_multi(
    target = 0.30, scenarios = scenarios,
    n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 42
  )

  for (i in seq_along(scenarios)) {
    single <- sim_comb_boin(
      target = 0.30, p_true = scenarios[[i]]$p_true,
      n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 42
    )
    expect_equal(multi$results[[i]]$sel_percent, single$sel_percent)
    expect_equal(multi$results[[i]]$n_pts_dose, single$n_pts_dose)
  }
})

test_that("a named list of matrices is accepted", {
  oc <- sim_comb_boin_multi(
    target = 0.30,
    scenarios = list(low = matrix(c(0.05, 0.10, 0.10, 0.30), 2),
                     high = matrix(c(0.30, 0.45, 0.45, 0.60), 2)),
    n_cohort = 6, cohort_size = 3, n_trials = 50, seed = 1
  )

  expect_equal(oc$scenario_names, c("low", "high"))
})

test_that("scenarios of different sizes or with repeated names are rejected", {
  expect_error(
    sim_comb_boin_multi(
      target = 0.30,
      scenarios = list(a = matrix(0.1, 2, 2), b = matrix(0.1, 2, 3)),
      n_cohort = 6, cohort_size = 3, n_trials = 10
    ),
    "same numbers of levels"
  )
  expect_error(
    sim_comb_boin_multi(
      target = 0.30,
      scenarios = list(list(name = "x", p_true = matrix(0.1, 2, 2)),
                       list(name = "x", p_true = matrix(0.2, 2, 2))),
      n_cohort = 6, cohort_size = 3, n_trials = 10
    ),
    "unique"
  )
  expect_error(
    sim_comb_boin_multi(target = 0.30, scenarios = list(c(0.1, 0.2)),
                        n_cohort = 6, cohort_size = 3, n_trials = 10),
    "matrix"
  )
})
