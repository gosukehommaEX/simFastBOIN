waterfall_scenarios_example <- function() {
  list(
    "Scenario 1" = zhang2016_p_true(1),
    "Scenario 4" = zhang2016_p_true(4),
    list(name = "All toxic",
         p_true = matrix(c(0.50, 0.60, 0.70,
                           0.60, 0.70, 0.80), nrow = 2, byrow = TRUE))
  )
}

test_that("sim_waterfall_multi runs every scenario", {
  oc <- sim_waterfall_multi(
    target = 0.30, scenarios = waterfall_scenarios_example(),
    cohort_size = 3, n_trials = 200, seed = 123
  )

  expect_s3_class(oc, "waterfall_oc_multi")
  expect_length(oc$results, 3L)
  expect_equal(oc$scenario_names, c("Scenario 1", "Scenario 4", "All toxic"))
  expect_equal(oc$n_rows, 2L)
  expect_equal(oc$n_cols, 3L)
  expect_true(all(vapply(oc$results, inherits, logical(1), "waterfall_oc")))
})

test_that("each scenario equals a separate call of sim_waterfall", {
  oc <- sim_waterfall_multi(
    target = 0.30, scenarios = waterfall_scenarios_example(),
    cohort_size = 3, n_trials = 200, seed = 123
  )
  single <- sim_waterfall(
    target = 0.30, p_true = zhang2016_p_true(4), cohort_size = 3,
    n_trials = 200, seed = 123
  )

  expect_equal(oc$results[["Scenario 4"]]$sel_percent, single$sel_percent)
  expect_equal(oc$results[["Scenario 4"]]$n_pts_dose, single$n_pts_dose)
  expect_equal(oc$results[["Scenario 4"]]$pct_contour_selected,
               single$pct_contour_selected)
})

test_that("the summary table has one row per scenario", {
  oc <- sim_waterfall_multi(
    target = 0.30, scenarios = waterfall_scenarios_example(),
    cohort_size = 3, n_trials = 200, seed = 123
  )
  tab <- oc$summary_table

  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 3L)
  expect_named(tab, c("Scenario", "MTD contour", "Contour selected (%)",
                      "Patients on contour (%)", "Patients above contour (%)",
                      "No MTD (%)", "Patients", "DLTs"))
  expect_equal(tab[["MTD contour"]], c("A1B3, A2B2", "A1B1", "none"))
  expect_equal(tab[["Contour selected (%)"]][1L],
               round(oc$results[[1L]]$pct_contour_selected, 1))
  expect_true(is.na(tab[["Contour selected (%)"]][3L]))
  expect_equal(tab[["Patients"]],
               round(vapply(oc$results, function(z) z$total_n_pts,
                            numeric(1)), 1), ignore_attr = TRUE)
})

test_that("sim_waterfall_multi validates the scenarios", {
  expect_error(
    sim_waterfall_multi(
      target = 0.30,
      scenarios = list(zhang2016_p_true(1), zhang2016_p_true(12)),
      cohort_size = 3, n_trials = 10
    ),
    "same numbers of levels of both drugs"
  )
  expect_error(
    sim_waterfall_multi(
      target = 0.30,
      scenarios = list(a = zhang2016_p_true(1), a = zhang2016_p_true(4)),
      cohort_size = 3, n_trials = 10
    ),
    "scenario names must be unique"
  )
  expect_error(
    sim_waterfall_multi(target = 0.30, scenarios = list(), cohort_size = 3),
    "'scenarios' must be a non-empty list"
  )
})
