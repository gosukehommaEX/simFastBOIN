uboin_scenarios_example <- function() {
  list(
    "Scenario 1" = zhou2019_scenario(1),
    "Scenario 4" = zhou2019_scenario(4),
    list(name = "Scenario 8", p_true = zhou2019_scenario(8)$p_true,
         p_resp = zhou2019_scenario(8)$p_resp)
  )
}

test_that("sim_uboin_multi runs every scenario", {
  oc <- sim_uboin_multi(tox_upper = 0.30, resp_lower = 0.20,
                        scenarios = uboin_scenarios_example(), n_cohort = 18,
                        cohort_size = 3, resp_assoc = 0.2, n_trials = 200,
                        seed = 123)
  expect_s3_class(oc, "uboin_oc_multi")
  expect_length(oc$results, 3L)
  expect_equal(oc$scenario_names, c("Scenario 1", "Scenario 4", "Scenario 8"))
  expect_equal(oc$n_doses, 5L)
  expect_true(all(vapply(oc$results, inherits, logical(1), "uboin_oc")))
})

test_that("each scenario equals a separate call of sim_uboin", {
  oc <- sim_uboin_multi(tox_upper = 0.30, resp_lower = 0.20,
                        scenarios = uboin_scenarios_example(), n_cohort = 18,
                        cohort_size = 3, resp_assoc = 0.2, n_trials = 200,
                        seed = 123)
  sc <- zhou2019_scenario(4)
  single <- sim_uboin(tox_upper = 0.30, resp_lower = 0.20,
                      p_true = sc$p_true, p_resp = sc$p_resp, n_cohort = 18,
                      cohort_size = 3, resp_assoc = 0.2, n_trials = 200,
                      seed = 123)
  expect_equal(oc$results[["Scenario 4"]]$sel_percent, single$sel_percent)
  expect_equal(oc$results[["Scenario 4"]]$n_pts_dose, single$n_pts_dose)
  expect_equal(oc$results[["Scenario 4"]]$percent_no_obd,
               single$percent_no_obd)
})

test_that("the summary table has six rows per scenario", {
  oc <- sim_uboin_multi(tox_upper = 0.30, resp_lower = 0.20,
                        scenarios = uboin_scenarios_example(), n_cohort = 18,
                        cohort_size = 3, n_trials = 100, seed = 1)
  tab <- oc$summary_table
  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 3L * 6L)
  expect_equal(tab$Scenario[c(1, 7, 13)],
               c("Scenario 1", "Scenario 4", "Scenario 8"))
  expect_equal(tab$Item[3], "OBD selected (%)")
  expect_equal(tab[["Total / No OBD"]][3],
               round(oc$results[[1L]]$percent_no_obd, 1))
})

test_that("sim_uboin_multi validates the scenarios", {
  expect_error(
    sim_uboin_multi(tox_upper = 0.30, resp_lower = 0.20,
                    scenarios = list(list(p_true = c(0.1, 0.2))),
                    n_cohort = 18, cohort_size = 3, n_trials = 10),
    "each scenario must be a list with 'p_true' and 'p_resp' elements"
  )
  expect_error(
    sim_uboin_multi(tox_upper = 0.30, resp_lower = 0.20, scenarios = list(),
                    n_cohort = 18, cohort_size = 3, n_trials = 10),
    "'scenarios' must be a non-empty list"
  )
})
