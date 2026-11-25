uboin_oc_args <- function(scenario, ...) {
  sc <- zhou2019_scenario(scenario)
  args <- list(tox_upper = 0.30, resp_lower = 0.20, p_true = sc$p_true,
               p_resp = sc$p_resp, n_cohort = 18, cohort_size = 3,
               resp_assoc = 0.2)
  given <- list(...)
  args[names(given)] <- given
  args
}

test_that("sim_uboin summarizes the trials and their OBDs", {
  oc <- do.call(sim_uboin, uboin_oc_args(4, n_trials = 400,
                                         keep_trials = TRUE, seed = 1))
  expect_s3_class(oc, "uboin_oc")
  expect_equal(oc$design, "U-BOIN")
  expect_lt(abs(sum(oc$sel_percent) + oc$percent_no_obd - 100), 1e-10)
  expect_lt(abs(sum(oc$stop_reason_percent) - 100), 1e-10)
  expect_equal(oc$total_n_pts, mean(rowSums(oc$trials$n_pts)))
  expect_equal(oc$total_n_resp, mean(rowSums(oc$trials$n_resp)))
  expect_equal(oc$percent_stage2, 100 * mean(oc$trials$stage2))

  trials <- do.call(uboin_simulate, uboin_oc_args(4, n_trials = 400,
                                                  seed = 1))
  sel <- uboin_select_obd(trials$n_outcome, tox_upper = 0.30,
                          resp_lower = 0.20, eliminated = trials$eliminated)
  expect_identical(oc$trials$n_outcome, trials$n_outcome)
  expect_identical(oc$obd, sel$obd)
  expect_true(all(is.na(oc$obd[trials$stop_reason %in%
                                 c("lowest_dose_eliminated",
                                   "no_admissible_dose")])))

  plain <- do.call(sim_uboin, uboin_oc_args(4, n_trials = 400, seed = 1))
  expect_null(plain$trials)
  expect_equal(plain$sel_percent, oc$sel_percent)
})

test_that("the true utilities are those of Table 3 of the article", {
  for (scenario in c(1, 4, 8)) {
    oc <- do.call(sim_uboin, uboin_oc_args(scenario, n_trials = 10,
                                           seed = 1))
    expect_equal(unname(round(oc$true_utility)),
                 zhou2019_table3(scenario)$utility,
                 info = paste("Scenario", scenario))
  }
})

test_that("sim_uboin agrees with Table 3 of the article", {
  # Four standard errors of the difference between the 2000 trials of the
  # article and the 4000 simulated here, plus half a unit of the last digit.
  factor_se <- sqrt(1 / 2000 + 1 / 4000)
  for (scenario in c(1, 4, 8)) {
    published <- zhou2019_table3(scenario)
    oc <- do.call(sim_uboin, uboin_oc_args(scenario, n_trials = 4000,
                                           keep_trials = TRUE, seed = 2019))
    info <- paste("Scenario", scenario)
    pub <- c(published$sel, published$stop)
    ours <- c(unname(oc$sel_percent), oc$percent_no_obd)
    p <- (pub + ours) / 200
    allowance <- 400 * sqrt(p * (1 - p)) * factor_se + 0.05
    expect_true(all(abs(pub - ours) <= allowance), info = info)
    allowance <- 4 * apply(oc$trials$n_pts, 2L, stats::sd) * factor_se + 0.05
    expect_true(all(abs(published$n_pts - unname(oc$n_pts_dose)) <=
                      allowance), info = info)
  }
})

test_that("the results of sim_uboin enter oc_results", {
  oc <- do.call(sim_uboin, uboin_oc_args(1, n_trials = 200, seed = 3))
  res <- oc_results(oc)
  expect_true("U-BOIN" %in% as.character(res$design))
  expect_equal(unique(as.character(res$target_type)), "OBD")
  sel <- res[res$metric == "sel_percent", ]
  expect_identical(sel$dose, 1:5)
  expect_identical(sel$value, unname(oc$sel_percent))
  no_sel <- res[res$metric == "percent_no_selection", ]
  expect_identical(no_sel$value, oc$percent_no_obd)
})

test_that("sim_uboin validates its arguments", {
  expect_error(do.call(sim_uboin, uboin_oc_args(1, keep_trials = NA)),
               "'keep_trials' must be TRUE or FALSE")
  expect_error(do.call(sim_uboin, uboin_oc_args(1, resp_lower = 1.2)),
               "'resp_lower' must lie strictly between 0 and 1")
})
