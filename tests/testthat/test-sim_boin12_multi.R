scenarios <- list(
  "Plateau" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                   p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
  "Increasing" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                      p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
)

test_that("each scenario matches its own simulation", {
  oc <- sim_boin12_multi(target = 0.35, target_resp = 0.25,
                         scenarios = scenarios, n_cohort = 12,
                         cohort_size = 3, n_trials = 200, seed = 1)
  expect_s3_class(oc, "boin12_oc_multi")
  expect_equal(oc$scenario_names, c("Plateau", "Increasing"))
  expect_equal(nrow(oc$summary_table), 12L)
  one <- sim_boin12(target = 0.35, target_resp = 0.25,
                    p_true = scenarios[[2]]$p_true,
                    p_resp = scenarios[[2]]$p_resp, n_cohort = 12,
                    cohort_size = 3, n_trials = 200, seed = 1)
  expect_equal(oc$results[[2]]$sel_percent, one$sel_percent)
  expect_equal(oc$results[[2]]$n_pts_dose, one$n_pts_dose)
})

test_that("invalid scenarios are rejected", {
  expect_error(sim_boin12_multi(0.35, 0.25, list(list(p_true = 0.1)), 10, 3),
               "p_resp")
  msgs <- capture_messages(
    sim_boin12_multi(0.35, 0.25, scenarios, 4, 3, n_trials = 5,
                     verbose = TRUE)
  )
  expect_identical(sum(grepl("Scenario [12] of 2", msgs)), 2L)
})

test_that("the rules for the next dose reach every scenario", {
  oc <- sim_boin12_multi(target = 0.35, target_resp = 0.25,
                         scenarios = scenarios, n_cohort = 12,
                         cohort_size = 3, n_trials = 100,
                         tie_break = "random", no_candidate = "stop",
                         explore_from = "admissible", seed = 1)
  one <- sim_boin12(target = 0.35, target_resp = 0.25,
                    p_true = scenarios[[2]]$p_true,
                    p_resp = scenarios[[2]]$p_resp, n_cohort = 12,
                    cohort_size = 3, n_trials = 100, tie_break = "random",
                    no_candidate = "stop", explore_from = "admissible",
                    seed = 1)
  expect_equal(oc$results[[2]]$sel_percent, one$sel_percent)
  expect_equal(oc$results[[2]]$n_pts_dose, one$n_pts_dose)
  for (res in oc$results) {
    expect_equal(res$settings$tie_break, "random")
    expect_equal(res$settings$no_candidate, "stop")
    expect_equal(res$settings$explore_from, "admissible")
  }
})
