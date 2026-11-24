tite_gboin_et_scenarios <- function() {
  list("Scenario 1" = takeda2023_scenario(1),
       "Scenario 4" = takeda2023_scenario(4))
}

# The design of the article; arguments given in ... replace these.
tite_takeda_multi <- function(...) {
  args <- list(
    target = 0.47 / 1.5, target_resp = 1.75 / 3, score = c(0, 0.5, 1, 1.5),
    score_resp = c(0, 0.25, 1, 3), n_cohort = 12, cohort_size = 3,
    window = 30, resp_window = 60, accrual_rate = 1 / 10
  )
  given <- list(...)
  args[names(given)] <- given
  do.call(sim_tite_gboin_et_multi, args)
}

test_that("sim_tite_gboin_et_multi runs every scenario with the same design", {
  scenarios <- tite_gboin_et_scenarios()
  res <- tite_takeda_multi(scenarios = scenarios, n_trials = 100, seed = 3)

  expect_s3_class(res, "boin_et_oc_multi")
  expect_identical(res$scenario_names, names(scenarios))
  expect_identical(res$n_doses, 6L)
  expect_equal(nrow(res$summary_table), 12L)
  one <- tite_takeda_oc(4, n_trials = 100, seed = 3)
  expect_identical(res$results[[2]]$sel_percent, one$sel_percent)
  expect_identical(res$results[[2]]$n_pts_dose, one$n_pts_dose)
  expect_identical(res$results[[2]]$duration, one$duration)
  expect_identical(res$results[[2]]$design, "TITE-gBOIN-ET")
  expect_true("True toxicity score (%)" %in% res$summary_table$Item)
  expect_output(print(res),
                "TITE-gBOIN-ET operating characteristics across 2")
})

test_that("the arguments reach every scenario", {
  scenarios <- tite_gboin_et_scenarios()
  res <- tite_takeda_multi(scenarios = scenarios, n_trials = 50, seed = 3,
                           min_completed_share = 0, resp_window = 90,
                           boinet_compatible = TRUE, keep_trials = TRUE)

  for (r in res$results) {
    expect_identical(r$settings$min_completed_share, 0)
    expect_identical(r$settings$resp_window, 90)
    expect_true(r$settings$boinet_compatible)
    expect_s3_class(r$trials, "boin_et_trials")
  }
})

test_that("scenario names come from the list or from the elements", {
  scenarios <- unname(tite_gboin_et_scenarios())
  scenarios[[2]]$name <- "Plateau"
  res <- tite_takeda_multi(scenarios = scenarios, n_trials = 20, seed = 3)

  expect_identical(res$scenario_names, c("Scenario 1", "Plateau"))
})

test_that("sim_tite_gboin_et_multi rejects invalid scenarios", {
  scenarios <- tite_gboin_et_scenarios()

  expect_error(tite_takeda_multi(scenarios = list()), "non-empty list")
  expect_error(tite_takeda_multi(scenarios = list(list(p_grade = 1))),
               "matrices")
  bad <- scenarios
  bad[[2]]$p_grade <- bad[[2]]$p_grade[, 1:5]
  expect_error(tite_takeda_multi(scenarios = bad), "same numbers")
  names(bad) <- c("A", "A")
  bad[[2]] <- scenarios[[2]]
  expect_error(tite_takeda_multi(scenarios = bad), "unique")
})

test_that("verbose reports each scenario", {
  msgs <- capture_messages(
    tite_takeda_multi(scenarios = tite_gboin_et_scenarios(), n_trials = 10,
                      verbose = TRUE, seed = 3)
  )
  expect_length(grep("^Scenario", msgs), 2L)
})
