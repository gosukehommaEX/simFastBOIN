gboin_scenarios <- function() {
  list("Scenario 1 of Mu" = mu_scenario(1), "Scenario 6 of Mu" = mu_scenario(6))
}

test_that("sim_gboin_multi runs every scenario with the same design", {
  res <- sim_gboin_multi(
    target = ets_target(), scenarios = gboin_scenarios(),
    score = ets_score(), n_cohort = 10, cohort_size = 3, n_trials = 200,
    seed = 17
  )

  expect_s3_class(res, "boin_oc_multi")
  expect_length(res$results, 2L)
  expect_identical(res$scenario_names,
                   c("Scenario 1 of Mu", "Scenario 6 of Mu"))
  expect_identical(res$n_doses, 6L)
  expect_equal(nrow(res$summary_table), 8L)
  expect_identical(res$summary_table$Item[c(1, 4)],
                   c("True score (%)", "Toxicity score"))

  for (i in 1:2) {
    single <- sim_gboin(
      target = ets_target(), p_grade = gboin_scenarios()[[i]],
      score = ets_score(), n_cohort = 10, cohort_size = 3, n_trials = 200,
      seed = 17
    )
    expect_identical(res$results[[i]]$sel_percent, single$sel_percent)
    expect_identical(res$results[[i]]$n_tox_dose, single$n_tox_dose)
    expect_s3_class(res$results[[i]], "gboin_oc")
  }
})

test_that("scenarios can be given with their names and options are passed", {
  res <- sim_gboin_multi(
    target = ets_target(),
    scenarios = list(list(name = "A", p_grade = mu_scenario(1)),
                     mu_scenario(6)),
    score = ets_score(), n_cohort = 10, cohort_size = 3, n_trials = 100,
    n_earlystop = 100, elim_rule = "count", keep_trials = TRUE, seed = 2
  )

  expect_identical(res$scenario_names, c("A", "Scenario 2"))
  expect_identical(res$results[[2]]$settings$elim_rule, "count")
  expect_false(is.null(res$results[[1]]$trials))
  reported <- capture_messages(
    sim_gboin_multi(ets_target(), gboin_scenarios(), ets_score(), 10, 3,
                    n_trials = 20, verbose = TRUE)
  )
  expect_true(any(grepl("Scenario 2 of 2", reported)))
})

test_that("the result prints as gBOIN and rejects inconsistent scenarios", {
  res <- sim_gboin_multi(
    target = ets_target(), scenarios = gboin_scenarios(),
    score = ets_score(), n_cohort = 10, cohort_size = 3, n_trials = 100,
    seed = 2
  )

  expect_output(print(res), "gBOIN operating characteristics across 2")
  expect_output(print(res), "target score    : ")
  expect_output(print(res), "True score \\(%\\)")
  expect_identical(levels(oc_results(res)$design), "gBOIN")

  expect_error(
    sim_gboin_multi(ets_target(), list(mu_scenario(1), mu_scenario(1)[, 1:5]),
                    ets_score(), 5, 3, n_trials = 10),
    "same numbers"
  )
  expect_error(
    sim_gboin_multi(ets_target(), list(c(0.5, 0.5)), ets_score(), 5, 3),
    "matrix of grade probabilities"
  )
  expect_error(
    sim_gboin_multi(ets_target(), list(A = mu_scenario(1), A = mu_scenario(6)),
                    ets_score(), 5, 3, n_trials = 10),
    "unique"
  )
  expect_error(sim_gboin_multi(ets_target(), list(), ets_score(), 5, 3),
               "non-empty")
})
