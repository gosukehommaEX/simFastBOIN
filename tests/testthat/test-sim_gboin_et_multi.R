takeda_multi <- function(...) {
  sim_gboin_et_multi(
    target = 0.47 / 1.5, target_resp = ees_target(),
    scenarios = list("Scenario 1" = takeda_scenario(1),
                     "Scenario 5" = takeda_scenario(5)),
    score = c(0, 0.5, 1, 1.5), score_resp = ees_score(), n_cohort = 10,
    cohort_size = 3, ...
  )
}

test_that("sim_gboin_et_multi collects the results of sim_gboin_et", {
  res <- takeda_multi(n_trials = 200, seed = 9, est_method = "fp_linear")

  expect_s3_class(res, "boin_et_oc_multi")
  expect_identical(res$scenario_names, c("Scenario 1", "Scenario 5"))
  expect_identical(res$n_doses, 6L)
  for (s in c(1, 5)) {
    single <- takeda_oc(s, n_trials = 200, seed = 9,
                        est_method = "fp_linear")
    multi <- res$results[[paste("Scenario", s)]]
    expect_equal(multi$sel_percent, single$sel_percent)
    expect_equal(multi$n_pts_dose, single$n_pts_dose)
    expect_equal(multi$n_resp_dose, single$n_resp_dose)
  }
})

test_that("the summary table has the labels of scores", {
  res <- takeda_multi(n_trials = 50, seed = 9)
  tab <- res$summary_table

  expect_identical(nrow(tab), 12L)
  expect_identical(tab$Item[1:6],
                   c("True toxicity score (%)", "True efficacy score (%)",
                     "OBD selected (%)", "Patients treated",
                     "Toxicity score", "Efficacy score"))
  first <- res$results[[1]]
  expect_equal(as.numeric(tab[3, paste0("DL", 1:6)]),
               unname(round(first$sel_percent, 1)))
  out <- capture.output(print(res))
  expect_true(any(grepl("gBOIN-ET operating characteristics across 2",
                        out, fixed = TRUE)))
  expect_true(any(grepl("target efficacy score", out, fixed = TRUE)))
})

test_that("scenario names come from the list or from the elements", {
  sc <- takeda_scenario(1)
  res <- sim_gboin_et_multi(
    0.47 / 1.5, ees_target(),
    scenarios = list(sc, c(sc, list(name = "Named"))),
    score = c(0, 0.5, 1, 1.5), score_resp = ees_score(), n_cohort = 5,
    cohort_size = 3, n_trials = 20
  )
  expect_identical(res$scenario_names, c("Scenario 1", "Named"))
})

test_that("sim_gboin_et_multi rejects invalid scenarios", {
  sc <- takeda_scenario(1)
  run <- function(scenarios) {
    sim_gboin_et_multi(0.47 / 1.5, ees_target(), scenarios,
                       c(0, 0.5, 1, 1.5), ees_score(), 5, 3, n_trials = 10)
  }
  expect_error(run(list()), "non-empty")
  expect_error(run(list(sc$p_grade)), "p_resp_grade")
  short <- list(p_grade = sc$p_grade[, 1:5],
                p_resp_grade = sc$p_resp_grade[, 1:5])
  expect_error(run(list(sc, short)), "same numbers")
  expect_error(run(list(A = sc, A = sc)), "unique")
})
