tite_scenarios_example <- function() {
  list(
    list(name = "MTD at dose 3", p_true = c(0.05, 0.15, 0.30, 0.45, 0.60)),
    list(name = "MTD at dose 1", p_true = c(0.30, 0.45, 0.55, 0.65, 0.75)),
    list(name = "All safe",      p_true = c(0.02, 0.04, 0.06, 0.08, 0.10))
  )
}

test_that("sim_tite_boin_multi runs every scenario", {
  oc <- sim_tite_boin_multi(
    target = 0.30, scenarios = tite_scenarios_example(),
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 200, seed = 123
  )

  expect_s3_class(oc, "tite_boin_oc_multi")
  expect_s3_class(oc, "boin_oc_multi")
  expect_length(oc$results, 3L)
  expect_equal(oc$scenario_names,
               c("MTD at dose 3", "MTD at dose 1", "All safe"))
  expect_equal(oc$n_doses, 5L)
  expect_true(all(vapply(oc$results, inherits, logical(1), "tite_boin_oc")))
})

test_that("the summary table has one block of six rows per scenario", {
  oc <- sim_tite_boin_multi(
    target = 0.30, scenarios = tite_scenarios_example(),
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    n_trials = 200, seed = 123
  )
  tab <- oc$summary_table

  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 18L)
  expect_equal(names(tab), c("Scenario", "Item", paste0("DL", 1:5),
                             "Total / No MTD"))
  expect_equal(tab$Scenario[c(1, 7, 13)], oc$scenario_names)
  expect_equal(unique(tab$Item),
               c("True DLT rate (%)", "MTD selected (%)",
                 "Patients treated", "Patients with DLT",
                 "Trial duration (mean)", "Accrual suspended (% trials)"))

  # The added rows carry the scenario values in the last column only.
  for (i in seq_along(oc$results)) {
    res <- oc$results[[i]]
    rows <- 6L * (i - 1L) + 5:6
    expect_equal(tab[["Total / No MTD"]][rows],
                 round(c(res$duration_mean, res$pct_trials_suspended), 1))
    expect_true(all(is.na(as.matrix(tab[rows, paste0("DL", 1:5)]))))
  }

  # The first four rows of each block are those of the BOIN table.
  base <- oc_multi_table(oc$results, oc$scenario_names, oc$n_doses)
  keep <- rep(c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 3)
  expect_equal(unname(as.list(tab[keep, ])), unname(as.list(base)))
})

test_that("a scenario matches the same scenario run on its own", {
  scenarios <- tite_scenarios_example()
  multi <- sim_tite_boin_multi(
    target = 0.30, scenarios = scenarios,
    n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
    method = "ess", n_trials = 300, seed = 42
  )

  for (i in seq_along(scenarios)) {
    single <- sim_tite_boin(
      target = 0.30, p_true = scenarios[[i]]$p_true,
      n_cohort = 10, cohort_size = 3, window = 3, accrual_rate = 2,
      method = "ess", n_trials = 300, seed = 42
    )
    expect_equal(multi$results[[i]]$sel_percent, single$sel_percent)
    expect_equal(multi$results[[i]]$n_pts_dose, single$n_pts_dose)
    expect_equal(multi$results[[i]]$duration_mean, single$duration_mean)
  }
})
