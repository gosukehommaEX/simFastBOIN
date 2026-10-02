backfill_scenarios_example <- function() {
  list(
    list(name = "Rising", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
         p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58)),
    list(name = "Plateau", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
         p_resp = c(0.30, 0.32, 0.35, 0.36, 0.36)),
    list(p_true = c(0.02, 0.05, 0.08, 0.11, 0.25),
         p_resp = c(0.05, 0.10, 0.15, 0.20, 0.30))
  )
}

test_that("a scenario matches the same scenario run on its own", {
  scenarios <- backfill_scenarios_example()
  multi <- sim_bf_boin_multi(
    target = 0.25, scenarios = scenarios, n_cohort = 10, cohort_size = 3,
    window = 1, accrual_rate = 3, n_trials = 300, seed = 42
  )

  expect_s3_class(multi, "bf_boin_oc_multi")
  expect_s3_class(multi, "backfill_oc_multi")
  expect_identical(multi$scenario_names, c("Rising", "Plateau", "Scenario 3"))
  for (i in seq_along(scenarios)) {
    single <- sim_bf_boin(
      target = 0.25, p_true = scenarios[[i]]$p_true,
      p_resp = scenarios[[i]]$p_resp, n_cohort = 10, cohort_size = 3,
      window = 1, accrual_rate = 3, n_trials = 300, seed = 42
    )
    expect_equal(multi$results[[i]]$sel_percent, single$sel_percent)
    expect_equal(multi$results[[i]]$n_bf_dose, single$n_bf_dose)
    expect_equal(multi$results[[i]]$duration_mean, single$duration_mean)
  }
})

test_that("the summary table has seven rows per scenario", {
  multi <- sim_bf_boin_multi(
    target = 0.25, scenarios = backfill_scenarios_example(), n_cohort = 10,
    cohort_size = 3, window = 1, accrual_rate = 3, n_trials = 200, seed = 1
  )
  tab <- multi$summary_table

  expect_equal(nrow(tab), 21L)
  expect_named(tab, c("Scenario", "Item", paste0("DL", 1:5), "Total / No MTD"))
  expect_identical(tab$Item[1:7], c(
    "True DLT rate (%)", "True response rate (%)", "MTD selected (%)",
    "Patients treated", "Backfilled patients", "Patients with DLT",
    "Trial duration (mean)"
  ))
  expect_identical(tab$Scenario[c(1, 8, 15)], c("Rising", "Plateau", "Scenario 3"))
  first <- multi$results[[1L]]
  expect_equal(tab[["Total / No MTD"]][5], round(first$total_n_bf, 1))
  expect_equal(tab[["Total / No MTD"]][7], round(first$duration_mean, 1))
  expect_equal(unlist(tab[2, paste0("DL", 1:5)], use.names = FALSE),
               c(10, 20, 30, 45, 58))
})

test_that("invalid scenarios are rejected", {
  call_with <- function(scenarios) {
    sim_bf_boin_multi(target = 0.25, scenarios = scenarios, n_cohort = 5,
                      cohort_size = 3, window = 1, accrual_rate = 3,
                      n_trials = 10)
  }

  expect_error(call_with(list(c(0.1, 0.25, 0.4))), "p_resp")
  expect_error(call_with(list(list(p_true = c(0.1, 0.25, 0.4)))), "p_resp")
  expect_error(call_with(list(list(p_true = c(0.1, 0.25, 0.4),
                                   p_resp = c(0.1, 0.2)))), "p_resp")
  expect_error(call_with(list()), "non-empty")
})
