# The operating characteristics of the drug-combination design must equal
# those of BOIN::get.oc.comb() exactly, not within a tolerance, because both
# implementations consume the same random variates in the same order and call
# the same bivariate isotonic regression. The reference rounds its matrices to
# two decimals and its totals to one, so the comparison rounds in the same way.
#
# The configurations are split across three tests so that none runs for long.
# The reference requires no more rows than columns.

compare_comb_with_reference <- function(config, n_trials = 500, seed = 6) {
  reference <- BOIN::get.oc.comb(
    target = config$target, p.true = config$p_true,
    ncohort = config$n_cohort, cohortsize = config$cohort_size,
    n.earlystop = config$n_earlystop, startdose = config$start_dose,
    titration = config$titration, cutoff.eli = 0.95,
    extrasafe = config$extrasafe, offset = 0.05, ntrial = n_trials,
    mtd.contour = FALSE, boundMTD = config$bound_mtd, seed = seed
  )

  ours <- sim_comb_boin(
    target = config$target, p_true = config$p_true,
    n_cohort = config$n_cohort, cohort_size = config$cohort_size,
    n_trials = n_trials, start_dose = config$start_dose,
    n_earlystop = config$n_earlystop, cutoff_eli = 0.95,
    extrasafe = config$extrasafe, offset = 0.05,
    titration = config$titration, bound_mtd = config$bound_mtd,
    seed = seed
  )

  info <- paste0("configuration: ", config$name)
  expect_equal(round(unname(ours$sel_percent), 2), unname(reference$selpercent),
               info = info)
  expect_equal(round(unname(ours$n_pts_dose), 2), unname(reference$npatients),
               info = info)
  expect_equal(round(unname(ours$n_tox_dose), 2), unname(reference$ntox),
               info = info)
  expect_equal(ours$percent_no_mtd, reference$percentstop,
               tolerance = 1e-8, info = info)
  expect_equal(round(ours$total_n_pts, 1), reference$totaln, info = info)
  expect_equal(round(ours$total_n_tox, 1), reference$totaltox, info = info)
  expect_equal(paste0(round(ours$pct_true_mtd_selected, 1), "%"),
               reference$pcs, info = info)
  expect_equal(paste0(round(ours$pct_patients_true_mtd, 1), "%"),
               reference$npercent, info = info)
  invisible(NULL)
}

comb_reference_matrix <- function(name) {
  switch(name,
    jss = matrix(c(0.02, 0.04, 0.08, 0.14,
                   0.08, 0.25, 0.42, 0.48,
                   0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE),
    # Scenario 1 of Lin and Yin (2017): rows are the levels of agent 2.
    lin_yin_1 = matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
                         0.10, 0.15, 0.30, 0.45, 0.55,
                         0.15, 0.30, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE),
    toxic = matrix(c(0.30, 0.45, 0.60, 0.70, 0.80,
                     0.45, 0.55, 0.65, 0.75, 0.85,
                     0.50, 0.60, 0.70, 0.80, 0.90), nrow = 3, byrow = TRUE),
    safe = matrix(c(0.01, 0.02, 0.03, 0.04,
                    0.02, 0.03, 0.04, 0.05,
                    0.03, 0.04, 0.05, 0.06), nrow = 3, byrow = TRUE),
    small = matrix(c(0.05, 0.20, 0.35,
                     0.15, 0.30, 0.50), nrow = 2, byrow = TRUE)
  )
}

test_that("sim_comb_boin reproduces BOIN::get.oc.comb for the core options", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  configs <- list(
    list(name = "plain", p_true = comb_reference_matrix("jss"),
         target = 0.25, n_cohort = 16, cohort_size = 3, n_earlystop = 100,
         start_dose = c(1, 1), titration = FALSE, extrasafe = FALSE,
         bound_mtd = FALSE),
    list(name = "early stopping", p_true = comb_reference_matrix("lin_yin_1"),
         target = 0.30, n_cohort = 20, cohort_size = 3, n_earlystop = 12,
         start_dose = c(1, 1), titration = FALSE, extrasafe = FALSE,
         bound_mtd = FALSE),
    list(name = "titration", p_true = comb_reference_matrix("jss"),
         target = 0.25, n_cohort = 12, cohort_size = 3, n_earlystop = 18,
         start_dose = c(1, 1), titration = TRUE, extrasafe = FALSE,
         bound_mtd = FALSE)
  )

  for (config in configs) compare_comb_with_reference(config)
})

test_that("sim_comb_boin reproduces BOIN::get.oc.comb for the safety options", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  configs <- list(
    list(name = "extrasafe", p_true = comb_reference_matrix("toxic"),
         target = 0.30, n_cohort = 20, cohort_size = 3, n_earlystop = 18,
         start_dose = c(1, 1), titration = FALSE, extrasafe = TRUE,
         bound_mtd = FALSE),
    list(name = "bounded MTD", p_true = comb_reference_matrix("lin_yin_1"),
         target = 0.30, n_cohort = 20, cohort_size = 3, n_earlystop = 18,
         start_dose = c(1, 1), titration = FALSE, extrasafe = FALSE,
         bound_mtd = TRUE),
    list(name = "every option", p_true = comb_reference_matrix("lin_yin_1"),
         target = 0.30, n_cohort = 15, cohort_size = 3, n_earlystop = 12,
         start_dose = c(2, 2), titration = TRUE, extrasafe = TRUE,
         bound_mtd = TRUE)
  )

  for (config in configs) compare_comb_with_reference(config)
})

test_that("sim_comb_boin reproduces BOIN::get.oc.comb for other cohort and dose settings", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  configs <- list(
    list(name = "cohort size one", p_true = comb_reference_matrix("small"),
         target = 0.30, n_cohort = 30, cohort_size = 1, n_earlystop = 12,
         start_dose = c(1, 1), titration = TRUE, extrasafe = FALSE,
         bound_mtd = FALSE),
    list(name = "cohort size two", p_true = comb_reference_matrix("jss"),
         target = 0.25, n_cohort = 20, cohort_size = 2, n_earlystop = 9,
         start_dose = c(1, 1), titration = TRUE, extrasafe = TRUE,
         bound_mtd = FALSE),
    list(name = "all combinations safe", p_true = comb_reference_matrix("safe"),
         target = 0.30, n_cohort = 20, cohort_size = 3, n_earlystop = 18,
         start_dose = c(1, 1), titration = TRUE, extrasafe = FALSE,
         bound_mtd = TRUE)
  )

  for (config in configs) compare_comb_with_reference(config)
})
