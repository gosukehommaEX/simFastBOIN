# The strongest evidence the package can offer: the operating characteristics
# must equal those of BOIN::get.oc() exactly, not within a tolerance, because
# both implementations consume the same random variates in the same order.
#
# The configurations are split across two tests so that neither runs for long.

compare_with_reference <- function(config, n_trials = 500, seed = 6) {
  reference <- BOIN::get.oc(
    target = config$target, p.true = config$p_true,
    ncohort = config$n_cohort, cohortsize = config$cohort_size,
    n.earlystop = config$n_earlystop, startdose = config$start_dose,
    titration = config$titration, cutoff.eli = 0.95,
    extrasafe = config$extrasafe, offset = 0.05,
    boundMTD = config$bound_mtd, ntrial = n_trials, seed = seed
  )

  ours <- sim_boin(
    target = config$target, p_true = config$p_true,
    n_cohort = config$n_cohort, cohort_size = config$cohort_size,
    n_trials = n_trials, start_dose = config$start_dose,
    n_earlystop = config$n_earlystop, cutoff_eli = 0.95,
    extrasafe = config$extrasafe, offset = 0.05,
    titration = config$titration, bound_mtd = config$bound_mtd,
    seed = seed
  )

  info <- paste0("configuration: ", config$name)
  expect_equal(unname(ours$sel_percent), reference$selpercent,
               tolerance = 1e-10, info = info)
  expect_equal(unname(ours$n_pts_dose), reference$npatients,
               tolerance = 1e-10, info = info)
  expect_equal(unname(ours$n_tox_dose), reference$ntox,
               tolerance = 1e-10, info = info)
  expect_equal(ours$percent_no_mtd, reference$percentstop,
               tolerance = 1e-10, info = info)
  expect_equal(ours$total_n_pts, reference$totaln,
               tolerance = 1e-10, info = info)
  expect_equal(ours$total_n_tox, reference$totaltox,
               tolerance = 1e-10, info = info)
  invisible(NULL)
}

test_that("sim_boin reproduces BOIN::get.oc for the core design options", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  configs <- list(
    list(name = "plain",
         target = 0.30, p_true = c(0.05, 0.15, 0.25, 0.45, 0.60),
         n_cohort = 20, cohort_size = 3, n_earlystop = 100, start_dose = 1,
         titration = FALSE, extrasafe = FALSE, bound_mtd = FALSE),
    list(name = "early stopping",
         target = 0.30, p_true = c(0.05, 0.15, 0.25, 0.45, 0.60),
         n_cohort = 20, cohort_size = 3, n_earlystop = 18, start_dose = 1,
         titration = FALSE, extrasafe = FALSE, bound_mtd = FALSE),
    list(name = "titration",
         target = 0.30, p_true = c(0.05, 0.15, 0.25, 0.45, 0.60),
         n_cohort = 20, cohort_size = 3, n_earlystop = 18, start_dose = 1,
         titration = TRUE, extrasafe = FALSE, bound_mtd = FALSE),
    list(name = "extrasafe",
         target = 0.30, p_true = c(0.35, 0.45, 0.55, 0.65, 0.75),
         n_cohort = 20, cohort_size = 3, n_earlystop = 18, start_dose = 1,
         titration = FALSE, extrasafe = TRUE, bound_mtd = FALSE),
    list(name = "bounded MTD",
         target = 0.30, p_true = c(0.10, 0.20, 0.28, 0.36, 0.50),
         n_cohort = 20, cohort_size = 3, n_earlystop = 18, start_dose = 1,
         titration = FALSE, extrasafe = FALSE, bound_mtd = TRUE),
    list(name = "every option",
         target = 0.30, p_true = c(0.15, 0.32, 0.45, 0.60, 0.75),
         n_cohort = 15, cohort_size = 3, n_earlystop = 12, start_dose = 1,
         titration = TRUE, extrasafe = TRUE, bound_mtd = TRUE)
  )

  for (config in configs) compare_with_reference(config)
})

test_that("sim_boin reproduces BOIN::get.oc for other cohort and dose settings", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  configs <- list(
    list(name = "start at dose 2",
         target = 0.25, p_true = c(0.05, 0.10, 0.26, 0.40, 0.55, 0.70),
         n_cohort = 12, cohort_size = 3, n_earlystop = 18, start_dose = 2,
         titration = FALSE, extrasafe = FALSE, bound_mtd = FALSE),
    list(name = "cohort size one",
         target = 0.25, p_true = c(0.05, 0.10, 0.26, 0.40, 0.55),
         n_cohort = 24, cohort_size = 1, n_earlystop = 12, start_dose = 1,
         titration = TRUE, extrasafe = FALSE, bound_mtd = FALSE),
    list(name = "cohort size two",
         target = 0.20, p_true = c(0.02, 0.06, 0.12, 0.21, 0.35),
         n_cohort = 15, cohort_size = 2, n_earlystop = 10, start_dose = 1,
         titration = TRUE, extrasafe = TRUE, bound_mtd = FALSE),
    list(name = "high target",
         target = 0.40, p_true = c(0.45, 0.55, 0.65, 0.75, 0.85),
         n_cohort = 10, cohort_size = 3, n_earlystop = 9, start_dose = 1,
         titration = FALSE, extrasafe = TRUE, bound_mtd = FALSE),
    list(name = "all doses safe",
         target = 0.30, p_true = c(0.02, 0.04, 0.06, 0.08, 0.10),
         n_cohort = 20, cohort_size = 3, n_earlystop = 18, start_dose = 1,
         titration = TRUE, extrasafe = FALSE, bound_mtd = FALSE)
  )

  for (config in configs) compare_with_reference(config)
})
