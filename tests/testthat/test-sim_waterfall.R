# The average number of patients as BOIN::get.oc.comb() prints it, with two
# decimals, and then rounded half up to the one decimal of the article.
waterfall_two_decimals <- function(x) {
  matrix(as.numeric(formatC(x, digits = 2, format = "f")), nrow = nrow(x))
}
waterfall_one_decimal <- function(x) {
  floor((round(waterfall_two_decimals(x) * 100) + 5) / 10) / 10
}

test_that("sim_waterfall returns a well formed object", {
  oc <- sim_waterfall(
    target = 0.30, p_true = zhang2016_p_true(12), cohort_size = 3,
    n_trials = 300, keep_trials = TRUE, seed = 1
  )

  expect_s3_class(oc, "waterfall_oc")
  expect_equal(dim(oc$sel_percent), c(3L, 5L))
  expect_equal(dimnames(oc$sel_percent), list(c("A1", "A2", "A3"),
                                              paste0("B", 1:5)))
  expect_equal(dim(oc$n_pts_dose), c(3L, 5L))
  expect_equal(sum(oc$n_pts_dose), oc$total_n_pts)
  expect_equal(sum(oc$n_tox_dose), oc$total_n_tox)
  expect_equal(oc$total_n_enrolled, oc$total_n_pts)
  # Each level of drug A selects at most one combination per trial.
  expect_true(all(rowSums(oc$sel_percent) <= 100 + 1e-8))
  # The third level of drug A has no true MTD: its lowest DLT rate is 0.49.
  expect_equal(unname(oc$true_mtd), cbind(1:2, c(3L, 1L)))
  expect_equal(colnames(oc$true_mtd), c("dose_a", "dose_b"))
  expect_equal(oc$pct_patients_contour + oc$pct_patients_above +
                 oc$pct_patients_below, 100)
  expect_equal(sum(oc$stop_reason_percent), 100)
  expect_equal(dim(oc$trials$mtd), c(300L, 3L))
  expect_equal(oc$settings$n_cohort, c(10L, 6L, 6L))

  # The summaries follow from the trial level data.
  mtd <- oc$trials$mtd
  expect_equal(oc$percent_no_mtd, mean(rowSums(!is.na(mtd)) == 0) * 100)
  expect_true(all(is.na(mtd) | (mtd >= 1L & mtd <= 5L)))
  expect_equal(unname(oc$sel_percent[1L, 3L]),
               mean(!is.na(mtd[, 1L]) & mtd[, 1L] == 3L) * 100)
  correct <- !is.na(mtd[, 1L]) & mtd[, 1L] == 3L &
    !is.na(mtd[, 2L]) & mtd[, 2L] == 1L & is.na(mtd[, 3L])
  expect_equal(oc$pct_contour_selected, mean(correct) * 100)

  plain <- sim_waterfall(
    target = 0.30, p_true = zhang2016_p_true(12), cohort_size = 3,
    n_trials = 300, seed = 1
  )
  expect_null(plain$trials)
  expect_equal(plain$sel_percent, oc$sel_percent)
})

test_that("sim_waterfall reproduces Zhang and Yuan (2016)", {
  # Tables IV, V and S1 of the article, which come from
  # BOIN::get.oc.comb(mtd.contour = TRUE) with its default seed of 6 and
  # 1000 trials: the cells are equal, not within a tolerance.
  for (scenario in c(1, 4, 12)) {
    published <- zhang2016_published(scenario)
    oc <- sim_waterfall(
      target = 0.30, p_true = zhang2016_p_true(scenario), cohort_size = 3,
      n_trials = 1000, boin_compatible = TRUE, seed = 6
    )
    info <- paste("Scenario", scenario)
    expect_equal(round(unname(oc$sel_percent), 1), published$sel, info = info)
    expect_equal(waterfall_one_decimal(oc$n_pts_dose), published$n_pts,
                 info = info)
    expect_equal(round(oc$pct_contour_selected, 1), published$pcs,
                 info = info)
    expect_equal(round(oc$pct_patients_above, 1), published$above,
                 info = info)
    expect_equal(round(oc$pct_patients_contour, 1), published$at,
                 info = info)
  }
})

compare_waterfall_with_reference <- function(config, n_trials = 200,
                                             seed = 6) {
  cohorts <- waterfall_simulate(
    target = config$target, p_true = config$p_true,
    cohort_size = config$cohort_size, n_trials = 1, seed = 1
  )$settings$n_cohort

  reference <- BOIN::get.oc.comb(
    target = config$target, p.true = config$p_true, ncohort = cohorts,
    cohortsize = config$cohort_size, n.earlystop = config$n_earlystop,
    titration = config$titration, cutoff.eli = 0.95,
    extrasafe = config$extrasafe, offset = 0.05, ntrial = n_trials,
    mtd.contour = TRUE, boundMTD = config$bound_mtd, seed = seed
  )

  ours <- sim_waterfall(
    target = config$target, p_true = config$p_true,
    cohort_size = config$cohort_size, n_cohort = cohorts,
    n_trials = n_trials, n_earlystop = config$n_earlystop,
    cutoff_eli = 0.95, extrasafe = config$extrasafe, offset = 0.05,
    titration = config$titration, bound_mtd = config$bound_mtd,
    boin_compatible = TRUE, seed = seed
  )

  one <- function(x) as.numeric(formatC(x, digits = 1, format = "f"))
  info <- paste0("configuration: ", config$name)
  expect_equal(round(unname(ours$sel_percent), 2), unname(reference$selpercent),
               info = info)
  expect_equal(waterfall_two_decimals(ours$n_pts_dose),
               unname(reference$npatients), info = info)
  expect_equal(waterfall_two_decimals(ours$n_tox_dose),
               unname(reference$ntox), info = info)
  expect_equal(one(sum(ours$n_pts_dose)), reference$totaln, info = info)
  expect_equal(one(sum(ours$n_tox_dose)), reference$totaltox, info = info)
  expect_equal(paste0(round(ours$pct_contour_selected, 1), "%"),
               reference$pcs.contour, info = info)
  expect_equal(paste0(round(ours$pct_patients_contour, 1), "%"),
               reference$npercent.contour, info = info)
  expect_equal(paste0(round(ours$pct_patients_above, 1), "%"),
               reference$npercent.above.contour, info = info)
  expect_equal(paste0(round(ours$pct_patients_below, 1), "%"),
               reference$npercent.below.contour, info = info)
  invisible(NULL)
}

test_that("sim_waterfall reproduces BOIN::get.oc.comb for the core options", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  configs <- list(
    list(name = "Scenario 4", p_true = zhang2016_p_true(4), target = 0.30,
         cohort_size = 3, n_earlystop = 12, titration = FALSE,
         extrasafe = FALSE, bound_mtd = FALSE),
    list(name = "Scenario 12, cohorts of 2", p_true = zhang2016_p_true(12),
         target = 0.30, cohort_size = 2, n_earlystop = 9, titration = FALSE,
         extrasafe = FALSE, bound_mtd = FALSE)
  )
  for (config in configs) compare_waterfall_with_reference(config)
})

test_that("sim_waterfall reproduces BOIN::get.oc.comb for the safety options", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  configs <- list(
    list(name = "steep, extrasafe", p_true = waterfall_steep_p_true(),
         target = 0.30, cohort_size = 3, n_earlystop = 12, titration = FALSE,
         extrasafe = TRUE, bound_mtd = FALSE),
    list(name = "Scenario 4, bound_mtd", p_true = zhang2016_p_true(4),
         target = 0.30, cohort_size = 3, n_earlystop = 12, titration = FALSE,
         extrasafe = FALSE, bound_mtd = TRUE)
  )
  for (config in configs) compare_waterfall_with_reference(config)
})

test_that("sim_waterfall reproduces BOIN::get.oc.comb with titration", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  configs <- list(
    list(name = "Scenario 1, titration", p_true = zhang2016_p_true(1),
         target = 0.30, cohort_size = 3, n_earlystop = 12, titration = TRUE,
         extrasafe = FALSE, bound_mtd = FALSE),
    list(name = "steep, all options", p_true = waterfall_steep_p_true(),
         target = 0.30, cohort_size = 2, n_earlystop = 12, titration = TRUE,
         extrasafe = TRUE, bound_mtd = TRUE)
  )
  for (config in configs) compare_waterfall_with_reference(config)
})

test_that("the default rules keep the patients of stopped subtrials", {
  args <- list(target = 0.30, p_true = zhang2016_p_true(4), cohort_size = 3,
               n_trials = 500, seed = 6)
  kept <- do.call(sim_waterfall, args)
  dropped <- do.call(sim_waterfall, c(args, boin_compatible = TRUE))

  expect_equal(kept$total_n_enrolled, dropped$total_n_enrolled)
  expect_equal(kept$total_n_pts, kept$total_n_enrolled)
  expect_lt(dropped$total_n_pts, dropped$total_n_enrolled)
  expect_equal(kept$stop_reason_percent, dropped$stop_reason_percent)
})

test_that("a matrix without a true MTD has no contour summaries", {
  p_true <- matrix(c(0.50, 0.60, 0.70,
                     0.60, 0.70, 0.80), nrow = 2, byrow = TRUE)
  oc <- sim_waterfall(target = 0.30, p_true = p_true, cohort_size = 3,
                      n_trials = 100, seed = 1)

  expect_equal(nrow(oc$true_mtd), 0L)
  expect_true(is.na(oc$pct_contour_selected))
  expect_true(is.na(oc$pct_patients_contour))
  expect_true(is.na(oc$pct_patients_above))
  expect_true(is.na(oc$pct_patients_below))
  expect_gt(oc$percent_no_mtd, 50)
})

test_that("the regions of the contour follow the reference implementation", {
  # Scenario 1: the true MTDs are A1B3 and A2B2. A1B1, A1B2 and A2B1 are
  # below the contour and A2B3 is above it.
  contour <- waterfall_contour(zhang2016_p_true(1), 0.30)
  expect_equal(unname(contour$true_mtd), cbind(1:2, c(3L, 2L)))
  expect_equal(contour$below, matrix(c(TRUE, TRUE, FALSE,
                                       TRUE, FALSE, FALSE),
                                     nrow = 2, byrow = TRUE))
  expect_equal(contour$above, matrix(c(FALSE, FALSE, FALSE,
                                       FALSE, FALSE, TRUE),
                                     nrow = 2, byrow = TRUE))
  expect_equal(contour$at, matrix(c(FALSE, FALSE, TRUE,
                                    FALSE, TRUE, FALSE),
                                  nrow = 2, byrow = TRUE))

  # Scenario 4: only A1B1 is a true MTD; every other combination is above.
  contour <- waterfall_contour(zhang2016_p_true(4), 0.30)
  expect_equal(unname(contour$true_mtd), cbind(1L, 1L))
  expect_equal(sum(contour$above), 5L)
  expect_equal(sum(contour$below), 0L)
})
