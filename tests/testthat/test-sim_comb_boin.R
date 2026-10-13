jss_p_true <- function() {
  matrix(c(0.02, 0.04, 0.08, 0.14,
           0.08, 0.25, 0.42, 0.48,
           0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE)
}

test_that("sim_comb_boin reproduces the published example of Yan et al. (2020)", {
  # Section 4 of Yan et al. (2020, Journal of Statistical Software 94(13)):
  # get.oc.comb(target = 0.25, p.true, ncohort = 16, cohortsize = 3,
  # ntrial = 1000), with the default n.earlystop = 100 and seed = 6, and the
  # values printed by summary(oc.comb).
  oc <- sim_comb_boin(
    target = 0.25, p_true = jss_p_true(), n_cohort = 16, cohort_size = 3,
    n_earlystop = 100, n_trials = 1000, seed = 6
  )

  # The selection percentages are multiples of 0.1 and must match exactly.
  expect_equal(round(unname(oc$sel_percent), 2),
               matrix(c(0.00, 1.10, 2.80, 19.80,
                        4.30, 37.60, 7.10, 1.40,
                        21.80, 3.90, 0.20, 0.00), nrow = 3, byrow = TRUE))

  # The averages are multiples of 0.001, and a published value with d decimals
  # is reproduced when the simulated value lies within half a unit of its last
  # decimal. A value exactly half-way can be printed either way: the average
  # number of patients at A2B4 is exactly 1.905, which the article prints as
  # 1.90 and R 4.6.0 as 1.91 (round() was changed in R 4.0.0).
  published_pts <- matrix(c(4.06, 3.73, 3.22, 4.18,
                            6.08, 10.04, 4.19, 1.90,
                            5.85, 3.65, 0.80, 0.30), nrow = 3, byrow = TRUE)
  published_tox <- matrix(c(0.10, 0.15, 0.25, 0.58,
                            0.47, 2.57, 1.80, 0.89,
                            1.48, 1.65, 0.39, 0.18), nrow = 3, byrow = TRUE)
  expect_lte(max(abs(unname(oc$n_pts_dose) - published_pts)), 0.005 + 1e-9)
  expect_lte(max(abs(unname(oc$n_tox_dose) - published_tox)), 0.005 + 1e-9)
  expect_lte(abs(oc$total_n_tox - 10.5), 0.05 + 1e-9)
  expect_lte(abs(oc$total_n_pts - 48.0), 0.05 + 1e-9)
  expect_equal(oc$n_pts_dose[2, 4], 1.905)

  expect_equal(round(oc$pct_true_mtd_selected, 1), 59.4)
  expect_equal(round(oc$pct_patients_true_mtd, 1), 33.1)
})

test_that("the example of Yan et al. (2020) agrees with a separate implementation", {
  # Sums over the 1,000 trials of the example, from a port of
  # BOIN::get.oc.comb() to Python run on a replica of R's random number stream.
  oc <- sim_comb_boin(
    target = 0.25, p_true = jss_p_true(), n_cohort = 16, cohort_size = 3,
    n_earlystop = 100, n_trials = 1000, seed = 6
  )

  expect_equal(unname(oc$n_pts_dose) * 1000,
               matrix(c(4056, 3735, 3216, 4179,
                        6078, 10041, 4194, 1905,
                        5853, 3648, 798, 297), nrow = 3, byrow = TRUE))
  expect_equal(unname(oc$n_tox_dose) * 1000,
               matrix(c(99, 146, 248, 585,
                        468, 2566, 1796, 886,
                        1478, 1649, 386, 177), nrow = 3, byrow = TRUE))
})

test_that("sim_comb_boin returns a well formed object", {
  oc <- sim_comb_boin(
    target = 0.25, p_true = jss_p_true(), n_cohort = 10, cohort_size = 3,
    n_trials = 300, seed = 1
  )

  expect_s3_class(oc, "comb_boin_oc")
  expect_equal(dim(oc$sel_percent), c(3L, 4L))
  expect_equal(dimnames(oc$sel_percent),
               list(c("A1", "A2", "A3"), c("B1", "B2", "B3", "B4")))
  expect_equal(dim(oc$n_pts_dose), c(3L, 4L))
  expect_equal(dim(oc$n_tox_dose), c(3L, 4L))
  expect_equal(sum(oc$sel_percent) + oc$percent_no_mtd, 100)
  expect_equal(sum(oc$n_pts_dose), oc$total_n_pts)
  expect_equal(sum(oc$n_tox_dose), oc$total_n_tox)
  expect_equal(sum(oc$stop_reason_percent), 100)
  expect_null(oc$trials)
  # A3B1 and A2B2, in column-major order.
  expect_equal(oc$true_mtd,
               matrix(c(3L, 2L, 1L, 2L), ncol = 2,
                      dimnames = list(NULL, c("dose_a", "dose_b"))))
})

test_that("the summaries agree with the trial data", {
  oc <- sim_comb_boin(
    target = 0.25, p_true = jss_p_true(), n_cohort = 10, cohort_size = 3,
    n_trials = 300, keep_trials = TRUE, seed = 2
  )
  tr <- oc$trials

  sel <- table(factor(tr$mtd_a, 1:3), factor(tr$mtd_b, 1:4))
  expect_equal(unname(oc$sel_percent), unname(matrix(sel, 3, 4)) / 300 * 100)
  expect_equal(oc$percent_no_mtd, mean(is.na(tr$mtd_a)) * 100)

  # The true MTD combinations are A2B2 and A3B1.
  at_mtd <- (tr$mtd_a == 2 & tr$mtd_b == 2) | (tr$mtd_a == 3 & tr$mtd_b == 1)
  expect_equal(oc$pct_true_mtd_selected, mean(at_mtd %in% TRUE) * 100)
  pts_mtd <- sum(tr$n_pts[, 2, 2]) + sum(tr$n_pts[, 3, 1])
  expect_equal(oc$pct_patients_true_mtd, pts_mtd / sum(tr$n_pts) * 100)

  # Combinations above the default cutoff, the target 0.25.
  above <- jss_p_true() > 0.25
  n_above <- apply(tr$n_pts, 1L, function(m) sum(m[above]))
  expect_equal(oc$overdose$pct_patients, sum(n_above) / sum(tr$n_pts) * 100)
  expect_equal(oc$overdose$pct_trials_any, mean(n_above > 0) * 100)
  expect_equal(nrow(oc$overdose$combinations), sum(above))

  stopped <- tr$stop_reason %in% c("lowest_dose_eliminated",
                                   "lowest_dose_too_toxic")
  expect_true(all(is.na(tr$mtd_a[stopped])))
})

test_that("a scenario where every combination is toxic mostly selects nothing", {
  p_true <- matrix(c(0.50, 0.60,
                     0.60, 0.70), nrow = 2, byrow = TRUE)
  oc <- sim_comb_boin(
    target = 0.25, p_true = p_true, n_cohort = 10, cohort_size = 3,
    n_trials = 300, extrasafe = TRUE, seed = 3
  )

  expect_true(oc$percent_no_mtd > 50)
  expect_equal(sum(oc$sel_percent) + oc$percent_no_mtd, 100)
  expect_equal(nrow(oc$overdose$combinations), 4L)
})

test_that("verbose reports progress and is silent otherwise", {
  # capture_messages() takes every message, so none escape to the console.
  reported <- capture_messages(
    sim_comb_boin(target = 0.25, p_true = jss_p_true(), n_cohort = 4,
                  cohort_size = 3, n_trials = 10, verbose = TRUE, seed = 4)
  )
  expect_true(any(grepl("Selecting the MTD", reported, fixed = TRUE)))
  expect_silent(
    sim_comb_boin(target = 0.25, p_true = jss_p_true(), n_cohort = 4,
                  cohort_size = 3, n_trials = 10, seed = 4)
  )
})

test_that("an invalid overdose cutoff is rejected before simulating", {
  expect_error(
    sim_comb_boin(target = 0.25, p_true = jss_p_true(), n_cohort = 4,
                  cohort_size = 3, n_trials = 10, overdose_cutoff = 1.5),
    "overdose_cutoff"
  )
})
