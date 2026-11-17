p <- c(0.10, 0.25, 0.40, 0.55)

# A result with chosen values, so that every quantity is known by hand.
fake_3p3 <- function(p_true) {
  structure(list(sel_percent = c(10, 50, 30, 5), percent_no_mtd = 5,
                 n_pts_dose = c(6, 9, 4, 1), n_tox_dose = c(1, 2, 1, 1),
                 total_n_pts = 20, total_n_tox = 5, p_true = p_true),
            class = "oc_3p3")
}

test_that("the quantities follow from the per-dose averages", {
  truth <- oc_scenarios(list(A = p, B = p / 2), target = 0.30)
  cmp <- oc_compare(Fake = list(A = fake_3p3(p), B = fake_3p3(p / 2)),
                    truth = truth)
  expect_s3_class(cmp, "oc_compare")
  expect_identical(attr(cmp, "overdose_cutoff"), 0.30)
  a <- cmp[cmp$scenario == "A", ]
  # True MTD dose 2; overdoses (above 0.30) doses 3 and 4.
  expect_equal(a$pct_correct, 50)
  expect_equal(a$n_pts_correct, 9)
  expect_equal(a$pct_patients_correct, 45)
  expect_equal(a$pct_overdose_selected, 35)
  expect_equal(a$n_pts_overdose, 5)
  expect_equal(a$pct_patients_overdose, 25)
  expect_equal(a$pct_no_selection, 5)
  expect_equal(a$total_n_pts, 20)
  expect_true(is.na(a$duration))
  b <- cmp[cmp$scenario == "B", ]
  # True MTD dose 4 (0.275); no overdose.
  expect_equal(b$pct_correct, 5)
  expect_equal(b$pct_overdose_selected, 0)
  expect_equal(b$n_pts_overdose, 0)

  # A lower cutoff makes dose 2 of Scenario A an overdose too.
  low <- oc_compare(Fake = list(A = fake_3p3(p)), truth = truth,
                    overdose_cutoff = 0.2)
  expect_equal(low$pct_overdose_selected, 85)
  expect_identical(levels(low$scenario), "A")
})

test_that("designs selecting an OBD are judged against the true OBD", {
  scen <- list(
    "Scenario 1" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                        p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
    "Scenario 16" = list(p_true = c(0.25, 0.55, 0.67, 0.87, 0.95),
                         p_resp = c(0.01, 0.55, 0.55, 0.55, 0.55))
  )
  truth <- oc_scenarios(scen, target = 0.35, target_resp = 0.25)
  b12 <- sim_boin12_multi(target = 0.35, target_resp = 0.25,
                          scenarios = scen, n_cohort = 12, cohort_size = 3,
                          n_trials = 100, seed = 4)
  boin <- sim_boin_multi(target = 0.35, scenarios = scen, n_cohort = 12,
                         cohort_size = 3, n_trials = 100, seed = 4)
  cmp <- oc_compare(BOIN12 = b12, BOIN = boin, truth = truth)
  r1 <- b12$results[["Scenario 1"]]
  row <- cmp[cmp$design == "BOIN12" & cmp$scenario == "Scenario 1", ]
  expect_equal(row$pct_correct, unname(r1$sel_percent[2]))
  expect_equal(row$n_pts_correct, unname(r1$n_pts_dose[2]))
  expect_equal(row$n_pts_overdose, unname(r1$n_pts_dose[5]))
  expect_identical(row$pct_no_selection, r1$percent_no_obd)
  # Scenario 16 has no true OBD; BOIN is judged against the MTD, dose 2.
  expect_true(is.na(cmp$pct_correct[cmp$design == "BOIN12" &
                                      cmp$scenario == "Scenario 16"]))
  expect_equal(
    cmp$pct_correct[cmp$design == "BOIN" & cmp$scenario == "Scenario 16"],
    unname(boin$results[["Scenario 16"]]$sel_percent[1])
  )
  expect_identical(as.character(cmp$target_type),
                   rep(c("OBD", "MTD"), each = 2))
})

test_that("inconsistent inputs are rejected", {
  truth <- oc_scenarios(list(A = p), target = 0.30)
  expect_error(oc_compare(Fake = list(A = fake_3p3(p))), "oc_scenarios")
  expect_error(oc_compare(Fake = list(B = fake_3p3(p)), truth = truth),
               "is not in 'truth'")
  expect_error(oc_compare(Fake = list(A = fake_3p3(p / 2)), truth = truth),
               "differ between")
  expect_error(oc_compare(Fake = list(A = fake_3p3(p)), truth = truth,
                          overdose_cutoff = 2), "strictly between")
  b12 <- sim_boin12(target = 0.30, target_resp = 0.25, p_true = p,
                    p_resp = c(0.2, 0.3, 0.4, 0.5), n_cohort = 4,
                    cohort_size = 3, n_trials = 20, seed = 1)
  expect_error(oc_compare(b12, truth = oc_scenarios(list(p), target = 0.3)),
               "needs response probabilities")
  m <- matrix(c(0.05, 0.10, 0.20, 0.10, 0.20, 0.30), nrow = 3)
  expect_error(oc_compare(Fake = list(A = fake_3p3(p)),
                          truth = oc_scenarios(list(A = m), target = 0.3)),
               "single agent")
})

test_that("Comb-BOIN12 is judged against the true OBD combinations", {
  p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), 2)
  p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), 2)
  psi00 <- 100 * 0.11 / 1.11
  psi11 <- 100 / 1.11
  oc <- sim_comb_boin12(target = 0.33, target_resp = 0.2, p_true = p_true,
                        p_resp = p_resp, n_cohort = 14, cohort_size = 3,
                        n_trials = 50, psi00 = psi00, psi11 = psi11,
                        est_method = "observed", seed = 1)
  truth <- oc_scenarios(list("Scenario 1" = list(p_true = p_true,
                                                 p_resp = p_resp)),
                        target = 0.33, psi00 = psi00, psi11 = psi11)
  cmp <- oc_compare(oc, truth = truth)
  expect_equal(cmp$pct_correct, oc$pct_true_obd_selected)
  expect_equal(cmp$pct_patients_correct, oc$pct_patients_true_obd)
  expect_equal(cmp$pct_overdose_selected, oc$overdose$pct_selected)
  expect_equal(cmp$pct_no_selection, oc$percent_no_obd)
})
