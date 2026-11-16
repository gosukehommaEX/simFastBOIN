sc1 <- list(p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
            p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70))

test_that("sim_boin_et summarizes the selection and the allocation", {
  oc <- sim_boin_et(target = 0.30, target_resp = 0.60, p_true = sc1$p_true,
                    p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
                    n_trials = 1000, seed = 1)

  expect_s3_class(oc, "boin_et_oc")
  expect_identical(names(oc$sel_percent), paste0("DL", 1:6))
  expect_equal(sum(oc$sel_percent) + oc$percent_no_obd, 100)
  expect_equal(sum(oc$stop_reason_percent), 100)
  expect_equal(oc$total_n_pts, sum(oc$n_pts_dose))
  expect_lte(oc$total_n_pts, 36)
  expect_lte(oc$percent_early_stop, oc$percent_no_obd)
  expect_equal(c(oc$lambda1, oc$lambda2, oc$eta1), c(0.14, 0.35, 0.48))
  expect_identical(oc$est_method, "observed")
  expect_identical(oc$obd_method, "max_resp")
  expect_identical(oc$design, "BOIN-ET")
  expect_null(oc$duration)
})

test_that("the summary is computed from the trials and their OBDs", {
  args <- list(target = 0.30, target_resp = 0.60, p_true = sc1$p_true,
               p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
               n_trials = 400, seed = 2)
  oc <- do.call(sim_boin_et, c(args, keep_trials = TRUE,
                               obd_method = "utility_scoring"))
  tr <- do.call(boin_et_simulate, args)
  sel <- boin_et_select_obd(tr$n_pts, tr$n_tox, tr$n_resp, target = 0.30,
                            target_resp = 0.60, obd_method = "utility_scoring",
                            early_stop = tr$early_stop)

  expect_identical(oc$trials$n_pts, tr$n_pts)
  expect_identical(oc$obd, sel$obd)
  expected <- as.numeric(table(factor(sel$obd, levels = 1:6))) / 400 * 100
  expect_equal(unname(oc$sel_percent), expected)
  expect_equal(unname(oc$n_resp_dose), unname(colMeans(tr$n_resp)))
  expect_equal(oc$percent_no_obd, mean(is.na(sel$obd)) * 100)
})

test_that("degenerate scenarios give the expected selections", {
  # No DLT and no response: every trial ends with 3, 3 and 12 patients (3, 3
  # and 30 under the rules of boinet) and selects the lowest dose under every
  # rule.
  for (m in c("max_resp", "utility_weighted", "utility_truncated_linear",
              "utility_scoring")) {
    oc <- sim_boin_et(target = 0.30, target_resp = 0.60, p_true = c(0, 0, 0),
                      p_resp = c(0, 0, 0), n_cohort = 12, cohort_size = 3,
                      n_trials = 50, obd_method = m, seed = 4)
    expect_equal(unname(oc$sel_percent), c(100, 0, 0))
    expect_equal(unname(oc$n_pts_dose), c(3, 3, 12))
    expect_equal(unname(oc$stop_reason_percent["top_dose_excluded"]), 100)
    oc <- sim_boin_et(target = 0.30, target_resp = 0.60, p_true = c(0, 0, 0),
                      p_resp = c(0, 0, 0), n_cohort = 12, cohort_size = 3,
                      n_trials = 50, obd_method = m, boinet_compatible = TRUE,
                      seed = 4)
    expect_equal(unname(oc$sel_percent), c(100, 0, 0))
    expect_equal(unname(oc$n_pts_dose), c(3, 3, 30))
  }

  oc <- sim_boin_et(target = 0.30, target_resp = 0.60, p_true = c(1, 1),
                    p_resp = c(0.5, 0.5), n_cohort = 6, cohort_size = 3,
                    n_trials = 50, seed = 4)
  expect_equal(oc$percent_no_obd, 100)
  expect_equal(oc$percent_early_stop, 100)
})

test_that("the multi_iso estimator can be chosen", {
  oc <- sim_boin_et(target = 0.30, target_resp = 0.60, p_true = sc1$p_true,
                    p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
                    n_trials = 300, est_method = "multi_iso", seed = 3)
  expect_identical(oc$est_method, "multi_iso")
  expect_equal(sum(oc$sel_percent) + oc$percent_no_obd, 100)
})

test_that("verbose reports progress", {
  msgs <- capture_messages(
    sim_boin_et(target = 0.30, target_resp = 0.60, p_true = sc1$p_true,
                p_resp = sc1$p_resp, n_cohort = 4, cohort_size = 3,
                n_trials = 20, verbose = TRUE, seed = 1)
  )
  expect_true(any(grepl("Simulating", msgs)))
  expect_true(any(grepl("Selecting the OBD", msgs)))
})
