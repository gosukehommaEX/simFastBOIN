sc7 <- list(
  p_true = matrix(c(0.01, 0.04, 0.11, 0.15, 0.30, 0.03, 0.05, 0.13, 0.30,
                    0.50, 0.07, 0.10, 0.30, 0.48, 0.54), 3, byrow = TRUE),
  p_resp = matrix(c(0.05, 0.10, 0.20, 0.30, 0.60, 0.10, 0.20, 0.30, 0.60,
                    0.65, 0.20, 0.30, 0.60, 0.65, 0.70), 3, byrow = TRUE)
)

test_that("sim_boin_etc summarizes the trials and their selection", {
  oc <- sim_boin_etc(target = 0.30, target_resp = 0.60, p_true = sc7$p_true,
                     p_resp = sc7$p_resp, n_cohort = 16, cohort_size = 3,
                     n_trials = 200, subtrials = "etc2", keep_trials = TRUE,
                     seed = 5)
  expect_s3_class(oc, "boin_etc_oc")
  expect_equal(oc$design, "BOIN-ETC2")
  trials <- boin_etc_simulate(target = 0.30, target_resp = 0.60,
                              p_true = sc7$p_true, p_resp = sc7$p_resp,
                              n_cohort = 16, cohort_size = 3, n_trials = 200,
                              subtrials = "etc2", seed = 5)
  expect_identical(oc$trials$n_pts, trials$n_pts)
  sel <- boin_etc_select_obd(trials$n_pts, trials$n_tox, trials$n_resp,
                             trials$n_tox_resp, target = 0.30,
                             excluded = trials$excluded,
                             early_stop = trials$early_stop)
  expect_identical(unname(oc$obd[, "dose_a"]), sel$obd_a)
  cell <- table(factor(sel$obd_a + 3L * (sel$obd_b - 1L), levels = 1:15))
  expect_equal(as.vector(oc$sel_percent), as.vector(cell) / 2)
  expect_equal(sum(oc$sel_percent) + oc$percent_no_obd, 100)
  expect_equal(unname(oc$n_pts_dose),
               apply(array(as.numeric(trials$n_pts), dim(trials$n_pts)),
                     c(2L, 3L), mean))
  expect_equal(oc$total_n_pts, mean(apply(trials$n_pts, 1L, sum)))
  expect_equal(oc$mean_subtrials, mean(trials$n_subtrials))
  expect_equal(sum(oc$stop_reason_percent), 100)
  expect_equal(c(oc$lambda1, oc$lambda2, oc$eta1), c(0.14, 0.35, 0.48))
  # The true OBDs by utility: A3B2 and A2B3 have the highest true utility
  # among the combinations with a DLT probability at most 0.30.
  expect_true(nrow(oc$true_obd) >= 1L)
  expect_equal(oc$pct_true_obd_selected,
               sum(oc$sel_percent[oc$true_obd]))
})

test_that("odc_tie = \"all\" counts every tied candidate", {
  args <- list(target = 0.30, target_resp = 0.60, p_true = sc7$p_true,
               p_resp = sc7$p_resp, n_cohort = 16, cohort_size = 3,
               n_trials = 300, subtrials = "etc1", etc1_start = "below",
               rule4 = "authors_code", seed = 8)
  low <- do.call(sim_boin_etc, args)
  all <- do.call(sim_boin_etc, c(args, list(odc_tie = "all",
                                             keep_trials = TRUE)))
  expect_true(all(all$sel_percent >= low$sel_percent - 1e-12))
  expect_gt(sum(all$sel_percent), sum(low$sel_percent))
  expect_equal(as.vector(all$sel_percent),
               as.vector(apply(array(as.numeric(all$tied), dim(all$tied)),
                               c(2L, 3L), mean)) * 100)
  expect_equal(all$percent_no_obd, low$percent_no_obd)
  expect_equal(all$settings$rule4, "authors_code")
})

test_that("oc_results() and oc_compare() accept sim_boin_etc() results", {
  oc <- sim_boin_etc(target = 0.30, target_resp = 0.60, p_true = sc7$p_true,
                     p_resp = sc7$p_resp, n_cohort = 16, cohort_size = 3,
                     n_trials = 100, subtrials = "etc3", seed = 2)
  res <- oc_results(ETC3 = oc)
  expect_s3_class(res, "oc_results")
  expect_true(all(res$design == "ETC3"))
  expect_true(all(res$target_type == "OBD"))
  sel <- res[res$metric == "sel_percent", , drop = FALSE]
  expect_equal(sum(sel$value), sum(oc$sel_percent))
})

test_that("sim_boin_etc checks its arguments", {
  expect_error(sim_boin_etc(target = 0.30, target_resp = 0.60,
                            p_true = sc7$p_true, p_resp = sc7$p_resp,
                            n_cohort = 16, cohort_size = 3, n_trials = 10,
                            odc_tie = "first"))
  expect_error(sim_boin_etc(target = 0.30, target_resp = 0.60,
                            p_true = sc7$p_true, p_resp = sc7$p_resp,
                            n_cohort = 16, cohort_size = 3, n_trials = 10,
                            keep_trials = NA), "keep_trials")
})
