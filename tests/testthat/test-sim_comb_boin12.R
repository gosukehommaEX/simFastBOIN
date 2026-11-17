sc16 <- list(p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), 2),
             p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), 2))
lu_args <- list(target = 0.33, target_resp = 0.2, p_true = sc16$p_true,
                p_resp = sc16$p_resp, n_cohort = 14, cohort_size = 3,
                psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11,
                benchmark_weight = 0.6, extrasafe = TRUE, offset = 0.10,
                dose_a = c(100, 200), dose_b = c(50, 100, 200))

test_that("sim_comb_boin12 returns consistent operating characteristics", {
  oc <- do.call(sim_comb_boin12, c(lu_args, list(n_trials = 200, seed = 3,
                                                 keep_trials = TRUE)))
  expect_s3_class(oc, "comb_boin12_oc")
  expect_equal(oc$design, "Comb-BOIN12")
  expect_equal(dim(oc$sel_percent), c(2L, 3L))
  expect_equal(sum(oc$sel_percent) + oc$percent_no_obd, 100)
  expect_equal(sum(oc$stop_reason_percent), 100)
  expect_equal(oc$total_n_pts, sum(oc$n_pts_dose))
  expect_equal(oc$n_pts_dose,
               apply(oc$trials$n_pts, c(2L, 3L), mean),
               ignore_attr = TRUE)
  # Scenario 16: the true OBDC is (1, 3) and only (2, 3) is above 0.33.
  expect_equal(unname(oc$true_obd), matrix(c(1L, 3L), 1))
  expect_equal(oc$pct_true_obd_selected, oc$sel_percent[1, 3])
  expect_equal(oc$overdose$pct_selected, oc$sel_percent[2, 3])
  expect_equal(oc$pct_patients_true_obd,
               oc$n_pts_dose[1, 3] / oc$total_n_pts * 100)
  # The kept selections are those of comb_boin12_select_obd().
  sel <- comb_boin12_select_obd(
    oc$trials$n_pts, oc$trials$n_tox, oc$trials$n_resp, target = 0.33,
    target_resp = 0.2, extrasafe = TRUE, offset = 0.10,
    psi00 = lu_args$psi00, psi11 = lu_args$psi11, dose_a = c(100, 200),
    dose_b = c(50, 100, 200), early_stop = oc$trials$early_stop)
  expect_identical(unname(oc$obd[, "dose_a"]), sel$obd_a)
  expect_identical(unname(oc$obd[, "dose_b"]), sel$obd_b)
  expect_equal(oc$est_method, "logistic")
  expect_equal(oc$obd_bound, "rectangle")
  expect_equal(oc$settings$dose_b, c(50, 100, 200))
})

test_that("the selection settings change only the selection", {
  args <- c(lu_args, list(n_trials = 200, seed = 4))
  obs <- do.call(sim_comb_boin12, c(args, list(est_method = "observed")))
  cont <- do.call(sim_comb_boin12, c(args, list(est_method = "observed",
                                                obd_bound = "contour")))
  expect_identical(obs$n_pts_dose, cont$n_pts_dose)
  expect_identical(obs$percent_early_stop, cont$percent_early_stop)
  expect_false(identical(obs$sel_percent, cont$sel_percent))
  expect_equal(cont$obd_bound, "contour")
  expect_null(obs$trials)
})

test_that("verbose prints progress messages", {
  msgs <- capture_messages(
    do.call(sim_comb_boin12, c(lu_args, list(n_trials = 20, seed = 1,
                                             est_method = "observed",
                                             verbose = TRUE)))
  )
  expect_true(any(grepl("Simulating", msgs)))
  expect_true(any(grepl("Selecting the OBD combination", msgs)))
  expect_error(do.call(sim_comb_boin12,
                       c(lu_args, list(n_trials = 20, keep_trials = NA))),
               "keep_trials")
})
