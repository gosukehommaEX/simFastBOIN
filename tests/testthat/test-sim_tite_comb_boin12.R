sc16 <- list(p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), 2),
             p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), 2))
tite_args <- list(target = 0.33, target_resp = 0.2, p_true = sc16$p_true,
                  p_resp = sc16$p_resp, n_cohort = 14, cohort_size = 3,
                  window = 1, resp_window = 2, accrual_rate = 2,
                  event_time = "uniform", psi00 = 100 * 0.11 / 1.11,
                  psi11 = 100 / 1.11, benchmark_weight = 0.6,
                  extrasafe = TRUE, offset = 0.10)

test_that("sim_tite_comb_boin12 summarizes its trials and their durations", {
  oc <- do.call(sim_tite_comb_boin12,
                c(tite_args, list(n_trials = 200, est_method = "observed",
                                  obd_bound = "contour", keep_trials = TRUE,
                                  seed = 5)))
  expect_s3_class(oc, "comb_boin12_oc")
  expect_equal(oc$design, "TITE-Comb-BOIN12")
  trials <- do.call(tite_comb_boin12_simulate,
                    c(tite_args, list(n_trials = 200, seed = 5)))
  expect_identical(oc$trials$n_pts, trials$n_pts)
  expect_equal(oc$duration, mean(trials$duration))
  sel <- comb_boin12_select_obd(
    trials$n_pts, trials$n_tox, trials$n_resp, target = 0.33,
    target_resp = 0.2, extrasafe = TRUE, offset = 0.10,
    psi00 = tite_args$psi00, psi11 = tite_args$psi11,
    est_method = "observed", obd_bound = "contour",
    early_stop = trials$early_stop)
  expect_identical(unname(oc$obd[, "dose_a"]), sel$obd_a)
  expect_identical(unname(oc$obd[, "dose_b"]), sel$obd_b)
  expect_equal(sum(oc$sel_percent) + oc$percent_no_obd, 100)
  expect_equal(oc$obd_bound, "contour")
  out <- utils::capture.output(res <- print(oc))
  expect_identical(res, oc)
  expect_true(any(grepl("TITE-Comb-BOIN12 operating characteristics", out,
                        fixed = TRUE)))
  expect_true(any(grepl("Average trial duration", out, fixed = TRUE)))
})

test_that("verbose prints the steps as messages", {
  msgs <- testthat::capture_messages(
    do.call(sim_tite_comb_boin12,
            c(tite_args, list(n_trials = 20, est_method = "observed",
                              verbose = TRUE)))
  )
  expect_true(any(grepl("Simulating", msgs)))
  expect_true(any(grepl("Selecting the OBD combination", msgs)))
})
