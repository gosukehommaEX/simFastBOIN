test_that("sim_tite_gboin_et gives the operating characteristics", {
  # Expected values from the C++ core run outside of R on R's uniform stream
  # for this seed, and from an implementation of the rules in Python
  oc <- tite_takeda_oc(n_trials = 1000, seed = 2023)

  expect_s3_class(oc, "boin_et_oc")
  expect_identical(oc$design, "TITE-gBOIN-ET")
  expect_equal(unname(oc$sel_percent), c(0.7, 1.5, 2.5, 6.0, 16.4, 72.9),
               tolerance = 1e-10)
  expect_equal(oc$percent_no_obd, 0, tolerance = 1e-10)
  expect_equal(oc$percent_early_stop, 0, tolerance = 1e-10)
  expect_equal(unname(oc$n_pts_dose),
               c(3.045, 3.207, 3.663, 4.554, 7.704, 13.827),
               tolerance = 1e-10)
  expect_lt(max(abs(unname(oc$n_tox_dose) -
                      c(0.066667, 0.164333, 0.333333, 0.682000, 1.515000,
                        3.946333))), 1e-6)
  expect_lt(max(abs(unname(oc$n_resp_dose) -
                      c(0.164167, 0.248000, 0.422833, 0.850833, 2.297750,
                        7.969417))), 1e-6)
  expect_equal(oc$total_n_pts, 36, tolerance = 1e-10)
  expect_lt(abs(oc$total_n_tox - 6.707667), 1e-6)
  expect_lt(abs(oc$total_n_resp - 11.953), 1e-6)
  expect_lt(abs(oc$duration - 560.232071), 1e-5)
  expect_identical(oc$est_method, "observed")
  expect_identical(oc$obd_method, "max_resp")
  expect_identical(oc$settings$resp_window, 60)
})

test_that("the settings of the article and a utility give their own result", {
  # Scenario 3 of Takeda et al. (2023) with the rounded targets and the
  # printed boundaries of the article, its efficacy cutoff, independent
  # outcomes and the rules of 'boinet'. Expected values as above.
  args <- list(3, target = 0.313, target_resp = 0.583, lambda1 = 0.14,
               lambda2 = 0.37, eta1 = 0.46, cutoff_eli_resp = 0.95,
               resp_cor = 0, boinet_compatible = TRUE, n_trials = 1000,
               seed = 5)
  oc <- do.call(tite_takeda_oc, c(args, obd_method = "utility_scoring"))

  expect_equal(unname(oc$sel_percent), c(5.5, 9.0, 9.5, 58.3, 14.7, 2.3),
               tolerance = 1e-10)
  expect_equal(oc$percent_no_obd, 0.7, tolerance = 1e-10)
  expect_equal(unname(oc$n_pts_dose),
               c(3.411, 5.262, 9.039, 14.151, 3.408, 0.597),
               tolerance = 1e-10)
  expect_equal(as.numeric(oc$stop_reason_percent["no_lower_admissible_dose"]),
               0.7, tolerance = 1e-10)
  expect_lt(abs(oc$duration - 519.197011), 1e-5)

  # The selection of the article on the same trials
  oc <- do.call(tite_takeda_oc, args)
  expect_equal(unname(oc$sel_percent), c(2.7, 8.4, 21.0, 61.0, 6.1, 0.1),
               tolerance = 1e-10)
  expect_equal(oc$percent_no_obd, 0.7, tolerance = 1e-10)
})

test_that("two categories with scores 0 and 1 give sim_tite_boin_et()", {
  p <- c(0.05, 0.15, 0.30, 0.45, 0.60)
  q <- c(0.10, 0.30, 0.55, 0.60, 0.62)
  new <- sim_tite_gboin_et(0.30, 0.60, rbind(1 - p, p), rbind(1 - q, q),
                           c(0, 1), c(0, 1), 10, 3, window = 30,
                           accrual_rate = 0.1, resp_window = 45,
                           n_trials = 400, seed = 77)
  ref <- sim_tite_boin_et(0.30, 0.60, p, q, 10, 3, window = 30,
                          accrual_rate = 0.1, resp_window = 45,
                          n_trials = 400, seed = 77)

  expect_equal(new$sel_percent, ref$sel_percent)
  expect_equal(new$percent_no_obd, ref$percent_no_obd)
  expect_equal(new$n_pts_dose, ref$n_pts_dose)
  expect_equal(new$n_tox_dose, ref$n_tox_dose)
  expect_equal(new$n_resp_dose, ref$n_resp_dose)
  expect_equal(new$stop_reason_percent, ref$stop_reason_percent)
  expect_equal(new$duration, ref$duration)
  expect_identical(new$lambda1, ref$lambda1)
  expect_identical(new$eta1, ref$eta1)
})

test_that("the categories are averaged by dose and the trials are kept", {
  oc <- tite_takeda_oc(n_trials = 100, seed = 4, keep_trials = TRUE)

  expect_identical(dim(oc$n_grade_dose), c(4L, 6L))
  expect_identical(dim(oc$n_resp_grade_dose), c(4L, 6L))
  expect_equal(colSums(oc$n_grade_dose), oc$n_pts_dose, tolerance = 1e-12)
  expect_equal(colSums(oc$n_resp_grade_dose), oc$n_pts_dose,
               tolerance = 1e-12)
  expect_equal(unname(colSums(oc$n_grade_dose * c(0, 1 / 3, 2 / 3, 1))),
               unname(oc$n_tox_dose), tolerance = 1e-12)
  expect_identical(length(oc$obd), 100L)
  expect_s3_class(oc$trials, "boin_et_trials")
  expect_equal(mean(is.na(oc$obd)) * 100, oc$percent_no_obd)
  expect_equal(oc$duration, mean(oc$trials$duration))

  sel <- gboin_et_select_obd(oc$trials$n_pts, oc$trials$n_tox,
                             oc$trials$n_resp, target = 0.47 / 1.5,
                             target_resp = 1.75 / 3,
                             early_stop = oc$trials$early_stop)
  expect_identical(oc$obd, sel$obd)
})

test_that("the result prints with the labels of scores and the duration", {
  oc <- tite_takeda_oc(n_trials = 50, seed = 4)
  out <- capture.output(print(oc))

  expect_true(any(grepl("TITE-gBOIN-ET operating characteristics", out,
                        fixed = TRUE)))
  expect_true(any(grepl("target toxicity score", out, fixed = TRUE)))
  expect_true(any(grepl("True efficacy score (%)", out, fixed = TRUE)))
  expect_true(any(grepl("Toxicity score", out, fixed = TRUE)))
  expect_true(any(grepl("observed efficacy scores", out, fixed = TRUE)))
  expect_true(any(grepl("Average trial duration", out, fixed = TRUE)))
  expect_false(any(grepl("DLT", out, fixed = TRUE)))

  out <- capture.output(print(tite_takeda_trials(n_trials = 20)))
  expect_true(any(grepl("Simulated TITE-gBOIN-ET trials", out,
                        fixed = TRUE)))
  expect_true(any(grepl("Efficacy score", out, fixed = TRUE)))
  expect_true(any(grepl("Average trial duration", out, fixed = TRUE)))
  expect_false(any(grepl("DLT", out, fixed = TRUE)))
})

test_that("the result enters the tables of several designs", {
  oc <- tite_takeda_oc(n_trials = 100, seed = 4)
  res <- oc_results(oc)

  expect_true("TITE-gBOIN-ET" %in% as.character(res$design))
})

test_that("sim_tite_gboin_et rejects invalid input", {
  expect_error(tite_takeda_oc(est_method = "multi_iso"), "should be one of")
  expect_error(tite_takeda_oc(keep_trials = NA), "keep_trials")
  expect_error(tite_takeda_oc(n_trials = 10, min_obd_sample = 0),
               "min_obd_sample")
  expect_message(
    expect_message(tite_takeda_oc(n_trials = 10, verbose = TRUE),
                   "Simulating"),
    "Selecting"
  )
})
