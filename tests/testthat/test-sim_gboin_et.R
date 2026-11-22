test_that("sim_gboin_et gives the operating characteristics of the rules", {
  # Expected values from the C++ core run outside of R on R's uniform stream
  # for this seed, and from an implementation of the rules in Python
  oc <- takeda_oc(n_trials = 1000, seed = 2022)

  expect_s3_class(oc, "boin_et_oc")
  expect_identical(oc$design, "gBOIN-ET")
  expect_equal(unname(oc$sel_percent), c(2.6, 9.6, 28.1, 56.0, 3.5, 0.0),
               tolerance = 1e-10)
  expect_equal(oc$percent_no_obd, 0.2, tolerance = 1e-10)
  expect_equal(oc$percent_early_stop, 0.2, tolerance = 1e-10)
  expect_equal(unname(oc$n_pts_dose),
               c(3.375, 5.325, 9.207, 10.578, 1.425, 0.078),
               tolerance = 1e-10)
  expect_lt(max(abs(unname(oc$n_tox_dose) -
                      c(0.243000, 0.682333, 1.998333, 3.381000, 0.690000,
                        0.055333))), 1e-6)
  expect_lt(max(abs(unname(oc$n_resp_dose) -
                      c(0.276333, 1.280333, 3.477667, 6.740333, 0.977667,
                        0.057667))), 1e-6)
  expect_equal(oc$total_n_pts, 29.988, tolerance = 1e-10)
  expect_lt(abs(oc$total_n_tox - 7.05), 1e-6)
  expect_lt(abs(oc$total_n_resp - 12.81), 1e-6)
  expect_identical(oc$est_method, "observed")
  expect_identical(oc$obd_method, "max_resp")
})

test_that("the OBD selection of the article gives its own result", {
  # Expected values as above. The trials are those of the first test.
  oc <- takeda_oc(n_trials = 1000, seed = 2022, est_method = "fp_linear",
                  min_obd_sample = 6)

  expect_equal(unname(oc$sel_percent), c(1.2, 3.8, 34.7, 56.4, 3.4, 0.0),
               tolerance = 1e-10)
  expect_equal(oc$percent_no_obd, 0.5, tolerance = 1e-10)
  expect_equal(unname(oc$n_pts_dose),
               c(3.375, 5.325, 9.207, 10.578, 1.425, 0.078),
               tolerance = 1e-10)
  expect_identical(oc$settings$min_obd_sample, 6)
})

test_that("the rules of boinet and a utility give their own result", {
  # Scenario 5 of Takeda et al. (2022). Expected values as above.
  oc <- takeda_oc(5, n_trials = 1000, seed = 5, boinet_compatible = TRUE,
                  resp_cor = 0, obd_method = "utility_scoring")

  expect_equal(unname(oc$sel_percent), c(3.9, 9.6, 26.6, 49.9, 7.8, 2.2),
               tolerance = 1e-10)
  expect_equal(oc$percent_no_obd, 0, tolerance = 1e-10)
  expect_equal(unname(oc$n_pts_dose),
               c(4.110, 6.615, 8.862, 8.313, 1.476, 0.624),
               tolerance = 1e-10)
  expect_equal(as.numeric(oc$stop_reason_percent["n_earlystop"]), 0.2,
               tolerance = 1e-10)
})

test_that("two categories with scores 0 and 1 give the result of sim_boin_et", {
  p <- c(0.05, 0.15, 0.30, 0.45, 0.60)
  q <- c(0.10, 0.30, 0.55, 0.60, 0.62)
  new <- sim_gboin_et(0.30, 0.60, rbind(1 - p, p), rbind(1 - q, q), c(0, 1),
                      c(0, 1), 10, 3, n_trials = 400, seed = 77)
  ref <- sim_boin_et(0.30, 0.60, p, q, 10, 3, n_trials = 400, seed = 77)

  expect_equal(new$sel_percent, ref$sel_percent)
  expect_equal(new$percent_no_obd, ref$percent_no_obd)
  expect_equal(new$n_pts_dose, ref$n_pts_dose)
  expect_equal(new$n_tox_dose, ref$n_tox_dose)
  expect_equal(new$n_resp_dose, ref$n_resp_dose)
  expect_equal(new$stop_reason_percent, ref$stop_reason_percent)
  expect_identical(new$lambda1, ref$lambda1)
  expect_identical(new$eta1, ref$eta1)
})

test_that("the categories are averaged by dose", {
  oc <- takeda_oc(n_trials = 100, seed = 4, keep_trials = TRUE)

  expect_identical(dim(oc$n_grade_dose), c(4L, 6L))
  expect_identical(dim(oc$n_resp_grade_dose), c(3L, 6L))
  expect_equal(colSums(oc$n_grade_dose), oc$n_pts_dose, tolerance = 1e-12)
  expect_equal(colSums(oc$n_resp_grade_dose), oc$n_pts_dose,
               tolerance = 1e-12)
  expect_equal(unname(colSums(oc$n_grade_dose * c(0, 1 / 3, 2 / 3, 1))),
               unname(oc$n_tox_dose), tolerance = 1e-12)
  expect_identical(length(oc$obd), 100L)
  expect_s3_class(oc$trials, "boin_et_trials")
  expect_equal(mean(is.na(oc$obd)) * 100, oc$percent_no_obd)
})

test_that("the result prints with the labels of scores", {
  oc <- takeda_oc(n_trials = 50, seed = 4)
  out <- capture.output(print(oc))

  expect_true(any(grepl("gBOIN-ET operating characteristics", out,
                        fixed = TRUE)))
  expect_true(any(grepl("target toxicity score", out, fixed = TRUE)))
  expect_true(any(grepl("True efficacy score (%)", out, fixed = TRUE)))
  expect_true(any(grepl("Toxicity score", out, fixed = TRUE)))
  expect_true(any(grepl("observed efficacy scores", out, fixed = TRUE)))
  expect_false(any(grepl("DLT", out, fixed = TRUE)))

  out <- capture.output(print(takeda_trials(n_trials = 20)))
  expect_true(any(grepl("Simulated gBOIN-ET trials", out, fixed = TRUE)))
  expect_true(any(grepl("Efficacy score", out, fixed = TRUE)))
  expect_false(any(grepl("DLT", out, fixed = TRUE)))
})

test_that("the result enters the tables of several designs", {
  oc <- takeda_oc(n_trials = 100, seed = 4)
  res <- oc_results(oc)

  expect_true("gBOIN-ET" %in% as.character(res$design))
})

test_that("sim_gboin_et rejects invalid input", {
  expect_error(takeda_oc(est_method = "multi_iso"), "should be one of")
  expect_error(takeda_oc(keep_trials = NA), "keep_trials")
  expect_error(takeda_oc(n_trials = 10, min_obd_sample = 0),
               "min_obd_sample")
  expect_message(
    expect_message(takeda_oc(n_trials = 10, verbose = TRUE), "Simulating"),
    "Selecting"
  )
})
