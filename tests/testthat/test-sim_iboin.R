zhou_oc <- function(...) {
  sim_iboin(
    target = 0.30, p_true = c(0.08, 0.15, 0.31, 0.45, 0.55),
    skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54), prior_ess = 3,
    n_cohort = 10, cohort_size = 3, ...
  )
}

test_that("sim_iboin gives the operating characteristics of the rules", {
  # Expected values from the C++ core run outside of R on R's uniform stream
  # for this seed, and from an implementation of the rules in Python
  oc <- zhou_oc(n_trials = 1000, n_earlystop = 100, seed = 2021)

  expect_s3_class(oc, "iboin_oc")
  expect_s3_class(oc, "boin_oc")
  expect_equal(unname(oc$sel_percent), c(0.6, 23.4, 60.8, 15.0, 0.2),
               tolerance = 1e-10)
  expect_equal(oc$percent_no_mtd, 0, tolerance = 1e-10)
  expect_true(oc$select_with_prior)
  expect_equal(unname(oc$n_pts_dose), c(3.651, 9.267, 12.675, 3.834, 0.573),
               tolerance = 1e-10)
  expect_equal(unname(oc$n_tox_dose), c(0.271, 1.374, 3.910, 1.714, 0.313),
               tolerance = 1e-10)
  expect_equal(oc$total_n_pts, 30, tolerance = 1e-10)
  expect_equal(oc$total_n_tox, 7.582, tolerance = 1e-10)

  # The MTD selected from the data alone: same trials, other selection
  data_only <- zhou_oc(n_trials = 1000, n_earlystop = 100, seed = 2021,
                       select_with_prior = FALSE)
  expect_equal(unname(data_only$sel_percent), c(0.6, 26.0, 54.8, 17.2, 1.4),
               tolerance = 1e-10)
  expect_identical(data_only$n_pts_dose, oc$n_pts_dose)
  expect_false(data_only$select_with_prior)

  early <- zhou_oc(n_trials = 1000, seed = 2021)
  expect_equal(unname(early$sel_percent), c(0.8, 24.4, 59.9, 14.1, 0.7),
               tolerance = 1e-10)
  expect_equal(early$percent_no_mtd, 0.1, tolerance = 1e-10)
  expect_equal(unname(early$n_pts_dose),
               c(3.633, 8.670, 12.375, 3.927, 0.564), tolerance = 1e-10)
  expect_equal(early$total_n_pts, 29.169, tolerance = 1e-10)
  reasons <- c("lowest_dose_eliminated", "n_earlystop", "max_sample_size")
  expect_equal(as.numeric(early$stop_reason_percent[reasons]),
               c(0.1, 19.8, 80.1), tolerance = 1e-10)
})

test_that("the robust prior helps when the prior MTD is too low", {
  # Prior MTD 3, true MTD 5. Expected values as above.
  args <- list(
    target = 0.30, p_true = c(0.05, 0.08, 0.10, 0.14, 0.30),
    skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54), prior_ess = 3,
    n_cohort = 10, cohort_size = 3, n_trials = 500, n_earlystop = 100,
    seed = 11
  )
  robust <- do.call(sim_iboin, c(args, list(robust = TRUE)))
  plain <- do.call(sim_iboin, args)

  expect_equal(unname(robust$sel_percent), c(0.0, 0.2, 3.6, 22.8, 73.4),
               tolerance = 1e-10)
  expect_equal(unname(robust$n_pts_dose),
               c(3.102, 4.038, 4.776, 6.936, 11.148), tolerance = 1e-10)
  expect_equal(unname(plain$sel_percent), c(0.0, 0.2, 3.8, 48.6, 47.4),
               tolerance = 1e-10)
  expect_equal(unname(plain$n_pts_dose),
               c(3.102, 4.044, 6.306, 8.526, 8.022), tolerance = 1e-10)
  expect_identical(robust$prior_ess, c(3L, 3L, 3L, 0L, 0L))
  expect_identical(plain$prior_ess, rep(3L, 5))
})

test_that("without prior information sim_iboin agrees with sim_boin", {
  p <- c(0.08, 0.15, 0.31, 0.45, 0.55)
  new <- sim_iboin(0.30, p, c(0.1, 0.2, 0.3, 0.4, 0.5), 0, 10, 3,
                   n_trials = 400, seed = 31, keep_trials = TRUE)
  ref <- sim_boin(0.30, p, 10, 3, n_trials = 400, seed = 31,
                  keep_trials = TRUE)

  expect_identical(new$sel_percent, ref$sel_percent)
  expect_identical(new$n_pts_dose, ref$n_pts_dose)
  expect_identical(new$n_tox_dose, ref$n_tox_dose)
  expect_identical(new$percent_no_mtd, ref$percent_no_mtd)
  expect_identical(new$overdose, ref$overdose)
  expect_identical(new$trials$mtd, ref$trials$mtd)
  expect_identical(new$lambda_e, ref$lambda_e)
  expect_identical(new$lambda_d, ref$lambda_d)
})

test_that("the trial level data hold the selected doses", {
  # Expected doses as above
  prior <- zhou_oc(n_trials = 6, n_earlystop = 100, seed = 7,
                   keep_trials = TRUE)
  data_only <- zhou_oc(n_trials = 6, n_earlystop = 100, seed = 7,
                       keep_trials = TRUE, select_with_prior = FALSE)

  expect_identical(prior$trials$mtd, c(2L, 3L, 3L, 3L, 3L, 4L))
  expect_identical(data_only$trials$mtd, c(2L, 4L, 3L, 3L, 4L, 4L))
  expect_identical(prior$trials$n_pts, data_only$trials$n_pts)
  expect_error(zhou_oc(n_trials = 6, select_with_prior = NA),
               "select_with_prior")
})

test_that("sim_iboin returns a coherent summary", {
  oc <- zhou_oc(n_trials = 300, seed = 4, keep_trials = TRUE)

  expect_equal(sum(oc$sel_percent) + oc$percent_no_mtd, 100,
               tolerance = 1e-10)
  expect_equal(sum(oc$n_pts_dose), oc$total_n_pts, tolerance = 1e-10)
  expect_equal(unname(oc$n_pts_dose), unname(colMeans(oc$trials$n_pts)))
  expect_identical(oc$design, "iBOIN")
  expect_identical(oc$skeleton, c(0.10, 0.19, 0.30, 0.42, 0.54))
  expect_s3_class(oc$boundary, "iboin_boundary")
  expect_identical(oc$overdose$doses, 3:5)
  expect_null(zhou_oc(n_trials = 50)$trials)
  reported <- capture_messages(zhou_oc(n_trials = 50, verbose = TRUE))
  expect_true(any(grepl("Simulating 50 trials", reported)))
})

test_that("the result prints and enters the comparison tables as iBOIN", {
  oc <- zhou_oc(n_trials = 200, seed = 4)

  expect_output(print(oc), "iBOIN operating characteristics")
  expect_output(print(oc), "skeleton        : 0.10 0.19 0.30 0.42 0.54")
  expect_output(print(oc), "prior ESS       : 3 3 3 3 3")
  expect_output(print(oc), "MTD selected")
  out <- capture.output(print(oc))
  expect_false(any(grepl("lambda_e", out)))

  boin <- sim_boin(0.30, c(0.08, 0.15, 0.31, 0.45, 0.55), 10, 3,
                   n_trials = 200, seed = 4)
  res <- oc_results(oc, boin)
  expect_identical(levels(res$design), c("iBOIN", "BOIN"))
})
