sc1 <- list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
            p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))

test_that("sim_boin12 summarizes the trials and their OBDs", {
  oc <- sim_boin12(target = 0.35, target_resp = 0.25, p_true = sc1$p_true,
                   p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
                   n_trials = 400, keep_trials = TRUE, seed = 1)
  expect_s3_class(oc, "boin12_oc")
  expect_lt(abs(sum(oc$sel_percent) + oc$percent_no_obd - 100), 1e-10)
  expect_lt(abs(sum(oc$stop_reason_percent) - 100), 1e-10)
  expect_lt(abs(oc$total_n_pts - mean(rowSums(oc$trials$n_pts))), 1e-10)
  trials <- boin12_simulate(target = 0.35, target_resp = 0.25,
                            p_true = sc1$p_true, p_resp = sc1$p_resp,
                            n_cohort = 12, cohort_size = 3, n_trials = 400,
                            seed = 1)
  sel <- boin12_select_obd(trials$n_pts, trials$n_tox, trials$n_resp,
                           trials$n_tox_resp, target = 0.35,
                           target_resp = 0.25, early_stop = trials$early_stop)
  expect_identical(oc$obd, sel$obd)
  expect_identical(oc$trials$n_pts, trials$n_pts)
  expect_equal(oc$est_method, "posterior_mean")
})

test_that("a perfect lowest dose is always the OBD", {
  oc <- sim_boin12(target = 0.35, target_resp = 0.25, p_true = rep(0, 5),
                   p_resp = rep(1, 5), n_cohort = 12, cohort_size = 3,
                   n_trials = 20, seed = 2)
  expect_equal(unname(oc$sel_percent), c(100, 0, 0, 0, 0))
  expect_equal(unname(oc$n_pts_dose), c(33, 3, 0, 0, 0))
})

test_that("the estimate of the utility can be chosen", {
  args <- list(target = 0.35, target_resp = 0.25, p_true = sc1$p_true,
               p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
               n_trials = 200, seed = 3)
  a <- do.call(sim_boin12, c(args, list(est_method = "observed")))
  b <- do.call(sim_boin12, c(args, list(est_method = "tail_probability")))
  expect_equal(a$est_method, "observed")
  expect_equal(b$settings$est_method, "tail_probability")
  expect_equal(a$n_pts_dose, b$n_pts_dose)
  msgs <- capture_messages(
    do.call(sim_boin12, c(args, list(verbose = TRUE)))
  )
  expect_true(any(grepl("Simulating", msgs)))
  expect_true(any(grepl("Selecting the OBD", msgs)))
  expect_error(do.call(sim_boin12, c(args, list(keep_trials = NA))),
               "keep_trials")
})

test_that("the rules for the next dose reach the trials", {
  args <- list(target = 0.35, target_resp = 0.25,
               p_true = c(0.10, 0.20, 0.40, 0.55, 0.70),
               p_resp = c(0.05, 0.10, 0.30, 0.35, 0.40), n_cohort = 12,
               cohort_size = 3, n_trials = 200, tie_break = "random",
               no_candidate = "stop", seed = 3)
  oc <- do.call(sim_boin12, c(args, list(keep_trials = TRUE)))
  trials <- do.call(boin12_simulate, args)
  expect_identical(oc$trials$n_pts, trials$n_pts)
  expect_equal(oc$settings$tie_break, "random")
  expect_equal(oc$settings$no_candidate, "stop")
  expect_named(oc$stop_reason_percent, boin12_stop_reasons())
  # 17 and 45 of the 200 trials stop with no admissible dose and with no
  # admissible candidate (Python with R's random numbers,
  # make_rules_fixture.py).
  expect_equal(unname(oc$stop_reason_percent["no_admissible_dose"]), 8.5)
  expect_equal(unname(oc$stop_reason_percent["no_admissible_candidate"]),
               22.5)
  expect_equal(oc$percent_early_stop, 31)
  expect_gte(oc$percent_no_obd, oc$percent_early_stop)
  expect_equal(oc$settings$explore_from, "any")
  adm <- do.call(sim_boin12, c(args, list(explore_from = "admissible",
                                          keep_trials = TRUE)))
  trials_adm <- do.call(boin12_simulate,
                        c(args, list(explore_from = "admissible")))
  expect_identical(adm$trials$n_pts, trials_adm$n_pts)
  expect_equal(adm$settings$explore_from, "admissible")
  expect_error(do.call(sim_boin12, modifyList(args, list(tie_break = "x"))),
               "should be one of")
  expect_error(do.call(sim_boin12, modifyList(args, list(explore_from = "x"))),
               "should be one of")
})
