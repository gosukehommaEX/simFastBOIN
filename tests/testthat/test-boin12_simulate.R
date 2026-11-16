sc1 <- list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
            p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))

test_that("boin12_simulate returns a well formed object", {
  trials <- boin12_simulate(target = 0.35, target_resp = 0.25,
                            p_true = sc1$p_true, p_resp = sc1$p_resp,
                            n_cohort = 12, cohort_size = 3, n_trials = 300,
                            seed = 1)
  expect_s3_class(trials, "boin12_trials")
  for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
    expect_equal(dim(trials[[nm]]), c(300L, 5L))
    expect_type(trials[[nm]], "integer")
  }
  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(trials$n_resp <= trials$n_pts))
  expect_true(all(trials$n_tox_resp <= pmin(trials$n_tox, trials$n_resp)))
  expect_true(all(rowSums(trials$n_pts) <= 36L))
  expect_true(all(trials$stop_reason %in% boin12_stop_reasons()))
  expect_identical(trials$early_stop,
                   trials$stop_reason == "no_admissible_dose")
  expect_equal(trials$settings$max_total_pts, 36)
  expect_equal(trials$design, "BOIN12")
})

test_that("a perfect lowest dose keeps the trial there", {
  # No DLT and a response in every patient: the lowest dose stays the most
  # desirable, and the exploration rule tries dose 2 once, after 9 patients.
  args <- list(target = 0.35, target_resp = 0.25, p_true = rep(0, 5),
               p_resp = rep(1, 5), n_cohort = 12, cohort_size = 3,
               n_trials = 3, seed = 2)
  trials <- do.call(boin12_simulate, args)
  for (i in 1:3) {
    expect_equal(unname(trials$n_pts[i, ]), c(33L, 3L, 0L, 0L, 0L))
  }
  expect_true(all(trials$n_resp == trials$n_pts))
  expect_true(all(trials$stop_reason == "max_cohorts"))
  early <- do.call(boin12_simulate, c(args, list(n_earlystop = 12)))
  expect_equal(unname(early$n_pts[1, ]), c(12L, 3L, 0L, 0L, 0L))
  expect_true(all(early$stop_reason == "n_earlystop"))
  simple <- do.call(boin12_simulate,
                    c(args, list(n_earlystop = 3,
                                 n_earlystop_rule = "simple")))
  expect_equal(unname(simple$n_pts[1, ]), c(3L, 0L, 0L, 0L, 0L))
})

test_that("a toxic lowest dose stops the trial", {
  trials <- boin12_simulate(target = 0.35, target_resp = 0.25,
                            p_true = rep(1, 4), p_resp = rep(0.5, 4),
                            n_cohort = 10, cohort_size = 3, n_trials = 4,
                            seed = 3)
  expect_true(all(trials$stop_reason == "no_admissible_dose"))
  expect_true(all(trials$early_stop))
  expect_true(all(rowSums(trials$n_pts) == 3L))
})

test_that("independent outcomes are the same under both links", {
  args <- list(target = 0.35, target_resp = 0.25, p_true = sc1$p_true,
               p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
               n_trials = 200, resp_cor = 0, seed = 4)
  a <- do.call(boin12_simulate, args)
  b <- do.call(boin12_simulate, c(args, list(resp_cor_type = "binary")))
  expect_identical(a$n_pts, b$n_pts)
  expect_identical(a$n_resp, b$n_resp)
  expect_identical(a$n_tox_resp, b$n_tox_resp)
})

test_that("the seed reproduces the trials and the session state is kept", {
  args <- list(target = 0.35, target_resp = 0.25, p_true = sc1$p_true,
               p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
               n_trials = 100, resp_cor = 0.3, seed = 5)
  set.seed(99)
  before <- .Random.seed
  a <- do.call(boin12_simulate, args)
  expect_identical(.Random.seed, before)
  b <- do.call(boin12_simulate, args)
  expect_identical(a$n_pts, b$n_pts)
  expect_identical(a$n_tox_resp, b$n_tox_resp)
})

test_that("invalid arguments are rejected", {
  args <- list(target = 0.35, target_resp = 0.25, p_true = sc1$p_true,
               p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3,
               n_trials = 10)
  expect_error(do.call(boin12_simulate,
                       modifyList(args, list(p_resp = c(0.2, 0.3)))),
               "p_resp")
  expect_error(do.call(boin12_simulate,
                       modifyList(args, list(start_dose = 6))),
               "start_dose")
  expect_error(do.call(boin12_simulate,
                       modifyList(args, list(n_star = 0))), "n_star")
  expect_error(do.call(boin12_simulate,
                       modifyList(args, list(resp_cor = 1))), "resp_cor")
  expect_error(do.call(boin12_simulate,
                       modifyList(args, list(n_earlystop = 0))),
               "n_earlystop")
  expect_error(do.call(boin12_simulate,
                       modifyList(args, list(resp_cor = 0.9,
                                             resp_cor_type = "binary"))),
               "attainable")
})
