uboin_check_trials <- function(trials, expected) {
  n_trials <- nrow(expected$n_outcome)
  n_doses <- dim(trials$n_outcome)[2L]
  counts <- matrix(aperm(trials$n_outcome, c(1L, 3L, 2L)), nrow = n_trials)
  expect_equal(counts, expected$n_outcome)
  expect_equal(unname(trials$eliminated), expected$eliminated)
  expect_equal(trials$stage2, expected$stage2)
  expect_equal(trials$cohorts_used, expected$cohorts_used)
  expect_equal(trials$stop_reason, expected$stop_reason)
  expect_equal(dim(trials$n_pts), c(n_trials, n_doses))
}

uboin_sim_args <- function(scenario, ...) {
  sc <- zhou2019_scenario(scenario)
  args <- list(tox_upper = 0.30, resp_lower = 0.20, p_true = sc$p_true,
               p_resp = sc$p_resp, n_cohort = 18, cohort_size = 3,
               resp_assoc = 0.2)
  given <- list(...)
  args[names(given)] <- given
  args
}

test_that("uboin_simulate returns a well formed object", {
  trials <- do.call(uboin_simulate,
                    uboin_sim_args(1, n_trials = 200, seed = 1))

  expect_s3_class(trials, "uboin_trials")
  expect_equal(dim(trials$n_outcome), c(200L, 5L, 4L))
  expect_equal(dimnames(trials$n_outcome)[[3L]],
               c("no_resp_dlt", "no_resp_no_dlt", "resp_dlt", "resp_no_dlt"))
  expect_equal(colnames(trials$n_pts), paste0("DL", 1:5))
  expect_type(trials$n_pts, "integer")
  expect_type(trials$eliminated, "logical")
  expect_length(trials$stage2, 200L)
  expect_equal(trials$design, "U-BOIN")
  expect_equal(trials$settings$target, 0.25)
  expect_equal(trials$settings$n_stop, 54L)
})

test_that("the counts of the trials are consistent", {
  trials <- do.call(uboin_simulate,
                    uboin_sim_args(4, n_trials = 500, seed = 2))
  totals <- rowSums(trials$n_pts)
  expect_identical(trials$n_pts, apply(trials$n_outcome, c(1L, 2L), sum))
  expect_identical(trials$n_tox,
                   trials$n_outcome[, , 1L] + trials$n_outcome[, , 3L])
  expect_identical(trials$n_resp,
                   trials$n_outcome[, , 3L] + trials$n_outcome[, , 4L])
  expect_equal(totals, 3 * trials$cohorts_used)
  expect_true(all(totals <= 54))
  expect_true(all(trials$stop_reason %in%
                    c("max_sample_size", "n_stop", "lowest_dose_eliminated",
                      "no_admissible_dose")))
  expect_true(all(totals[trials$stop_reason == "max_sample_size"] == 54))
  # A trial ended by its admissible doses has reached stage II.
  expect_true(all(trials$stage2[trials$stop_reason == "no_admissible_dose"]))
  # A trial stopped in stage I has the lowest dose and all doses eliminated.
  low <- trials$stop_reason == "lowest_dose_eliminated"
  expect_true(all(trials$eliminated[low, ]))
  expect_false(any(trials$stage2[low]))
})

test_that("the trials agree with a separate implementation in Python", {
  uboin_check_trials(
    do.call(uboin_simulate, uboin_sim_args(1, n_trials = 5, seed = 1)),
    uboin_expected_scenario1()
  )
  uboin_check_trials(
    do.call(uboin_simulate, uboin_sim_args(8, n_trials = 8, seed = 2)),
    uboin_expected_scenario8()
  )
  uboin_check_trials(
    uboin_simulate(tox_upper = 0.30, resp_lower = 0.20,
                   p_true = c(0.45, 0.55, 0.65, 0.75, 0.85),
                   p_resp = c(0.30, 0.40, 0.50, 0.60, 0.70), n_cohort = 18,
                   cohort_size = 3, n_trials = 6, seed = 4),
    uboin_expected_toxic()
  )
  uboin_check_trials(
    do.call(uboin_simulate,
            uboin_sim_args(1, n_cohort = 15, cohort_size = 2, n_stage1 = 6,
                           n_stop = 12, start_dose = 2, n_trials = 5,
                           seed = 3)),
    uboin_expected_options()
  )
})

test_that("p_outcome gives the trials of the marginal rates", {
  sc <- zhou2019_scenario(4)
  probs <- uboin_outcome_probs(sc$p_true, sc$p_resp, 0.2)
  expect_equal(rowSums(probs), rep(1, 5))
  expect_equal(probs[, 1] + probs[, 3], sc$p_true)
  expect_equal(probs[, 3] + probs[, 4], sc$p_resp)

  a <- do.call(uboin_simulate, uboin_sim_args(4, n_trials = 300, seed = 5))
  b <- uboin_simulate(tox_upper = 0.30, resp_lower = 0.20, p_outcome = probs,
                      n_cohort = 18, cohort_size = 3, n_trials = 300,
                      seed = 5)
  expect_identical(a$n_outcome, b$n_outcome)
  expect_equal(b$settings$p_true, sc$p_true)
})

test_that("the seed makes the trials reproducible and is restored", {
  set.seed(99)
  before <- .Random.seed
  a <- do.call(uboin_simulate, uboin_sim_args(1, n_trials = 50, seed = 7))
  expect_identical(.Random.seed, before)
  b <- do.call(uboin_simulate, uboin_sim_args(1, n_trials = 50, seed = 7))
  expect_identical(a$n_outcome, b$n_outcome)
})

test_that("the boundaries of stage I are those of Table S2 of the article", {
  trials <- do.call(uboin_simulate, uboin_sim_args(1, n_trials = 2, seed = 1))
  b <- trials$boundary
  expect_equal(b$b_esc[1:24], c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
                                3, 3, 3, 3, 3, 4, 4, 4, 4))
  expect_equal(b$b_deesc[1:24], c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5,
                                  5, 5, 6, 6, 6, 6, 7, 7, 7, 8))
  expect_equal(b$b_elim[1:24], c(NA, NA, 3, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8,
                                 8, 8, 9, 9, 9, 10, 10, 11, 11, 11))
})

test_that("uboin_simulate validates its arguments", {
  sc <- zhou2019_scenario(1)
  expect_error(
    uboin_simulate(tox_upper = 0.3, resp_lower = 0.2, p_true = sc$p_true,
                   n_cohort = 18, cohort_size = 3, n_trials = 2),
    "give 'p_true' and 'p_resp', or 'p_outcome'"
  )
  expect_error(
    uboin_simulate(tox_upper = 0.3, resp_lower = 0.2, p_true = sc$p_true,
                   p_resp = sc$p_resp, p_outcome = diag(4), n_cohort = 18,
                   cohort_size = 3, n_trials = 2),
    "not both"
  )
  expect_error(
    uboin_simulate(tox_upper = 0.3, resp_lower = 0.2,
                   p_outcome = matrix(0.3, 2, 4), n_cohort = 18,
                   cohort_size = 3, n_trials = 2),
    "each row of 'p_outcome' must add up to 1"
  )
  expect_error(
    uboin_simulate(tox_upper = 0.3, resp_lower = 0.2, p_true = sc$p_true,
                   p_resp = sc$p_resp[1:4], n_cohort = 18, cohort_size = 3,
                   n_trials = 2),
    "'p_resp' must hold one response probability"
  )
  expect_error(
    do.call(uboin_simulate, uboin_sim_args(1, n_stop = 9, n_trials = 2)),
    "'n_stop' must not be smaller than 'n_stage1'"
  )
  expect_error(
    do.call(uboin_simulate, uboin_sim_args(1, start_dose = 6, n_trials = 2)),
    "'start_dose' must not exceed the number of doses"
  )
  expect_error(
    do.call(uboin_simulate, uboin_sim_args(1, target = 0.35, n_trials = 2)),
    "'target' must not exceed 'tox_upper'"
  )
})
