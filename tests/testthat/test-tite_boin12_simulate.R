sc1 <- list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
            p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))
sc2 <- list(p_true = c(0.05, 0.10, 0.15, 0.20, 0.25),
            p_resp = c(0.05, 0.08, 0.10, 0.12, 0.15))

# Trials of the fixture boin12_rules.csv, see test-boin12_simulate.R.
rules_fixture <- function(design, tie_break, no_candidate) {
  fx <- utils::read.csv(test_path("fixtures", "boin12_rules.csv"),
                        stringsAsFactors = FALSE)
  fx[fx$design == design & fx$tie_break == tie_break &
       fx$no_candidate == no_candidate, , drop = FALSE]
}

rules_counts <- function(fx, prefix) {
  unname(as.matrix(fx[, paste0(prefix, 1:5)]))
}

tite_args <- function(...) {
  utils::modifyList(
    list(target = 0.35, target_resp = 0.25, p_true = sc1$p_true,
         p_resp = sc1$p_resp, n_cohort = 12, cohort_size = 3, window = 45,
         resp_window = 60, accrual_rate = 0.1, n_trials = 200, seed = 1),
    list(...)
  )
}

# Interim quantities of the patients of Table 2 of Zhou et al. (2022)
# (TITE-BOIN12 columns) at day `day`; utilities 100, 40, 60 and 0, windows of
# 45 and 60 days, upper DLT limit 0.35 and lower response limit 0.25.
zhou_interim <- function(day, pending) {
  dose <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 3, 3, 3, 3, 3, 3)
  yt <- c(0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1)
  ye <- c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0)
  enter <- c(1, 11, 21, 72, 82, 92, 143, 153, 163, 214, 224, 234, 285, 295,
             305, 316, 326, 336)
  xt <- c(0, 0, 0, 0, 30, 0, 0, 25, 0, 0, 40, 35, 25, 0, 0, 0, 0, 30)
  xe <- c(0, 0, 0, 40, 0, 0, 40, 0, 0, 0, 0, 0, 30, 55, 0, 50, 50, 0)
  des <- boin12_check_design(0.35, 0.25)
  tite_boin12_interim_cpp(
    dose = as.integer(dose), enter = enter,
    end_tox = enter + ifelse(yt == 1, xt, 45),
    end_resp = enter + ifelse(ye == 1, xe, 60),
    dlt = as.integer(yt), resp = as.integer(ye), time = day, n_doses = 5L,
    window = 45, resp_window = 60, psi00 = 40, psi11 = 60, alpha = 1,
    beta = 1, u_benchmark = des$u_benchmark / 100, target = 0.35,
    target_resp = 0.25, cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
    pending = if (pending == "marginal") 1L else 0L
  )
}

test_that("tite_boin12_simulate returns a well formed object", {
  trials <- do.call(tite_boin12_simulate, tite_args())
  expect_s3_class(trials, "boin12_trials")
  for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
    expect_equal(dim(trials[[nm]]), c(200L, 5L))
    expect_type(trials[[nm]], "integer")
  }
  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(trials$n_resp <= trials$n_pts))
  expect_true(all(trials$n_tox_resp <= pmin(trials$n_tox, trials$n_resp)))
  expect_true(all(rowSums(trials$n_pts) <= 36L))
  expect_true(all(trials$stop_reason %in% boin12_stop_reasons()))
  expect_identical(trials$early_stop,
                   trials$stop_reason == "no_admissible_dose")
  expect_length(trials$duration, 200L)
  expect_true(all(trials$duration > 0))
  expect_equal(trials$design, "TITE-BOIN12")
  expect_equal(trials$settings$pending_estimate, "conditional")
  expect_equal(trials$settings$min_completed_share, 0.5)
  expect_equal(trials$settings$tie_break, "higher")
  expect_equal(trials$settings$no_candidate, "nearest")
})

test_that("waiting for every outcome gives the BOIN12 trials", {
  # With min_completed_share = 1 and fixed arrivals no outcome is pending at a
  # decision and no variate is drawn for the arrivals, so the trials are
  # those of boin12_simulate() with the same seed.
  settings <- list(
    list(sc = sc1, extra = list()),
    list(sc = sc2, extra = list(n_earlystop = 9, resp_cor = 0.3)),
    list(sc = sc1, extra = list(resp_cor = 0.1, resp_cor_type = "binary",
                                start_dose = 3)),
    list(sc = sc2, extra = list(psi00 = 30, psi11 = 50, n_earlystop = 12,
                                n_earlystop_rule = "simple")),
    list(sc = sc1, extra = list(tie_break = "random", no_candidate = "stop",
                                start_dose = 3)),
    list(sc = sc2, extra = list(tie_break = "random", no_candidate = "stop",
                                n_earlystop = 9))
  )
  for (s in settings) {
    common <- c(list(target = 0.35, target_resp = 0.25,
                     p_true = s$sc$p_true, p_resp = s$sc$p_resp,
                     n_cohort = 12, cohort_size = 3, n_trials = 300,
                     seed = 11), s$extra)
    a <- do.call(boin12_simulate, common)
    b <- do.call(tite_boin12_simulate,
                 c(common, list(window = 45, resp_window = 60,
                                accrual_rate = 0.1, accrual = "fixed",
                                min_completed_share = 1)))
    expect_identical(b$n_pts, a$n_pts)
    expect_identical(b$n_tox, a$n_tox)
    expect_identical(b$n_resp, a$n_resp)
    expect_identical(b$n_tox_resp, a$n_tox_resp)
    expect_identical(b$cohorts_used, a$cohorts_used)
    expect_identical(b$stop_reason, a$stop_reason)
  }
})

test_that("the rules for the next dose agree with a separate implementation", {
  # The fixture of test-boin12_simulate.R, with windows of 45 and 60 days and
  # 0.1 patients a day; the durations include the trials stopped by
  # no_candidate = "stop" at the time of the decision.
  combos <- list(c("higher", "nearest"), c("random", "nearest"),
                 c("higher", "stop"), c("random", "stop"))
  for (cmb in combos) {
    tr <- do.call(tite_boin12_simulate,
                  tite_args(p_true = c(0.10, 0.20, 0.40, 0.55, 0.70),
                            p_resp = c(0.05, 0.10, 0.30, 0.35, 0.40),
                            n_trials = 40, start_dose = 2,
                            tie_break = cmb[1], no_candidate = cmb[2],
                            seed = 20261117))
    fx <- rules_fixture("TITE-BOIN12", cmb[1], cmb[2])
    expect_identical(unname(tr$n_pts), rules_counts(fx, "n"))
    expect_identical(unname(tr$n_tox), rules_counts(fx, "y"))
    expect_identical(unname(tr$n_resp), rules_counts(fx, "x"))
    expect_identical(unname(tr$n_tox_resp), rules_counts(fx, "z"))
    expect_identical(tr$cohorts_used, fx$cohorts_used)
    expect_identical(tr$stop_reason, fx$stop_reason)
    expect_equal(tr$duration, fx$duration, tolerance = 1e-10)
    expect_identical(tr$early_stop,
                     tr$stop_reason %in% c("no_admissible_dose",
                                           "no_admissible_candidate"))
  }
})

test_that("no_candidate = \"stop\" ends a trial as in BOIN12", {
  # The two cases of test-boin12_simulate.R, waiting for every outcome.
  timing <- list(window = 3, resp_window = 2, accrual_rate = 1,
                 accrual = "fixed", event_time = "uniform",
                 min_completed_share = 1, no_candidate = "stop")
  one <- do.call(tite_boin12_simulate,
                 c(list(target = 0.35, target_resp = 0.25,
                        p_true = c(0, 1, 1), p_resp = c(0, 1, 1),
                        n_cohort = 4, cohort_size = c(9, 1), n_trials = 3,
                        seed = 1), timing))
  two <- do.call(tite_boin12_simulate,
                 c(list(target = 0.35, target_resp = 0.25,
                        p_true = rep(0, 3), p_resp = rep(0, 3), n_cohort = 6,
                        cohort_size = c(9, 9, 3), start_dose = 2,
                        n_trials = 3, seed = 1), timing))
  for (i in 1:3) {
    expect_equal(unname(one$n_pts[i, ]), c(9L, 1L, 0L))
    expect_equal(unname(two$n_pts[i, ]), c(0L, 9L, 9L))
  }
  expect_true(all(one$stop_reason == "no_admissible_candidate"))
  expect_true(all(two$stop_reason == "no_admissible_candidate"))
  expect_true(all(one$early_stop & two$early_stop))
  # Cohort 1 arrives at days 0 to 8 and is complete at day 11, cohort 2
  # arrives at days 11 to 19 and the trial stops at day 22.
  expect_equal(two$duration, rep(22, 3))
  expect_equal(two$settings$no_candidate, "stop")
})

test_that("the interim desirabilities of Table S.3 of Zhou et al. (2022)", {
  # Expected values computed independently with Python (scipy). Analysis 2
  # (day 142), dose 2: one patient has a pending response, w = 50 / 60, and
  # pi_E = 1 / (2 + 50 / 60) = 6 / 17. Replacing the pending response by
  # pi_E gives 0.211 as in the table; Equation (11) gives 0.168.
  marg <- zhou_interim(142, "marginal")
  cond <- zhou_interim(142, "conditional")
  expect_equal(marg$n, c(3L, 3L, 0L, 0L, 0L))
  expect_equal(marg$pi_tox[2], 1 / 3, tolerance = 1e-12)
  expect_equal(marg$pi_resp[2], 6 / 17, tolerance = 1e-12)
  expect_equal(marg$quasi_events[2], 1.6117647058823528, tolerance = 1e-12)
  expect_equal(cond$quasi_events[2], 1.45, tolerance = 1e-12)
  expect_equal(marg$desirability[2], 0.21097652928898158, tolerance = 1e-10)
  expect_equal(cond$desirability[2], 0.1678137799050769, tolerance = 1e-10)
  expect_equal(round(marg$desirability[1:2], 3), c(0.113, 0.211))
  # Untried doses: Pr(Beta(1, 1) > 0.705)
  expect_equal(marg$desirability[3:5], rep(0.295, 3), tolerance = 1e-12)
  expect_true(all(marg$admissible))
  # Analysis 5 (day 315), dose 3: the table reports 0.193, which neither
  # replacement gives from Table 2 (see inst/reproduce/reproduce_Zhou2022.R).
  marg5 <- zhou_interim(315, "marginal")
  cond5 <- zhou_interim(315, "conditional")
  expect_equal(marg5$pi_tox[3], 3 / 7, tolerance = 1e-12)
  expect_equal(marg5$pi_resp[3], 4 / 9, tolerance = 1e-12)
  expect_equal(marg5$desirability[3], 0.11821463548407772, tolerance = 1e-10)
  expect_equal(cond5$desirability[3], 0.11720076710331795, tolerance = 1e-10)
  # Analysis 6 (day 396): every outcome is observed, as in BOIN12
  final <- zhou_interim(396, "conditional")
  expect_equal(round(final$desirability[1:4], 3),
               c(0.113, 0.156, 0.208, 0.022))
  expect_equal(final$desirability, zhou_interim(396, "marginal")$desirability)
})

test_that("a late response does not exclude a dose", {
  # Nine patients at dose 1: five completed both windows without event, four
  # arrived 30 days before the decision. Admissibility uses
  # Beta(1, 1 + 5 + 4 * 0.5) for the response, Pr(p_E < 0.25) = 0.8999 < 0.90,
  # whereas counting the four as non-responders would give 0.9437.
  enter <- c(rep(0, 5), rep(70, 4))
  res <- tite_boin12_interim_cpp(
    dose = rep(1L, 9), enter = enter, end_tox = enter + 45,
    end_resp = enter + 60, dlt = rep(0L, 9), resp = rep(0L, 9), time = 100,
    n_doses = 2L, window = 45, resp_window = 60, psi00 = 40, psi11 = 60,
    alpha = 1, beta = 1, u_benchmark = 0.705, target = 0.35,
    target_resp = 0.25, cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
    pending = 0L
  )
  expect_true(res$admissible[1])
  expect_equal(res$pi_resp[1], 0)
  expect_equal(res$quasi_events[1], 3.6, tolerance = 1e-12)
  expect_equal(res$desirability[1], 0.024903959669505104, tolerance = 1e-10)
  expect_equal(stats::pbeta(0.25, 1, 1 + 5 + 2), 0.8998870849609375,
               tolerance = 1e-12)
})

test_that("the duration counts the waits and the windows", {
  # No DLT and a response in every patient, one arrival per unit of time,
  # DLT window 3 and response window 2. Without waiting the 36 patients
  # arrive at times 0 to 35 and the last DLT window ends at 38; waiting for
  # every outcome, cohort k arrives at 5 k and the trial ends at 57 + 3.
  args <- list(target = 0.35, target_resp = 0.25, p_true = rep(0, 3),
               p_resp = rep(1, 3), n_cohort = 12, cohort_size = 3,
               window = 3, resp_window = 2, accrual_rate = 1,
               accrual = "fixed", event_time = "uniform", n_trials = 5,
               seed = 2)
  a <- do.call(tite_boin12_simulate, c(args, list(min_completed_share = 0)))
  b <- do.call(tite_boin12_simulate, c(args, list(min_completed_share = 1)))
  expect_equal(a$duration, rep(38, 5))
  expect_equal(b$duration, rep(60, 5))
  for (i in 1:5) {
    expect_equal(unname(a$n_pts[i, ]), c(33L, 3L, 0L))
  }
  expect_true(all(a$stop_reason == "max_cohorts"))
})

test_that("a toxic lowest dose stops the trial", {
  trials <- tite_boin12_simulate(target = 0.35, target_resp = 0.25,
                                 p_true = rep(1, 4), p_resp = rep(0.5, 4),
                                 n_cohort = 10, cohort_size = 3, window = 30,
                                 accrual_rate = 0.1, event_time = "uniform",
                                 n_trials = 20, seed = 3)
  expect_true(all(trials$stop_reason == "no_admissible_dose"))
  expect_true(all(trials$early_stop))
  expect_true(all(trials$n_pts[, 2:4] == 0L))
  expect_true(all(trials$n_tox == trials$n_pts))
})

test_that("the replacement of a pending outcome changes the decisions", {
  a <- do.call(tite_boin12_simulate, tite_args(n_trials = 300))
  b <- do.call(tite_boin12_simulate,
               tite_args(n_trials = 300, pending_estimate = "marginal"))
  expect_false(identical(a$n_pts, b$n_pts))
  expect_equal(b$settings$pending_estimate, "marginal")
})

test_that("the seed reproduces the trials and the session state is kept", {
  args <- tite_args(resp_cor = 0.3, accrual = "exponential", n_trials = 100)
  set.seed(99)
  before <- .Random.seed
  a <- do.call(tite_boin12_simulate, args)
  expect_identical(.Random.seed, before)
  b <- do.call(tite_boin12_simulate, args)
  expect_identical(a$n_pts, b$n_pts)
  expect_identical(a$n_tox_resp, b$n_tox_resp)
  expect_identical(a$duration, b$duration)
})

test_that("invalid arguments are rejected", {
  args <- tite_args(n_trials = 10)
  expect_error(do.call(tite_boin12_simulate,
                       modifyList(args, list(min_completed_share = 2))),
               "min_completed_share")
  expect_error(do.call(tite_boin12_simulate,
                       modifyList(args, list(window = 0))), "window")
  expect_error(do.call(tite_boin12_simulate,
                       modifyList(args, list(accrual_rate = -1))),
               "accrual_rate")
  expect_error(do.call(tite_boin12_simulate,
                       modifyList(args, list(p_true = c(sc1$p_true[1:4],
                                                        1)))),
               "below 1")
  expect_error(do.call(tite_boin12_simulate,
                       modifyList(args, list(pending_estimate = "x"))),
               "should be one of")
  expect_error(do.call(tite_boin12_simulate,
                       modifyList(args, list(tie_break = "lower"))),
               "should be one of")
  expect_error(do.call(tite_boin12_simulate,
                       modifyList(args, list(no_candidate = "wait"))),
               "should be one of")
  expect_error(do.call(tite_boin12_simulate,
                       modifyList(args, list(p_resp = c(0.2, 0.3)))),
               "p_resp")
  expect_error(do.call(tite_boin12_simulate,
                       modifyList(args, list(start_dose = 6))),
               "start_dose")
  expect_error(do.call(tite_boin12_simulate,
                       modifyList(args, list(n_earlystop = 0))),
               "n_earlystop")
})
