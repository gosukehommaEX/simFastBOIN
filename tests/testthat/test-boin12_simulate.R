sc1 <- list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
            p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))

# Trials of the fixture boin12_rules.csv, computed with the Python
# implementation of the rules fed with R's random numbers (see
# Validation/tite_boin12/make_rules_fixture.py of the development folder):
# the scenarios below, with set.seed(20261117), split by the combination of
# tie_break, no_candidate and explore_from. Scenario F is Scenario 16 of Lin
# et al. (2020).
rules_scenarios <- list(
  E = list(p_true = c(0.10, 0.20, 0.40, 0.55, 0.70),
           p_resp = c(0.05, 0.10, 0.30, 0.35, 0.40), start_dose = 2,
           n_trials = 40),
  F = list(p_true = c(0.25, 0.55, 0.67, 0.87, 0.95),
           p_resp = c(0.01, 0.55, 0.55, 0.55, 0.55), start_dose = 1,
           n_trials = 60)
)

rules_fixture <- function(design) {
  fx <- utils::read.csv(test_path("fixtures", "boin12_rules.csv"),
                        stringsAsFactors = FALSE)
  fx <- fx[fx$design == design, , drop = FALSE]
  split(fx, paste(fx$scenario, fx$tie_break, fx$no_candidate,
                  fx$explore_from))
}

rules_counts <- function(fx, prefix) {
  unname(as.matrix(fx[, paste0(prefix, 1:5)]))
}

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
  expect_equal(trials$settings$tie_break, "higher")
  expect_equal(trials$settings$no_candidate, "nearest")
  expect_equal(trials$settings$explore_from, "any")
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

test_that("no_candidate = \"stop\" ends a trial with no admissible candidate", {
  # Nine patients without response exclude dose 1 by efficacy and the
  # exploration rule escalates; one DLT in the one patient at dose 2 then
  # calls for a de-escalation with no admissible dose below. "nearest" stays
  # at dose 2, where the next DLT excludes every dose; dose 3 is never tried.
  args <- list(target = 0.35, target_resp = 0.25, p_true = c(0, 1, 1),
               p_resp = c(0, 1, 1), n_cohort = 4, cohort_size = c(9, 1),
               n_trials = 3, seed = 1)
  near <- do.call(boin12_simulate, args)
  stop <- do.call(boin12_simulate, c(args, list(no_candidate = "stop")))
  for (i in 1:3) {
    expect_equal(unname(near$n_pts[i, ]), c(9L, 2L, 0L))
    expect_equal(unname(stop$n_pts[i, ]), c(9L, 1L, 0L))
  }
  expect_true(all(near$stop_reason == "no_admissible_dose"))
  expect_true(all(stop$stop_reason == "no_admissible_candidate"))
  expect_true(all(stop$early_stop))
  expect_equal(stop$cohorts_used, rep(2L, 3))
  expect_equal(stop$settings$no_candidate, "stop")
  # Doses 2 and 3 are excluded by efficacy after nine patients each, so that
  # at the highest dose neither candidate is admissible; "nearest" goes to
  # dose 1, untried, and treats nine patients there.
  args <- list(target = 0.35, target_resp = 0.25, p_true = rep(0, 3),
               p_resp = rep(0, 3), n_cohort = 6, cohort_size = c(9, 9, 3),
               start_dose = 2, n_trials = 3, seed = 1)
  near <- do.call(boin12_simulate, args)
  stop <- do.call(boin12_simulate, c(args, list(no_candidate = "stop")))
  for (i in 1:3) {
    expect_equal(unname(near$n_pts[i, ]), c(9L, 9L, 9L))
    expect_equal(unname(stop$n_pts[i, ]), c(0L, 9L, 9L))
  }
  expect_equal(near$cohorts_used, rep(5L, 3))
  expect_true(all(near$stop_reason == "no_admissible_dose"))
  expect_true(all(stop$stop_reason == "no_admissible_candidate"))
})

test_that("tie_break = \"random\" draws one of the tied doses", {
  # After three patients without DLT and response at dose 2, the untried
  # doses 1 and 3 are tied (desirability 0.295, against 0.113 for dose 2).
  # With seed 1 the second cohort goes to dose 1 in 213 of the 400 trials
  # (Python with R's random numbers, make_rules_fixture.py).
  args <- list(target = 0.35, target_resp = 0.25, p_true = rep(0, 3),
               p_resp = rep(0, 3), n_cohort = 2, cohort_size = 3,
               start_dose = 2, n_trials = 400, seed = 1)
  high <- do.call(boin12_simulate, args)
  rand <- do.call(boin12_simulate, c(args, list(tie_break = "random")))
  expect_true(all(high$n_pts[, 1] == 0L))
  expect_true(all(high$n_pts[, 3] == 3L))
  expect_true(all(rand$n_pts[, 1] + rand$n_pts[, 3] == 3L))
  expect_identical(sum(rand$n_pts[, 1] == 3L), 213L)
  expect_equal(rand$settings$tie_break, "random")
})

test_that("explore_from = \"admissible\" explores from admissible doses only", {
  # Nine patients who all respond at dose 1, then the exploration rule
  # escalates; nine patients without response exclude dose 2 by efficacy.
  # "any" still explores dose 3, "admissible" goes back to dose 1.
  args <- list(target = 0.35, target_resp = 0.25, p_true = rep(0, 3),
               p_resp = c(1, 0, 0), n_cohort = 5,
               cohort_size = c(3, 3, 3, 9, 3), n_trials = 3, seed = 1)
  from_any <- do.call(boin12_simulate, args)
  from_adm <- do.call(boin12_simulate,
                      c(args, list(explore_from = "admissible")))
  for (i in 1:3) {
    expect_equal(unname(from_any$n_pts[i, ]), c(9L, 9L, 3L))
    expect_equal(unname(from_adm$n_pts[i, ]), c(12L, 9L, 0L))
  }
  expect_true(all(from_adm$stop_reason == "max_cohorts"))
  expect_equal(from_adm$settings$explore_from, "admissible")
})

test_that("the rules for the next dose agree with a separate implementation", {
  # The defaults give the trials of the rules before the three arguments were
  # added; in scenario F explore_from changes some of the trials.
  groups <- rules_fixture("BOIN12")
  expect_length(groups, 6L)
  n_pts <- vector("list", length(groups))
  for (k in seq_along(groups)) {
    fx <- groups[[k]]
    sc <- rules_scenarios[[fx$scenario[1]]]
    tr <- boin12_simulate(target = 0.35, target_resp = 0.25,
                          p_true = sc$p_true, p_resp = sc$p_resp,
                          n_cohort = 12, cohort_size = 3,
                          n_trials = sc$n_trials, start_dose = sc$start_dose,
                          tie_break = fx$tie_break[1],
                          no_candidate = fx$no_candidate[1],
                          explore_from = fx$explore_from[1], seed = 20261117)
    expect_identical(unname(tr$n_pts), rules_counts(fx, "n"))
    expect_identical(unname(tr$n_tox), rules_counts(fx, "y"))
    expect_identical(unname(tr$n_resp), rules_counts(fx, "x"))
    expect_identical(unname(tr$n_tox_resp), rules_counts(fx, "z"))
    expect_identical(tr$cohorts_used, fx$cohorts_used)
    expect_identical(tr$stop_reason, fx$stop_reason)
    expect_identical(tr$early_stop,
                     tr$stop_reason %in% c("no_admissible_dose",
                                           "no_admissible_candidate"))
    n_pts[[k]] <- tr$n_pts
  }
  for (a in 1:5) {
    for (b in (a + 1):6) {
      expect_false(identical(n_pts[[a]], n_pts[[b]]))
    }
  }
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
                       modifyList(args, list(tie_break = "lower"))),
               "should be one of")
  expect_error(do.call(boin12_simulate,
                       modifyList(args, list(no_candidate = "wait"))),
               "should be one of")
  expect_error(do.call(boin12_simulate,
                       modifyList(args, list(explore_from = "all"))),
               "should be one of")
  expect_error(do.call(boin12_simulate,
                       modifyList(args, list(resp_cor = 0.9,
                                             resp_cor_type = "binary"))),
               "attainable")
})
