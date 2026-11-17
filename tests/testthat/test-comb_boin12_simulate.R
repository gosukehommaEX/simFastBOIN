# Scenario 16 of Lu et al. (2025), with their utility (wT = 0.11)
sc16 <- list(p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), 2),
             p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), 2))
lu_psi00 <- 100 * 0.11 / 1.11
lu_psi11 <- 100 / 1.11

# Trials of the fixture comb_boin12_rules.csv, computed with the Python
# implementation of the rules fed with R's random numbers and checked
# against the C++ core (see Validation/comb_boin12/make_cb12_fixture.py of
# the development folder): scenario G below with set.seed(20261117), eight
# cohorts of three, the utility of Lu et al. (2025), split by tie_break,
# no_candidate and extrasafe.
scenario_g <- list(
  p_true = matrix(c(0.05, 0.10, 0.30, 0.10, 0.30, 0.50, 0.30, 0.50, 0.70),
                  3),
  p_resp = matrix(c(0.05, 0.10, 0.40, 0.10, 0.40, 0.50, 0.40, 0.50, 0.55),
                  3)
)

test_that("comb_boin12_simulate returns a well formed object", {
  trials <- comb_boin12_simulate(target = 0.33, target_resp = 0.2,
                                 p_true = sc16$p_true, p_resp = sc16$p_resp,
                                 n_cohort = 14, cohort_size = 3,
                                 n_trials = 200, seed = 1)
  expect_s3_class(trials, "comb_boin12_trials")
  for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
    expect_equal(dim(trials[[nm]]), c(200L, 2L, 3L))
    expect_type(trials[[nm]], "integer")
    expect_equal(dimnames(trials[[nm]])[[2L]], c("A1", "A2"))
    expect_equal(dimnames(trials[[nm]])[[3L]], c("B1", "B2", "B3"))
  }
  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(trials$n_resp <= trials$n_pts))
  expect_true(all(trials$n_tox_resp <= pmin(trials$n_tox, trials$n_resp)))
  expect_true(all(apply(trials$n_pts, 1L, sum) <= 42L))
  expect_true(all(trials$stop_reason %in% boin12_stop_reasons()))
  expect_identical(trials$early_stop,
                   trials$stop_reason == "no_admissible_dose")
  expect_equal(trials$settings$max_total_pts, 42)
  expect_equal(trials$settings$tie_break, "higher")
  expect_equal(trials$settings$no_candidate, "nearest")
  expect_false(trials$settings$extrasafe)
  expect_equal(trials$design, "Comb-BOIN12")
  lam <- boin_lambda(0.33)
  expect_equal(trials$boundary$lambda_e, lam$lambda_e)
  expect_equal(trials$boundary$lambda_d, lam$lambda_d)
})

test_that("one row or one column of combinations gives the BOIN12 trials", {
  # Comb-BOIN12 has no exploration rule; with n_explore above the sample
  # size, BOIN12 follows the same rules along a single drug.
  p_true <- c(0.05, 0.12, 0.27, 0.35, 0.50)
  p_resp <- c(0.20, 0.35, 0.36, 0.37, 0.38)
  settings <- expand.grid(tie_break = c("higher", "random"),
                          no_candidate = c("nearest", "stop"),
                          stringsAsFactors = FALSE)
  for (i in seq_len(nrow(settings))) {
    args <- list(target = 0.35, target_resp = 0.25, n_cohort = 12,
                 cohort_size = 3, n_trials = 200, n_earlystop = 9,
                 tie_break = settings$tie_break[i],
                 no_candidate = settings$no_candidate[i], seed = 11)
    single <- do.call(boin12_simulate,
                      c(args, list(p_true = p_true, p_resp = p_resp,
                                   n_explore = 100)))
    for (shape in c("row", "column")) {
      dims <- if (shape == "row") c(1L, 5L) else c(5L, 1L)
      comb <- do.call(comb_boin12_simulate,
                      c(args, list(p_true = matrix(p_true, dims[1L]),
                                   p_resp = matrix(p_resp, dims[1L]))))
      for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
        expect_identical(unname(matrix(comb[[nm]], 200L)),
                         unname(single[[nm]]))
      }
      expect_identical(comb$cohorts_used, single$cohorts_used)
      expect_identical(comb$stop_reason, single$stop_reason)
    }
  }
})

test_that("ties keep the combinations with the higher sum of the levels", {
  # Three patients without DLT and response at (2, 2): its four untried
  # neighbors are tied. "higher" keeps (3, 2) and (2, 3) and draws one of
  # them; "random" draws one of the four. With seed 1 (Python with R's
  # random numbers, test_expect_cb12.py): 208 and 192 trials for "higher",
  # 213 trials to (1, 2) or (2, 1), 103 to (3, 2) and 84 to (2, 3) for
  # "random".
  args <- list(target = 0.35, target_resp = 0.25,
               p_true = matrix(0, 3, 3), p_resp = matrix(0, 3, 3),
               n_cohort = 2, cohort_size = 3, start_dose = c(2, 2),
               n_trials = 400, seed = 1)
  high <- do.call(comb_boin12_simulate, args)
  rand <- do.call(comb_boin12_simulate, c(args, list(tie_break = "random")))
  expect_true(all(high$n_pts[, 1, 2] + high$n_pts[, 2, 1] == 0L))
  expect_identical(sum(high$n_pts[, 3, 2] == 3L), 208L)
  expect_identical(sum(high$n_pts[, 2, 3] == 3L), 192L)
  expect_identical(sum(rand$n_pts[, 1, 2] + rand$n_pts[, 2, 1] == 3L), 213L)
  expect_identical(sum(rand$n_pts[, 3, 2] == 3L), 103L)
  expect_identical(sum(rand$n_pts[, 2, 3] == 3L), 84L)
  expect_equal(rand$settings$tie_break, "random")
})

test_that("no_candidate = \"stop\" ends a trial with no admissible candidate", {
  # Nine patients without response exclude (1, 2), (2, 1) and (3, 1) by
  # efficacy; nine DLTs at (2, 2) exclude it with (3, 2). At (3, 1) none of
  # the candidates is admissible: "nearest" goes to (1, 1), "stop" ends the
  # trial (test_expect_cb12.py).
  args <- list(target = 0.35, target_resp = 0.25,
               p_true = matrix(c(0, 0, 0, 0, 1, 1), 3),
               p_resp = matrix(0, 3, 2), n_cohort = 5,
               cohort_size = c(9, 9, 9, 9, 3), start_dose = c(1, 2),
               n_trials = 3, seed = 1)
  near <- do.call(comb_boin12_simulate, args)
  stop <- do.call(comb_boin12_simulate, c(args, list(no_candidate = "stop")))
  for (i in 1:3) {
    expect_equal(unname(near$n_pts[i, , ]), matrix(c(3, 9, 9, 9, 9, 0), 3))
    expect_equal(unname(stop$n_pts[i, , ]), matrix(c(0, 9, 9, 9, 9, 0), 3))
  }
  expect_true(all(near$stop_reason == "max_cohorts"))
  expect_true(all(stop$stop_reason == "no_admissible_candidate"))
  expect_true(all(stop$early_stop))
  expect_equal(stop$cohorts_used, rep(4L, 3))
  expect_equal(stop$settings$no_candidate, "stop")
})

test_that("extrasafe lowers the cutoff at the lowest combination", {
  # One DLT in one patient: Pr(pT > 0.35) = 0.8775, above 0.95 - 0.10 but
  # not above 0.95. Without extrasafe the second patient stops the trial.
  args <- list(target = 0.35, target_resp = 0.25, p_true = matrix(1, 2, 2),
               p_resp = matrix(0.5, 2, 2), n_cohort = 4, cohort_size = 1,
               n_trials = 3, seed = 1)
  plain <- do.call(comb_boin12_simulate, args)
  safe <- do.call(comb_boin12_simulate,
                  c(args, list(extrasafe = TRUE, offset = 0.10)))
  expect_true(all(plain$n_pts[, 1, 1] == 2L))
  expect_true(all(safe$n_pts[, 1, 1] == 1L))
  expect_true(all(plain$stop_reason == "no_admissible_dose"))
  expect_true(all(safe$stop_reason == "no_admissible_dose"))
  expect_true(safe$settings$extrasafe)
  expect_equal(safe$settings$offset, 0.10)
})

test_that("n_earlystop ends the trial at the current combination", {
  args <- list(target = 0.35, target_resp = 0.25, p_true = matrix(0, 2, 2),
               p_resp = matrix(1, 2, 2), n_cohort = 12, cohort_size = 3,
               n_trials = 3, seed = 2)
  full <- do.call(comb_boin12_simulate, args)
  early <- do.call(comb_boin12_simulate, c(args, list(n_earlystop = 12)))
  simple <- do.call(comb_boin12_simulate,
                    c(args, list(n_earlystop = 3,
                                 n_earlystop_rule = "simple")))
  expect_true(all(full$n_pts[, 1, 1] == 36L))
  expect_true(all(early$n_pts[, 1, 1] == 12L))
  expect_true(all(simple$n_pts[, 1, 1] == 3L))
  expect_true(all(early$stop_reason == "n_earlystop"))
  expect_true(all(simple$stop_reason == "n_earlystop"))
})

test_that("the rules agree with a separate implementation", {
  fx <- utils::read.csv(test_path("fixtures", "comb_boin12_rules.csv"),
                        stringsAsFactors = FALSE)
  groups <- split(fx, paste(fx$tie_break, fx$no_candidate, fx$extrasafe))
  expect_length(groups, 8L)
  counts <- function(f, prefix) {
    array(as.integer(unlist(f[, paste0(prefix, 1:9)])), c(nrow(f), 3L, 3L))
  }
  n_pts <- list()
  for (g in groups) {
    tr <- comb_boin12_simulate(target = 0.30, target_resp = 0.20,
                               p_true = scenario_g$p_true,
                               p_resp = scenario_g$p_resp, n_cohort = 8,
                               cohort_size = 3, n_trials = 30,
                               psi00 = lu_psi00, psi11 = lu_psi11,
                               benchmark_weight = 0.6,
                               extrasafe = as.logical(g$extrasafe[1]),
                               offset = 0.10, tie_break = g$tie_break[1],
                               no_candidate = g$no_candidate[1],
                               seed = 20261117)
    expect_identical(unname(tr$n_pts), counts(g, "n"))
    expect_identical(unname(tr$n_tox), counts(g, "y"))
    expect_identical(unname(tr$n_resp), counts(g, "x"))
    expect_identical(unname(tr$n_tox_resp), counts(g, "z"))
    expect_identical(tr$cohorts_used, g$cohorts_used)
    expect_identical(tr$stop_reason, g$stop_reason)
    n_pts[[g$tie_break[1]]] <- tr$n_pts
  }
  expect_false(identical(n_pts$higher, n_pts$random))
})

test_that("independent outcomes are the same under both links", {
  args <- list(target = 0.33, target_resp = 0.2, p_true = sc16$p_true,
               p_resp = sc16$p_resp, n_cohort = 14, cohort_size = 3,
               n_trials = 200, resp_cor = 0, seed = 4)
  a <- do.call(comb_boin12_simulate, args)
  b <- do.call(comb_boin12_simulate, c(args, list(resp_cor_type = "binary")))
  expect_identical(a$n_pts, b$n_pts)
  expect_identical(a$n_resp, b$n_resp)
  expect_identical(a$n_tox_resp, b$n_tox_resp)
})

test_that("the seed reproduces the trials and the session state is kept", {
  args <- list(target = 0.33, target_resp = 0.2, p_true = sc16$p_true,
               p_resp = sc16$p_resp, n_cohort = 14, cohort_size = 3,
               n_trials = 100, resp_cor = 0.3, seed = 5)
  set.seed(99)
  before <- .Random.seed
  a <- do.call(comb_boin12_simulate, args)
  expect_identical(.Random.seed, before)
  b <- do.call(comb_boin12_simulate, args)
  expect_identical(a$n_pts, b$n_pts)
  expect_identical(a$n_tox_resp, b$n_tox_resp)
})

test_that("invalid arguments are rejected", {
  args <- list(target = 0.33, target_resp = 0.2, p_true = sc16$p_true,
               p_resp = sc16$p_resp, n_cohort = 14, cohort_size = 3,
               n_trials = 10)
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(p_resp = c(0.1, 0.2)))),
               "p_resp")
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(p_resp = matrix(0.2, 3, 2)))),
               "p_resp")
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(p_true = c(0.1, 0.2)))),
               "p_true")
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(start_dose = c(3, 1)))),
               "start_dose")
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(n_star = 0))), "n_star")
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(offset = 0.95))), "offset")
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(extrasafe = NA))), "extrasafe")
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(resp_cor = 1))), "resp_cor")
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(n_earlystop = 0))),
               "n_earlystop")
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(tie_break = "lower"))),
               "should be one of")
  expect_error(do.call(comb_boin12_simulate,
                       modifyList(args, list(no_candidate = "wait"))),
               "should be one of")
})
