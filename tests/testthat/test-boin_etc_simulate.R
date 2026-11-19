# Scenarios 3 and 7 of Kakizume et al. (2024) and three more (rows are the
# levels of drug A, from A1).
etc_scenarios <- list(
  Sc3 = list(
    p_true = matrix(c(0.02, 0.06, 0.12, 0.18, 0.04, 0.11, 0.15, 0.30,
                      0.05, 0.13, 0.30, 0.50, 0.10, 0.30, 0.48, 0.54),
                    4, byrow = TRUE),
    p_resp = matrix(c(0.05, 0.10, 0.20, 0.30, 0.10, 0.20, 0.30, 0.60,
                      0.20, 0.30, 0.60, 0.65, 0.30, 0.60, 0.65, 0.70),
                    4, byrow = TRUE)),
  Sc7 = list(
    p_true = matrix(c(0.01, 0.04, 0.11, 0.15, 0.30, 0.03, 0.05, 0.13, 0.30,
                      0.50, 0.07, 0.10, 0.30, 0.48, 0.54), 3, byrow = TRUE),
    p_resp = matrix(c(0.05, 0.10, 0.20, 0.30, 0.60, 0.10, 0.20, 0.30, 0.60,
                      0.65, 0.20, 0.30, 0.60, 0.65, 0.70), 3, byrow = TRUE)),
  toxic = list(
    p_true = matrix(c(0.30, 0.50, 0.60, 0.45, 0.65, 0.80), 2, byrow = TRUE),
    p_resp = matrix(c(0.20, 0.40, 0.50, 0.40, 0.50, 0.60), 2, byrow = TRUE)),
  "4x3" = list(
    p_true = matrix(c(0.02, 0.05, 0.10, 0.05, 0.10, 0.20, 0.10, 0.20, 0.35,
                      0.20, 0.35, 0.50), 4, byrow = TRUE),
    p_resp = matrix(c(0.10, 0.20, 0.40, 0.20, 0.40, 0.60, 0.40, 0.60, 0.60,
                      0.60, 0.60, 0.60), 4, byrow = TRUE)),
  flat = list(
    p_true = matrix(c(0.18, 0.20, 0.22, 0.24, 0.20, 0.22, 0.24, 0.26, 0.22,
                      0.24, 0.26, 0.28), 3, byrow = TRUE),
    p_resp = matrix(0.33, 3, 4))
)

# Settings of the fixture boin_etc_rules.csv, computed with the Python
# implementation of the rules fed with R's random numbers and checked against
# the C++ core (see Validation/boin_etc/make_etc_fixture.py of the
# development folder).
etc_fixture_settings <- data.frame(
  setting = paste0("F", 1:7),
  scenario = c("Sc3", "Sc3", "Sc7", "Sc7", "toxic", "4x3", "flat"),
  subtrials = c("etc1", "etc1", "etc2", "etc3", "etc1", "etc3", "etc2"),
  etc1_start = c("diagonal", "below", "diagonal", "diagonal", "diagonal",
                 "diagonal", "diagonal"),
  rule4 = c("paper", "authors_code", "paper", "paper", "paper", "paper",
            "paper"),
  n_earlystop = c(6, 6, 9, 6, 6, 6, 6),
  n_cohort = c(16, 16, 24, 16, 8, 12, 16),
  seed = c(20261120, 20261120, 7, 7, 11, 3, 5),
  n_trials = c(40, 40, 30, 40, 60, 40, 40),
  stringsAsFactors = FALSE
)

etc_fixture <- function() {
  utils::read.csv(test_path("fixtures", "boin_etc_rules.csv"),
                  stringsAsFactors = FALSE, colClasses = "character")
}

etc_run_setting <- function(s) {
  sc <- etc_scenarios[[s$scenario]]
  boin_etc_simulate(target = 0.30, target_resp = 0.60, p_true = sc$p_true,
                    p_resp = sc$p_resp, n_cohort = s$n_cohort,
                    cohort_size = 3, n_trials = s$n_trials,
                    subtrials = s$subtrials, n_earlystop = s$n_earlystop,
                    etc1_start = s$etc1_start, rule4 = s$rule4,
                    seed = s$seed)
}

etc_split <- function(x) as.integer(strsplit(x, ";", fixed = TRUE)[[1L]])

test_that("boin_etc_simulate returns a well formed object", {
  sc <- etc_scenarios$Sc3
  trials <- boin_etc_simulate(target = 0.30, target_resp = 0.60,
                              p_true = sc$p_true, p_resp = sc$p_resp,
                              n_cohort = 16, cohort_size = 3, n_trials = 100,
                              seed = 1)
  expect_s3_class(trials, "boin_etc_trials")
  for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
    expect_equal(dim(trials[[nm]]), c(100L, 4L, 4L))
    expect_type(trials[[nm]], "integer")
    expect_equal(dimnames(trials[[nm]])[[2L]], paste0("A", 1:4))
    expect_equal(dimnames(trials[[nm]])[[3L]], paste0("B", 1:4))
  }
  expect_type(trials$excluded, "logical")
  expect_equal(dim(trials$excluded), c(100L, 4L, 4L))
  expect_true(all(trials$n_tox <= trials$n_pts))
  expect_true(all(trials$n_tox_resp <= pmin(trials$n_tox, trials$n_resp)))
  expect_true(all(apply(trials$n_pts, 1L, sum) == 3L * trials$cohorts_used))
  expect_true(all(trials$cohorts_used <= 16L))
  expect_true(all(trials$n_subtrials >= 1L & trials$n_subtrials <= 4L))
  expect_true(all(trials$stop_reason %in% boin_etc_stop_reasons()))
  expect_identical(trials$early_stop, trials$stop_reason == "lowest_excluded")
  expect_equal(trials$design, "BOIN-ETC1")
  bd <- boin_et_boundary(0.30, 0.60)
  expect_equal(trials$boundary$lambda1, bd$lambda1)
  expect_equal(trials$boundary$lambda2, bd$lambda2)
  expect_equal(trials$boundary$eta1, bd$eta1)
  expect_equal(c(bd$lambda1, bd$lambda2, bd$eta1), c(0.14, 0.35, 0.48))
  expect_equal(trials$settings$max_total_pts, 48)
  expect_equal(trials$settings$p_eli, 0.30)
})

test_that("boin_etc_simulate reproduces the trials of the Python rules", {
  fx <- etc_fixture()
  for (i in seq_len(nrow(etc_fixture_settings))) {
    s <- etc_fixture_settings[i, ]
    trials <- etc_run_setting(s)
    rows <- fx[fx$setting == s$setting, , drop = FALSE]
    expect_equal(nrow(rows), s$n_trials)
    n_cells <- prod(dim(trials$n_pts)[2:3])
    for (t in seq_len(nrow(rows))) {
      for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
        expect_identical(as.vector(trials[[nm]][t, , ]),
                         etc_split(rows[[nm]][t]))
      }
      expect_identical(as.integer(as.vector(trials$excluded[t, , ])),
                       etc_split(rows$excluded[t]))
      expect_equal(length(etc_split(rows$n_pts[t])), n_cells)
    }
    expect_identical(trials$cohorts_used, as.integer(rows$cohorts_used))
    expect_identical(trials$n_subtrials, as.integer(rows$n_subtrials))
    expect_identical(trials$stop_reason, rows$stop_reason)
  }
})

test_that("the subtrials start where the rules say", {
  # No DLT and no response except at A3B1: the first subtrial escalates to
  # A3B1 and ends there (n_stop = 3), so the rows start at A2 (expected
  # counts from the Python rules, make_etc_fixture.py of the development
  # folder). etc1 starts row A1 one level of drug B to the left of the
  # candidate A2B4 (diagonal) or below it; etc2 starts every row at B2; etc3
  # runs the L-shapes from their lowest combination.
  p_resp <- matrix(0, 3, 4)
  p_resp[3, 1] <- 1
  args <- list(target = 0.30, target_resp = 0.60, p_true = matrix(0, 3, 4),
               p_resp = p_resp, n_cohort = 12, cohort_size = 3,
               n_earlystop = 3, n_trials = 2, seed = 1)
  run <- function(...) do.call(boin_etc_simulate, c(args, list(...)))
  diag <- run(subtrials = "etc1")
  below <- run(subtrials = "etc1", etc1_start = "below")
  etc2 <- run(subtrials = "etc2")
  etc3 <- run(subtrials = "etc3")
  full <- c(3L, 3L, 3L, 3L)
  # rbind() would name the rows after the variable full.
  rows <- function(...) rbind(..., deparse.level = 0)
  expect_equal(unname(diag$n_pts[1, , ]),
               rows(c(3L, 0L, 3L, 3L), full, c(3L, 0L, 0L, 0L)))
  expect_equal(unname(below$n_pts[1, , ]),
               rows(c(3L, 0L, 0L, 3L), full, c(3L, 0L, 0L, 0L)))
  expect_equal(unname(etc2$n_pts[1, , ]),
               rows(full, full, c(3L, 0L, 0L, 0L)))
  expect_equal(unname(etc3$n_pts[1, , ]),
               rows(full, full, c(3L, 0L, 0L, 0L)))
  expect_equal(diag$cohorts_used, c(8L, 8L))
  expect_equal(below$cohorts_used, c(7L, 7L))
  expect_equal(etc3$n_subtrials, c(3L, 3L))
  # Every response certain: the trial stays where every subtrial starts.
  # After the candidate A1B1 of the first subtrial only the row A1 is run.
  args$p_resp <- matrix(1, 3, 4)
  etc2 <- run(subtrials = "etc2")
  etc3 <- run(subtrials = "etc3")
  expect_equal(unname(etc2$n_pts[1, , ]), rbind(c(3L, 3L, 0L, 0L), 0L, 0L))
  expect_equal(unname(etc3$n_pts[1, , ]), rbind(c(3L, 3L, 3L, 0L), 0L, 0L))
  expect_equal(etc2$n_subtrials, c(2L, 2L))
  expect_equal(etc3$n_subtrials, c(3L, 3L))
  expect_true(all(etc2$stop_reason == "subtrials_completed"))
})

test_that("the lowest combination excluded ends the trial", {
  trials <- boin_etc_simulate(target = 0.30, target_resp = 0.60,
                              p_true = matrix(1, 2, 3),
                              p_resp = matrix(0.5, 2, 3), n_cohort = 8,
                              cohort_size = 3, n_trials = 3, seed = 1)
  expect_true(all(trials$stop_reason == "lowest_excluded"))
  expect_true(all(trials$early_stop))
  expect_equal(trials$cohorts_used, rep(1L, 3))
  expect_true(all(trials$excluded[, 2, 3]))
  expect_true(all(trials$excluded[, , 1]))
})

test_that("boin_etc_simulate checks its arguments", {
  sc <- etc_scenarios$Sc7
  args <- list(target = 0.30, target_resp = 0.60, p_true = sc$p_true,
               p_resp = sc$p_resp, n_cohort = 16, cohort_size = 3,
               n_trials = 10)
  expect_error(do.call(boin_etc_simulate,
                       c(args, list(subtrials = "etc4"))))
  expect_error(do.call(boin_etc_simulate, c(args, list(rule4 = "code"))))
  expect_error(do.call(boin_etc_simulate, c(args, list(n_earlystop = 0))))
  expect_error(do.call(boin_etc_simulate, c(args, list(cutoff_eli = 2))))
  args$p_true <- matrix(c(0.1, 0.2, 0.3), 1)
  args$p_resp <- matrix(c(0.3, 0.4, 0.5), 1)
  expect_error(do.call(boin_etc_simulate, args),
               "at least two levels of each drug")
})

test_that("the binary link and the seed are handled", {
  sc <- etc_scenarios$Sc7
  args <- list(target = 0.30, target_resp = 0.60, p_true = sc$p_true,
               p_resp = sc$p_resp, n_cohort = 16, cohort_size = 3,
               n_trials = 50, subtrials = "etc3")
  a <- do.call(boin_etc_simulate,
               c(args, list(resp_cor = 0.3, resp_cor_type = "binary",
                            seed = 9)))
  b <- do.call(boin_etc_simulate,
               c(args, list(resp_cor = 0.3, resp_cor_type = "binary",
                            seed = 9)))
  expect_identical(a$n_pts, b$n_pts)
  expect_equal(a$settings$resp_cor_type, "binary")
  set.seed(42)
  before <- stats::runif(1)
  set.seed(42)
  do.call(boin_etc_simulate, c(args, list(seed = 1)))
  expect_identical(stats::runif(1), before)
})
