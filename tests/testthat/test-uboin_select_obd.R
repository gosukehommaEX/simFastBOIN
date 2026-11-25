# Table S3 of the supporting information of the article: the counts of the
# outcomes (no response, DLT), (no response, no DLT), (response, DLT) and
# (response, no DLT) at five doses, with the utilities 0, 15, 25, 100.
uboin_table_s3 <- function() {
  rbind(c(0, 3, 0, 0),
        c(0, 2, 0, 1),
        c(0, 4, 1, 1),
        c(1, 3, 1, 4),
        c(1, 2, 1, 2))
}

test_that("uboin_select_obd reproduces the example of Table S3", {
  # The utilities printed in Table S3 are 41.2, 31.4, 50.5 and 41.4 for the
  # doses 2 to 5; dose 1 is inadmissible there, which is the efficacy prior
  # Beta(0.05, 0.05), and the next cohort goes to dose 4.
  sel <- uboin_select_obd(uboin_table_s3(), tox_upper = 0.30,
                          resp_lower = 0.20, utility = c(0, 15, 25, 100),
                          prior_resp = c(0.05, 0.05))
  expect_equal(sel$obd, 4L)
  expect_equal(sel$reason, "selected")
  utilities <- unname(unlist(sel[1L, paste0("utility_DL", 2:5)]))
  expect_true(all(abs(utilities - c(41.2, 31.4, 50.5, 41.4)) <= 0.05 + 1e-9))
  expect_equal(utilities[1L], 41.25)
  admissible <- unlist(sel[1L, paste0("admissible_DL", 1:5)])
  expect_equal(unname(admissible), c(FALSE, TRUE, TRUE, TRUE, TRUE))

  # With the default prior dose 1 (0 responses of 3) is admissible.
  sel <- uboin_select_obd(uboin_table_s3(), tox_upper = 0.30,
                          resp_lower = 0.20, utility = c(0, 15, 25, 100))
  expect_true(sel$admissible_DL1)
  expect_equal(sel$utility_DL1, (15 * 3.25 + 25 * 0.25 + 100 * 0.25) / 4)
  expect_equal(sel$obd, 4L)
})

test_that("eliminated and inadmissible doses are not selected", {
  sel <- uboin_select_obd(uboin_table_s3(), tox_upper = 0.30,
                          resp_lower = 0.20, utility = c(0, 15, 25, 100),
                          eliminated = c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(sel$obd, 2L)
  expect_false(sel$admissible_DL4)

  toxic <- rbind(c(3, 0, 0, 0), c(0, 0, 0, 0))
  sel <- uboin_select_obd(toxic, tox_upper = 0.30, resp_lower = 0.20)
  expect_true(is.na(sel$obd))
  expect_equal(sel$reason, "no_admissible_dose")
  # An untried dose is not admissible.
  expect_false(sel$admissible_DL2)
})

test_that("ties go to the lower dose", {
  counts <- rbind(c(0, 2, 0, 1), c(0, 2, 0, 1), c(0, 3, 0, 0))
  sel <- uboin_select_obd(counts, tox_upper = 0.30, resp_lower = 0.20)
  expect_equal(sel$utility_DL1, sel$utility_DL2)
  expect_equal(sel$obd, 1L)
})

test_that("several trials are selected as each trial alone", {
  trials <- uboin_simulate(tox_upper = 0.30, resp_lower = 0.20,
                           p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
                           p_resp = c(0.08, 0.46, 0.25, 0.20, 0.10),
                           n_cohort = 18, cohort_size = 3, n_trials = 40,
                           seed = 6)
  all_trials <- uboin_select_obd(trials$n_outcome, tox_upper = 0.30,
                                 resp_lower = 0.20,
                                 eliminated = trials$eliminated)
  expect_equal(nrow(all_trials), 40L)
  for (t in c(1L, 17L, 40L)) {
    one <- uboin_select_obd(trials$n_outcome[t, , ], tox_upper = 0.30,
                            resp_lower = 0.20,
                            eliminated = trials$eliminated[t, ])
    expect_equal(one$obd, all_trials$obd[t])
    expect_equal(unlist(one[1L, -1L]), unlist(all_trials[t, -1L]))
  }
})

test_that("uboin_select_obd validates its arguments", {
  expect_error(uboin_select_obd(matrix(1.5, 2, 4), tox_upper = 0.3,
                                resp_lower = 0.2),
               "'n_outcome' must hold non-negative whole numbers")
  expect_error(uboin_select_obd(matrix(1, 2, 3), tox_upper = 0.3,
                                resp_lower = 0.2),
               "'n_outcome' must have four outcomes")
  expect_error(uboin_select_obd(matrix(1, 2, 4), tox_upper = 0.3,
                                resp_lower = 0.2, eliminated = c(1, 0)),
               "'eliminated' must be logical")
  expect_error(uboin_select_obd(matrix(1, 2, 4), tox_upper = 0.3,
                                resp_lower = 0.2, eliminated = c(TRUE)),
               "'eliminated' must have one row per trial")
})
