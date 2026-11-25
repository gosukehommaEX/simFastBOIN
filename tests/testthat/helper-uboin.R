# Scenarios and published values of Zhou, Lee and Yuan (2019) and trial
# level expectations for U-BOIN. Written by
# 02_Rpackage/Validation/uboin/make_test_expect.py of the development files:
# the scenarios and the table are read from the text of the article and of
# its supporting information, and the trial level expectations come from an
# implementation in Python on a replica of R's Mersenne-Twister stream,
# which agrees with the C++ core trial by trial.

# Scenarios of Table 3 (Simulation A).
zhou2019_scenario <- function(scenario) {
  switch(as.character(scenario),
    "1" = list(p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
               p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65)),
    "4" = list(p_true = c(0.15, 0.25, 0.40, 0.45, 0.50),
               p_resp = c(0.15, 0.45, 0.30, 0.25, 0.20)),
    "8" = list(p_true = c(0.22, 0.45, 0.55, 0.65, 0.70),
               p_resp = c(0.03, 0.10, 0.20, 0.35, 0.40))
  )
}

# U-BOIN rows of Table 3: selection percentages, early stopping, patients
# and the true utilities.
zhou2019_table3 <- function(scenario) {
  switch(as.character(scenario),
    "1" = list(sel = c(1.7, 72.9, 22.4, 2.8, 0.0), stop = 0.2,
               n_pts = c(6.2, 29.9, 13.8, 3.5, 0.5),
               utility = c(43, 69, 63, 56, 50)),
    "4" = list(sel = c(11.9, 74.1, 3.7, 0.4, 0.0), stop = 9.9,
               n_pts = c(15.0, 29.5, 5.1, 1.1, 0.2),
               utility = c(36, 52, 36, 32, 27)),
    "8" = list(sel = c(0.8, 5.5, 1.7, 0.0, 0.0), stop = 92.0,
               n_pts = c(14.3, 9.7, 1.6, 0.1, 0.0),
               utility = c(25, 23, 25, 30, 31))
  )
}

# The admissible rows of Table S4 (3 patients, utilities 0, 15, 25, 100):
# responses, DLTs, responses without a DLT and the printed utility.
zhou2019_table_s4 <- function() {
  data.frame(
    n_resp = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 1, 2),
    n_tox = c(0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
    n_resp_no_tox = c(1, 2, 3, 0, 1, 2, 1, 2, 1, 0, 0, 1, 1),
    utility = c(41.2, 62.5, 83.8, 22.5, 43.8, 65.0, 37.5,
                58.8, 46.2, 18.8, 25.0, 33.8, 40.0)
  )
}

# Scenario 1 of Table 3, resp_assoc 0.2, seed 1, 5 trials.
# One row per trial; n_outcome holds the counts dose by dose, four
# outcomes per dose.
uboin_expected_scenario1 <- function() {
  list(
    n_outcome = rbind(
      c(0, 3, 0, 0, 0, 2, 1, 3, 2, 6, 10, 21, 1, 2, 1, 2, 0, 0, 0, 0),
      c(0, 5, 0, 1, 0, 13, 6, 23, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0),
      c(1, 8, 0, 3, 1, 10, 3, 22, 0, 2, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0),
      c(0, 7, 0, 2, 0, 16, 3, 23, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
      c(0, 4, 0, 2, 1, 5, 2, 4, 1, 12, 10, 13, 0, 0, 0, 0, 0, 0, 0, 0)
    ),
    eliminated = rbind(
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE)
    ),
    stage2 = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    cohorts_used = c(18, 18, 18, 18, 18),
    stop_reason = c("max_sample_size", "max_sample_size",
                    "max_sample_size", "max_sample_size",
                    "max_sample_size"),
    obd = c(3, 2, 2, 2, 3)
  )
}

# Scenario 8 of Table 3, resp_assoc 0.2, seed 2, 8 trials.
# One row per trial; n_outcome holds the counts dose by dose, four
# outcomes per dose.
uboin_expected_scenario8 <- function() {
  list(
    n_outcome = rbind(
      c(3, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(5, 12, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(7, 15, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(3, 9, 0, 0, 8, 6, 0, 1, 3, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0),
      c(7, 15, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(4, 13, 0, 1, 4, 3, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(1, 11, 0, 0, 4, 4, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(4, 13, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    ),
    eliminated = rbind(
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, TRUE, TRUE, TRUE, TRUE),
      c(FALSE, TRUE, TRUE, TRUE, TRUE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE)
    ),
    stage2 = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    cohorts_used = c(4, 6, 8, 12, 8, 9, 8, 6),
    stop_reason = c("no_admissible_dose", "no_admissible_dose",
                    "no_admissible_dose", "no_admissible_dose",
                    "no_admissible_dose", "no_admissible_dose",
                    "no_admissible_dose", "no_admissible_dose"),
    obd = c(NA, NA, NA, NA, NA, NA, NA, NA)
  )
}

# Toxic doses (DLT 0.45 to 0.85), independent outcomes, seed 4, 6 trials.
# One row per trial; n_outcome holds the counts dose by dose, four
# outcomes per dose.
uboin_expected_toxic <- function() {
  list(
    n_outcome = rbind(
      c(4, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(2, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(6, 11, 2, 8, 4, 7, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(11, 12, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(6, 9, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(3, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    ),
    eliminated = rbind(
      c(TRUE, TRUE, TRUE, TRUE, TRUE),
      c(TRUE, TRUE, TRUE, TRUE, TRUE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(TRUE, TRUE, TRUE, TRUE, TRUE)
    ),
    stage2 = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
    cohorts_used = c(2, 2, 18, 11, 7, 2),
    stop_reason = c("lowest_dose_eliminated", "lowest_dose_eliminated",
                    "max_sample_size", "no_admissible_dose",
                    "no_admissible_dose", "lowest_dose_eliminated"),
    obd = c(NA, NA, 1, NA, NA, NA)
  )
}

# Scenario 1 of Table 3, resp_assoc 0.2, cohorts of 2, 15 cohorts,
# n_stage1 6, n_stop 12, start dose 2, seed 3, 5 trials.
# One row per trial; n_outcome holds the counts dose by dose, four
# outcomes per dose.
uboin_expected_options <- function() {
  list(
    n_outcome = rbind(
      c(0, 0, 0, 0, 0, 3, 0, 9, 0, 3, 1, 2, 0, 0, 2, 0, 0, 0, 0, 0),
      c(0, 2, 0, 0, 1, 5, 0, 6, 0, 3, 1, 2, 0, 0, 2, 0, 0, 0, 0, 0),
      c(0, 0, 0, 0, 1, 5, 0, 6, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(0, 0, 0, 0, 0, 3, 0, 9, 0, 3, 0, 1, 2, 1, 2, 3, 0, 0, 0, 0),
      c(0, 0, 0, 0, 1, 2, 1, 8, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    ),
    eliminated = rbind(
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE, FALSE),
      c(FALSE, FALSE, TRUE, TRUE, TRUE)
    ),
    stage2 = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    cohorts_used = c(10, 11, 7, 12, 8),
    stop_reason = c("n_stop", "n_stop",
                    "n_stop", "n_stop",
                    "n_stop"),
    obd = c(2, 2, 2, 2, 2)
  )
}
