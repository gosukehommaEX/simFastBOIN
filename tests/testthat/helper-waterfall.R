# Scenarios and published cells of Zhang and Yuan (2016), and trial level
# expectations for the waterfall design. Written by
# 02_Rpackage/Validation/waterfall/make_test_expect.py of the development
# files: the scenarios and the published cells are read from the text of
# the article and of its supporting information, and the trial level
# expectations come from an implementation in Python on a replica of R's
# Mersenne-Twister stream, which agrees with the C++ core trial by trial.

# True DLT probabilities of a scenario of Table III.
zhang2016_p_true <- function(scenario) {
  switch(as.character(scenario),
    "1" = matrix(c(0.03, 0.10, 0.28,
                   0.10, 0.30, 0.50),
                 nrow = 2, byrow = TRUE),
    "4" = matrix(c(0.30, 0.40, 0.50,
                   0.42, 0.49, 0.55),
                 nrow = 2, byrow = TRUE),
    "12" = matrix(c(0.01, 0.03, 0.30, 0.45, 0.52,
                    0.30, 0.41, 0.52, 0.61, 0.73,
                    0.49, 0.51, 0.57, 0.64, 0.77),
                  nrow = 3, byrow = TRUE)
  )
}

# Waterfall cells of Table IV (selection percentages), Table S1 (average
# numbers of patients) and Table V (percentage of trials selecting the
# contour, percentages of patients above and on the contour).
zhang2016_published <- function(scenario) {
  switch(as.character(scenario),
    "1" = list(
      sel = matrix(c(0.3, 15.5, 84.2,
                     21.1, 59.8, 18.5),
                   nrow = 2, byrow = TRUE),
      n_pts = matrix(c(3.4, 1.7, 7.2,
                       5.4, 6.6, 2.5),
                     nrow = 2, byrow = TRUE),
      pcs = 50.4, above = 9.4, at = 51.3
    ),
    "4" = list(
      sel = matrix(c(56.5, 24.4, 6.9,
                     23.4, 5.9, 0.3),
                   nrow = 2, byrow = TRUE),
      n_pts = matrix(c(8.3, 2.5, 0.8,
                       4.4, 1.0, 0.1),
                     nrow = 2, byrow = TRUE),
      pcs = 48.5, above = 51.7, at = 48.3
    ),
    "12" = list(
      sel = matrix(c(4.7, 14.5, 54.7, 21.8, 4.3,
                     64.8, 20.1, 2.1, 0.2, 0.0,
                     14.7, 2.1, 0.2, 0.0, 0.0),
                   nrow = 3, byrow = TRUE),
      n_pts = matrix(c(5.8, 4.4, 7.7, 3.7, 0.7,
                       10.2, 3.0, 0.9, 0.1, 0.0,
                       4.6, 0.9, 0.2, 0.0, 0.0),
                     nrow = 3, byrow = TRUE),
      pcs = 35.9, above = 33.2, at = 42.6
    )
  )
}

# A matrix with steep DLT probabilities, in which subtrials often end
# without a candidate MTD.
waterfall_steep_p_true <- function() {
  matrix(c(0.05, 0.25, 0.45, 0.60, 0.70,
           0.20, 0.40, 0.55, 0.70, 0.80,
           0.30, 0.50, 0.65, 0.80, 0.90),
         nrow = 3, byrow = TRUE)
}

# Scenario 4, cohorts of 3, seed 3, boin_compatible = TRUE, 8 trials.
# One row per trial, the combinations in column-major order.
waterfall_expected_compatible <- function() {
  list(
    n_pts = rbind(
      c(9, 9, 0, 0, 0, 0),
      c(12, 3, 0, 0, 0, 0),
      c(0, 0, 0, 0, 0, 0),
      c(12, 3, 0, 0, 0, 0),
      c(6, 9, 9, 3, 0, 0),
      c(12, 0, 0, 0, 0, 0),
      c(12, 0, 0, 0, 0, 0),
      c(3, 9, 0, 6, 0, 0)
    ),
    n_tox = rbind(
      c(1, 3, 0, 0, 0, 0),
      c(3, 3, 0, 0, 0, 0),
      c(0, 0, 0, 0, 0, 0),
      c(3, 3, 0, 0, 0, 0),
      c(1, 1, 3, 2, 0, 0),
      c(4, 0, 0, 0, 0, 0),
      c(6, 0, 0, 0, 0, 0),
      c(0, 1, 0, 4, 0, 0)
    ),
    eliminated = rbind(
      c(0, 0, 0, 1, 0, 1),
      c(0, 1, 0, 1, 0, 1),
      c(0, 0, 0, 0, 0, 0),
      c(0, 1, 0, 1, 0, 1),
      c(0, 0, 0, 1, 0, 1),
      c(0, 1, 0, 1, 0, 1),
      c(0, 1, 0, 1, 0, 1),
      c(0, 0, 0, 1, 0, 1)
    ),
    n_subtrials = c(2, 1, 1, 1, 2, 1, 1, 2),
    n_enrolled = c(24, 15, 9, 15, 27, 12, 12, 24),
    stop_reason = c("later_subtrial_no_mtd", "completed",
                    "first_subtrial_no_mtd", "completed",
                    "completed", "completed",
                    "completed", "later_subtrial_no_mtd")
  )
}

# Scenario 4, cohorts of 3, seed 3, default rules, 8 trials.
# One row per trial, the combinations in column-major order.
waterfall_expected_default <- function() {
  list(
    n_pts = rbind(
      c(9, 9, 6, 0, 0, 0),
      c(12, 3, 0, 0, 0, 0),
      c(9, 0, 0, 0, 0, 0),
      c(12, 3, 0, 0, 0, 0),
      c(6, 9, 9, 3, 0, 0),
      c(12, 0, 0, 0, 0, 0),
      c(12, 0, 0, 0, 0, 0),
      c(3, 9, 6, 6, 0, 0)
    ),
    n_tox = rbind(
      c(1, 3, 5, 0, 0, 0),
      c(3, 3, 0, 0, 0, 0),
      c(6, 0, 0, 0, 0, 0),
      c(3, 3, 0, 0, 0, 0),
      c(1, 1, 3, 2, 0, 0),
      c(4, 0, 0, 0, 0, 0),
      c(6, 0, 0, 0, 0, 0),
      c(0, 1, 4, 4, 0, 0)
    ),
    eliminated = rbind(
      c(0, 0, 1, 1, 1, 1),
      c(0, 1, 0, 1, 0, 1),
      c(1, 1, 0, 1, 0, 1),
      c(0, 1, 0, 1, 0, 1),
      c(0, 0, 0, 1, 0, 1),
      c(0, 1, 0, 1, 0, 1),
      c(0, 1, 0, 1, 0, 1),
      c(0, 0, 1, 1, 1, 1)
    ),
    n_subtrials = c(2, 1, 1, 1, 2, 1, 1, 2),
    n_enrolled = c(24, 15, 9, 15, 27, 12, 12, 24),
    stop_reason = c("later_subtrial_no_mtd", "completed",
                    "first_subtrial_no_mtd", "completed",
                    "completed", "completed",
                    "completed", "later_subtrial_no_mtd")
  )
}

# The steep matrix, cohorts of 2, all options, seed 11, 6 trials.
# One row per trial, the combinations in column-major order.
waterfall_expected_options <- function() {
  list(
    n_pts = rbind(
      c(5, 12, 0, 12, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0),
      c(5, 12, 8, 12, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0),
      c(1, 1, 13, 0, 12, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(1, 5, 13, 2, 12, 2, 10, 2, 0, 4, 0, 0, 0, 0, 0),
      c(3, 9, 12, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(1, 1, 1, 12, 12, 13, 4, 4, 12, 0, 0, 0, 0, 0, 0)
    ),
    n_tox = rbind(
      c(0, 3, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
      c(0, 3, 3, 3, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0),
      c(0, 0, 3, 0, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(0, 0, 4, 1, 4, 2, 4, 1, 0, 4, 0, 0, 0, 0, 0),
      c(0, 2, 3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      c(0, 0, 0, 4, 1, 5, 3, 3, 5, 0, 0, 0, 0, 0, 0)
    ),
    eliminated = rbind(
      c(0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1),
      c(0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1),
      c(0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1),
      c(0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
      c(0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1),
      c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    ),
    n_subtrials = c(2, 2, 2, 3, 2, 3),
    n_enrolled = c(31, 41, 31, 51, 28, 60),
    stop_reason = c("completed", "completed",
                    "later_subtrial_no_mtd", "later_subtrial_no_mtd",
                    "later_subtrial_no_mtd", "completed")
  )
}

waterfall_expected_seed <- function() 3L

