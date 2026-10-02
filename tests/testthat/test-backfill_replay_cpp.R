# backfill_replay_cpp() runs the engine of bf_boin_simulate() and
# be_boin_simulate() on scripted random numbers, so that the trial examples of
# the articles can be replayed patient by patient. Every dose has a DLT and a
# response probability of 0.5 with uniform times, so that a variate u below
# 0.5 gives a DLT, or a response, at day window * u / 0.5. Arrivals use
# uniform times with an accrual rate of 0.02, so that a variate g gives a gap
# of 100 * g days.

replay_example <- function(patients, design, window, n_cohort, next_arrival,
                           stay_on_1_of_3 = FALSE, max_pending_ratio = 1,
                           min_follow_up = 0, conflict_dose = 0L,
                           n_earlystop = 18L) {
  n_cap <- 12L
  bound <- boin_boundary(0.25, 3L * n_cohort + n_cap * 5L,
                         stay_on_1_of_3 = stay_on_1_of_3)
  b_elim <- bound$b_elim
  b_elim[is.na(b_elim)] <- 0L

  dlt_u <- ifelse(patients$dlt_day > 0, 0.5 * patients$dlt_day / window, 0.99)
  resp_u <- ifelse(patients$response_day > 0,
                   0.5 * patients$response_day / window, 0.99)
  de_u <- dlt_u[!patients$backfill]
  extra_u <- unlist(lapply(seq_len(nrow(patients)), function(k) {
    if (patients$backfill[k]) c(dlt_u[k], resp_u[k]) else resp_u[k]
  }))
  # The gap after a patient runs from the day the patient was enrolled.
  arrivals <- c(patients$arrival[-1L], next_arrival)
  gap_u <- (arrivals - patients$enrolled[seq_along(arrivals)]) / 100

  backfill_replay_cpp(
    de_u = de_u, gap_u = gap_u, extra_u = extra_u,
    estimate = if (design == "be_boin") 1L else 0L,
    p_true = rep(0.5, 5), p_resp = rep(0.5, 5),
    cohort_size = rep(3L, n_cohort), start_dose = 1L,
    n_earlystop = as.integer(n_earlystop), early_stop_simple = FALSE,
    extrasafe = FALSE, target = 0.25, cutoff_eli = 0.95, offset = 0.05,
    b_esc = bound$b_esc, b_deesc = bound$b_deesc, b_elim = as.integer(b_elim),
    max_total_pts = 3L * n_cohort, lambda_e = bound$lambda_e,
    lambda_d = bound$lambda_d, max_pending_ratio = max_pending_ratio,
    min_completed = 0L, min_follow_up = min_follow_up, window = window,
    accrual = 1L, accrual_rate = 0.02, dlt_time = 1L, late_fraction = 0.5,
    resp_window = window, resp_late_fraction = 0.5, resp_cor = 0,
    n_cap = n_cap, backfill_dose = 0L, conflict_dose = conflict_dose,
    no_slot = 0L
  )
}

selected_mtd <- function(res) {
  n_pts <- matrix(tabulate(res$dose, nbins = 5L), nrow = 1L)
  n_tox <- matrix(tabulate(res$dose[res$dlt], nbins = 5L), nrow = 1L)
  boin_select_mtd(n_pts = n_pts, n_tox = n_tox, target = 0.25)$mtd
}

test_that("the trial example of Zhao et al. (2024) is replayed", {
  # Figure 1 and Supplementary Section B: target 0.25, n_stop 9, n_cap 12,
  # cohorts of three, one DLT out of three stays. The figure gives the order
  # of the patients, not their days; the days below, with a window of 30, are
  # chosen to match the narrative (for example, patients 10 and 11 have
  # completed the assessment when patient 14 arrives, patient 21 has not when
  # patient 24 arrives, and patient 24 has not when patient 25 arrives).
  patients <- data.frame(
    arrival = c(0, 1, 2, 33, 34, 35, 66, 67, 68, 70, 72, 80, 90, 103, 104,
                105, 110, 115, 136, 137, 138, 145, 150, 155, 170, 171, 172,
                180, 185),
    dlt_day = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0,
                0, 0, 1, 2, 3, 0, 1, 1, 0, 0, 4, 0, 0, 5),
    response_day = c(0, 0, 0, 0, 5, 0, 0, 0, 5, 5, 0, 0, 0, 0, 5,
                     5, 0, 5, 5, 5, 0, 5, 0, 0, 0, 0, 5, 0, 0),
    backfill = c(rep(FALSE, 9), rep(TRUE, 4), rep(FALSE, 3), TRUE, TRUE,
                 rep(FALSE, 3), TRUE, TRUE, TRUE, rep(FALSE, 3), TRUE, TRUE)
  )
  patients$enrolled <- patients$arrival
  res <- replay_example(patients, design = "bf_boin", window = 30,
                        n_cohort = 10L, next_arrival = 205,
                        stay_on_1_of_3 = TRUE, n_earlystop = 9L)

  # Doses of the figure: escalation to 4 after pooling doses 2 and 3 for
  # patient 14, dose 4 closed for patient 24, de-escalation from 5 to 3 after
  # pooling doses 4 and 5 for patient 25.
  expect_identical(res$dose, c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 2L, 2L, 2L,
                               2L, 4L, 4L, 4L, 3L, 3L, 5L, 5L, 5L, 4L, 4L, 3L,
                               3L, 3L, 3L, 2L, 2L))
  expect_identical(res$backfill, patients$backfill)
  # The days are sums of scripted gaps, equal up to rounding.
  expect_equal(res$entry, patients$arrival)
  expect_identical(res$dlt, patients$dlt_day > 0)
  # Patient 30 arrives once nine patients at dose 3 have 2 DLTs: the decision
  # is to stay, so n_stop ends the trial after six cohorts.
  expect_identical(res$stop_code, 2L)
  expect_identical(res$cohorts_used, 6L)
  expect_identical(selected_mtd(res), 3L)
  # Every scripted number was used: none was left over or missing.
  expect_identical(res$used, c(18L, 29L, 40L))
})

test_that("the trial example of Chen et al. (2026) is replayed", {
  # Supplementary Section S2: target 0.25, window of 60 days, seven cohorts of
  # three, rules 1 (A = 51) and 2 (B = 25). The days follow the narrative where
  # it states them; where the text and the figure differ, the text is followed
  # (patient 10, not 11, has the DLT), and decision days are moved by a day or
  # two so that no event coincides with an arrival.
  patients <- data.frame(
    arrival = c(0, 10, 35, 45, 85, 100, 150, 170, 180, 190, 200, 210, 241, 250,
                260, 270, 280, 290, 302, 310, 318, 337, 339, 351, 360, 380,
                390, 405, 410, 420),
    dlt_day = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 40, 0, 0, 0, 48, 0, 30, 0, 0, 0, 0,
                47, 0, 0, 47, 0, 0, 0, 0, 40, 0),
    response_day = c(0, 0, 0, 0, 30, 0, 0, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    backfill = c(rep(FALSE, 9), rep(TRUE, 3), rep(FALSE, 3), rep(TRUE, 3),
                 rep(FALSE, 3), TRUE, TRUE, rep(FALSE, 3), TRUE,
                 rep(FALSE, 3))
  )
  # Patient 4 arrives during the suspension and is enrolled on day 70, when
  # patients 1 and 2 have completed the assessment and patient 3 has been
  # followed for more than a quarter of the window.
  patients$enrolled <- replace(patients$arrival, 4L, 70)
  res <- replay_example(patients, design = "be_boin", window = 60,
                        n_cohort = 7L, next_arrival = numeric(0),
                        max_pending_ratio = 0.49, min_follow_up = 0.25,
                        conflict_dose = 1L)

  # Cohorts at doses 1, 2, 3, 4, 3, 4 and 3; backfilling to dose 2 while
  # cohort 3 is pending, to dose 3 while cohort 4 is, to dose 2 while the
  # escalation from dose 3 waits for rules 1 and 2, and to dose 3 again.
  expect_identical(res$dose, c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 2L, 2L, 2L,
                               4L, 4L, 4L, 3L, 3L, 3L, 3L, 3L, 3L, 2L, 2L, 4L,
                               4L, 4L, 3L, 3L, 3L, 3L))
  expect_identical(res$backfill, patients$backfill)
  expect_equal(res$entry, patients$enrolled)
  expect_identical(res$dlt, patients$dlt_day > 0)
  # The seventh cohort ends the trial; 30 patients, 9 of them backfilled.
  expect_identical(res$stop_code, 3L)
  expect_identical(res$cohorts_used, 7L)
  expect_identical(sum(res$backfill), 9L)
  expect_identical(selected_mtd(res), 3L)
  expect_identical(res$used, c(21L, 29L, 39L))
})
