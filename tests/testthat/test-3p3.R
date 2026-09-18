test_that("the selection probabilities of the exact 3+3 add up", {
  scenarios <- list(
    c(0.30, 0.48, 0.67),
    c(0.16, 0.30, 0.44),
    c(0.42, 0.43, 0.44),
    c(0.08, 0.20, 0.30),
    c(0.05, 0.12, 0.20, 0.30, 0.45)
  )

  for (p in scenarios) {
    for (rule in c("previous", "expand")) {
      oc <- oc_3p3(p, mtd_rule = rule)
      expect_equal(sum(oc$sel_percent) + oc$percent_no_mtd, 100)
      expect_true(all(oc$sel_percent >= 0))
      expect_equal(oc$total_n_pts, sum(oc$n_pts_dose))
      expect_equal(oc$total_n_tox, sum(oc$n_tox_dose))
    }
  }
})

test_that("the exact 3+3 matches a hand calculation at a single dose", {
  q <- 0.1
  no_dlt <- (1 - q)^3
  one_dlt <- 3 * q * (1 - q)^2
  escalate <- no_dlt + one_dlt * no_dlt

  # With one dose the trial either escalates out of the top dose, which selects
  # it, or declares it too toxic, which leaves no MTD.
  plain <- oc_3p3(q)
  expect_equal(plain$sel_percent, escalate * 100)
  expect_equal(plain$percent_no_mtd, (1 - escalate) * 100)
  expect_equal(plain$n_pts_dose, 3 + 3 * one_dlt)
  expect_equal(plain$n_tox_dose, q * (3 + 3 * one_dlt))

  # The other rule expands the dose to six patients whenever only three were
  # treated, and then needs at most one DLT out of the six.
  only_three <- no_dlt / escalate
  pass <- no_dlt + one_dlt
  expanded <- oc_3p3(q, mtd_rule = "expand")
  expect_equal(expanded$sel_percent,
               escalate * (only_three * pass + (1 - only_three)) * 100)
  expect_equal(expanded$n_pts_dose,
               3 + 3 * one_dlt + escalate * only_three * 3)
})

test_that("the exact 3+3 is right in the deterministic cases", {
  # A harmless dose followed by a certainly toxic one.
  plain <- oc_3p3(c(0, 1))
  expect_equal(plain$sel_percent, c(100, 0))
  expect_equal(plain$percent_no_mtd, 0)
  expect_equal(plain$n_pts_dose, c(3, 3))
  expect_equal(plain$n_tox_dose, c(0, 3))

  # The same trial under the other rule expands the lower dose to six.
  expanded <- oc_3p3(c(0, 1), mtd_rule = "expand")
  expect_equal(expanded$sel_percent, c(100, 0))
  expect_equal(expanded$n_pts_dose, c(6, 3))
  expect_equal(expanded$total_n_pts, 9)

  # A certainly toxic lowest dose leaves no MTD and stops after one cohort.
  for (rule in c("previous", "expand")) {
    hopeless <- oc_3p3(c(1, 1), mtd_rule = rule)
    expect_equal(hopeless$percent_no_mtd, 100)
    expect_equal(hopeless$sel_percent, c(0, 0))
    expect_equal(hopeless$n_pts_dose, c(3, 0))
    expect_equal(hopeless$n_tox_dose, c(3, 0))
  }
})

test_that("start_dose leaves the doses below it untouched", {
  for (rule in c("previous", "expand")) {
    oc <- oc_3p3(c(0, 0, 1), mtd_rule = rule, start_dose = 2)
    expect_equal(oc$sel_percent, c(0, 100, 0))
    expect_equal(oc$n_pts_dose[1], 0)
    expect_equal(oc$n_pts_dose[3], 3)
    expect_equal(oc$n_pts_dose[2], if (rule == "previous") 3 else 6)
  }
  expect_error(oc_3p3(c(0.1, 0.2), start_dose = 3), "must not exceed")
})

test_that("simulation agrees with the exact 3+3", {
  p <- c(0.16, 0.30, 0.44)

  for (rule in c("previous", "expand")) {
    exact <- oc_3p3(p, mtd_rule = rule)
    simulated <- sim_3p3(p, n_trials = 20000, mtd_rule = rule, seed = 11)

    # Three Monte Carlo standard errors at 20000 trials is about one point for a
    # percentage near a half, so 1.5 points is generous without being vacuous.
    # The differences are compared on their own scale rather than through the
    # relative tolerance of expect_equal(), which is not meaningful for a
    # percentage that can be close to zero.
    worst <- function(a, b) max(abs(a - b))

    expect_lt(worst(simulated$sel_percent, exact$sel_percent), 1.5)
    expect_lt(worst(simulated$percent_no_mtd, exact$percent_no_mtd), 1.5)
    expect_lt(worst(simulated$n_pts_dose, exact$n_pts_dose), 0.2)
    expect_lt(worst(simulated$n_tox_dose, exact$n_tox_dose), 0.2)
    expect_lt(worst(simulated$overdose$pct_patients,
                    exact$overdose$pct_patients), 1.5)
    expect_lt(worst(simulated$overdose$pct_trials_any,
                    exact$overdose$pct_trials_any), 1.5)
    expect_lt(worst(simulated$overdose$pct_trials_mtd_above,
                    exact$overdose$pct_trials_mtd_above), 1.5)
    expect_identical(names(simulated), names(exact))
    expect_identical(names(simulated$overdose), names(exact$overdose))
  }
})

test_that("the overdose component is consistent with the dose level figures", {
  p <- c(0.20, 0.40, 0.60)

  for (rule in c("previous", "expand")) {
    oc <- oc_3p3(p, mtd_rule = rule, overdose_cutoff = 1 / 3)
    expect_equal(oc$overdose$doses, c(2L, 3L))
    expect_equal(oc$overdose$pct_patients,
                 sum(oc$n_pts_dose[2:3]) / oc$total_n_pts * 100)
    expect_equal(oc$overdose$avg_n_patients, sum(oc$n_pts_dose[2:3]))
    expect_equal(oc$overdose$pct_trials_mtd_above, sum(oc$sel_percent[2:3]))
    expect_equal(oc$overdose$pct_mtd_above_when_selected,
                 sum(oc$sel_percent[2:3]) / (100 - oc$percent_no_mtd) * 100)
  }

  # No dose above the cutoff means nothing to report rather than a missing value.
  none <- oc_3p3(c(0.08, 0.20, 0.30), overdose_cutoff = 1 / 3)
  expect_length(none$overdose$doses, 0L)
  expect_equal(none$overdose$pct_patients, 0)
  expect_equal(none$overdose$pct_trials_any, 0)
  expect_equal(none$overdose$pct_trials_mtd_above, 0)
})

test_that("the two MTD rules differ only where they should", {
  p <- c(0.16, 0.30, 0.44)
  plain <- oc_3p3(p)
  expanded <- oc_3p3(p, mtd_rule = "expand")

  # Escalation is untouched, so exposure to the doses above the cutoff is the
  # same proportion of trials.
  expect_equal(plain$overdose$pct_trials_any, expanded$overdose$pct_trials_any)

  # The expansion only ever adds patients, and it can only make the selection
  # more conservative.
  expect_gt(expanded$total_n_pts, plain$total_n_pts)
  expect_gte(expanded$percent_no_mtd, plain$percent_no_mtd)
})

test_that("3+3 arguments are checked", {
  expect_error(oc_3p3(c(0.1, 0.2), mtd_rule = "other"), "arg")
  expect_error(oc_3p3(c(0.1, 0.2), overdose_cutoff = 2), "between 0 and 1")
  expect_error(sim_3p3(c(0.1, 0.2), n_trials = 0), "at least 1")
  expect_error(oc_3p3("not a probability"), "numeric vector")
})

test_that("print.oc_3p3 shows the design and the table", {
  oc <- oc_3p3(c(0.30, 0.48, 0.67))

  expect_output(print(oc), "3\\+3 operating characteristics")
  expect_output(print(oc), "exact enumeration")
  expect_output(print(oc), "MTD selected")
  expect_output(print(oc), "DL1")
  expect_output(print(oc), "Trials selecting an MTD there")

  expect_output(print(oc, percent = TRUE), "Patients treated \\(%\\)")
  expect_error(print(oc, percent = "yes"), "TRUE or FALSE")

  simulated <- sim_3p3(c(0.30, 0.48, 0.67), n_trials = 200, seed = 1)
  expect_output(print(simulated), "200 simulated trials")

  returned <- quiet_print(oc)
  expect_false(returned$visible)
  expect_identical(returned$value, oc)
})

test_that("sim_3p3 leaves the random seed of the session alone", {
  set.seed(99)
  before <- .Random.seed
  invisible(sim_3p3(c(0.1, 0.3, 0.5), n_trials = 100, seed = 5))
  expect_identical(.Random.seed, before)
})
