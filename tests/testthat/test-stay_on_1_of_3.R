test_that("the modification raises the de-escalation boundary at three patients", {
  plain <- boin_boundary(target = 0.25, max_n = 18)
  stayed <- boin_boundary(target = 0.25, max_n = 18, stay_on_1_of_3 = TRUE)

  expect_equal(plain$b_deesc[3], 1L)
  expect_equal(stayed$b_deesc[3], 2L)
  expect_true(stayed$stay_on_1_of_3_applied)

  # Nothing else in the boundaries moves.
  expect_equal(plain$b_esc, stayed$b_esc)
  expect_equal(plain$b_elim, stayed$b_elim)
  expect_equal(plain$b_deesc[-3], stayed$b_deesc[-3])
})

test_that("exactly one cell of the decision table changes", {
  plain <- boin_decision_table(target = 0.25, max_n = 30)
  stayed <- boin_decision_table(target = 0.25, max_n = 30,
                                stay_on_1_of_3 = TRUE)

  differs <- which(plain != stayed | xor(is.na(plain), is.na(stayed)),
                   arr.ind = TRUE)
  expect_equal(nrow(differs), 1L)
  expect_equal(rownames(plain)[differs[1, "row"]], "1")
  expect_equal(colnames(plain)[differs[1, "col"]], "3")
  expect_equal(plain["1", "3"], "D")
  expect_equal(stayed["1", "3"], "S")
})

test_that("the modification applies over the expected range of targets", {
  # Below about 0.098 one DLT out of three already eliminates the dose, and the
  # safety rule must not be overridden.
  for (target in c(0.05, 0.09, 0.097)) {
    bd <- boin_boundary(target, max_n = 12, stay_on_1_of_3 = TRUE)
    expect_false(bd$stay_on_1_of_3_applied)
    expect_equal(boin_decision_table(target, 12, stay_on_1_of_3 = TRUE)["1", "3"],
                 "DE")
  }

  # In between it de-escalates by default and now stays.
  for (target in c(0.098, 0.15, 0.25, 0.279)) {
    bd <- boin_boundary(target, max_n = 12, stay_on_1_of_3 = TRUE)
    expect_true(bd$stay_on_1_of_3_applied)
  }

  # Above it the design already stays, so there is nothing to change.
  for (target in c(0.28, 0.30, 0.40)) {
    bd <- boin_boundary(target, max_n = 12, stay_on_1_of_3 = TRUE)
    expect_false(bd$stay_on_1_of_3_applied)
    expect_equal(boin_decision_table(target, 12, stay_on_1_of_3 = TRUE)["1", "3"],
                 "S")
  }

  # And where one DLT out of three escalates, an escalation is not overridden.
  for (target in c(0.45, 0.50)) {
    bd <- boin_boundary(target, max_n = 12, stay_on_1_of_3 = TRUE)
    expect_false(bd$stay_on_1_of_3_applied)
    expect_equal(boin_decision_table(target, 12, stay_on_1_of_3 = TRUE)["1", "3"],
                 "E")
  }
})

test_that("a design that already stays is simulated identically", {
  args <- list(target = 0.30, p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
               n_cohort = 20, cohort_size = 3, n_trials = 500, seed = 11)

  plain <- do.call(sim_boin, args)
  stayed <- do.call(sim_boin, c(args, list(stay_on_1_of_3 = TRUE)))

  expect_identical(plain$sel_percent, stayed$sel_percent)
  expect_identical(plain$n_pts_dose, stayed$n_pts_dose)
  expect_identical(plain$n_tox_dose, stayed$n_tox_dose)
})

test_that("a design that de-escalates on one of three is simulated differently", {
  args <- list(target = 0.25, p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
               n_cohort = 20, cohort_size = 3, n_trials = 1000, seed = 12)

  plain <- do.call(sim_boin, args)
  stayed <- do.call(sim_boin, c(args, list(stay_on_1_of_3 = TRUE)))

  expect_false(isTRUE(all.equal(plain$sel_percent, stayed$sel_percent)))
  # Staying rather than stepping down keeps patients on the higher doses.
  expect_gt(sum(stayed$n_pts_dose[3:5]), sum(plain$n_pts_dose[3:5]))
})

test_that("tightening the de-escalation boundary makes the option relevant", {
  # With lambda_d at 0.33 one DLT out of three de-escalates even at a target of
  # 0.30, which is the situation the modification is meant for.
  p_tox <- boin_p_tox(target = 0.30, lambda_d = 0.33)

  expect_equal(boin_decision_table(0.30, 18, p_tox = p_tox)["1", "3"], "D")
  expect_equal(
    boin_decision_table(0.30, 18, p_tox = p_tox, stay_on_1_of_3 = TRUE)["1", "3"],
    "S"
  )
})

test_that("the option is reported and validated", {
  bd <- boin_boundary(target = 0.25, max_n = 18, stay_on_1_of_3 = TRUE)

  expect_true(bd$stay_on_1_of_3)
  expect_output(print(bd), "1 DLT of 3")
  expect_output(print(bd), "stay")
  expect_error(boin_boundary(0.25, 18, stay_on_1_of_3 = NA), "TRUE or FALSE")
  expect_error(boin_simulate(target = 0.25, p_true = c(0.1, 0.3),
                             n_cohort = 5, cohort_size = 3, n_trials = 10,
                             stay_on_1_of_3 = "yes"), "TRUE or FALSE")
})
