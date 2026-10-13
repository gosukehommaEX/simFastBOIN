zero_matrix <- function() matrix(0, nrow = 3, ncol = 4)

test_that("print.comb_boin_next shows the data, the boundaries and the decision", {
  n_pts <- rbind(c(6, 6, 0, 0), c(3, 3, 0, 0), c(0, 0, 0, 0))
  n_tox <- rbind(c(0, 1, 0, 0), c(1, 2, 0, 0), c(0, 0, 0, 0))
  nxt <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                        current = c(2, 2))
  out <- capture.output(print(nxt))

  expect_identical(out[1L],
                   "BOIN drug-combination design: dose for the next cohort")
  expect_true(any(grepl("A2B2, 2 DLTs out of 3 patients (0.667)", out,
                        fixed = TRUE)))
  expect_true(any(grepl(paste("escalate if DLTs <= 0, de-escalate if DLTs >= 1,",
                              "eliminate if DLTs >= 3 (3 patients)"),
                        out, fixed = TRUE)))
  expect_true(any(out == "Decision: de-escalate to A1B2"))
  expect_true(any(grepl("Pr(in interval)", out, fixed = TRUE)))
  expect_true(any(grepl("drug A +A1B2 +6 +1 +0.218", out)))
  expect_true(any(grepl("drug B +A2B1 +3 +1 +0.167", out)))
  expect_true(any(out == "Eliminated combinations: none"))

  result <- quiet_print(nxt)
  expect_false(result$visible)
  expect_identical(result$value, nxt)
})

test_that("print.comb_boin_next reports ties, notes and stops", {
  n_pts <- zero_matrix()
  n_pts[1, 1] <- 3
  out <- capture.output(print(
    comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = zero_matrix(),
                   current = c(1, 1))
  ))
  expect_true(any(out == "Decision: escalate to A2B1 or A1B2"))
  expect_true(any(grepl("tie", out, fixed = TRUE)))

  n_pts <- zero_matrix()
  n_pts[1, 1] <- 2
  out <- capture.output(print(
    comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = zero_matrix(),
                   current = c(1, 1))
  ))
  expect_true(any(grepl("no elimination yet (2 patients)", out, fixed = TRUE)))
  expect_true(any(grepl("0 DLTs out of 2 patients", out, fixed = TRUE)))

  top <- comb_boin_next(target = 0.25,
                        n_pts = rbind(c(3, 3, 3, 3), c(3, 3, 3, 3),
                                      c(3, 3, 3, 6)),
                        n_tox = zero_matrix(), current = c(3, 4))
  out <- capture.output(print(top))
  expect_true(any(out == "Decision: stay at A3B4"))
  expect_true(any(grepl("no combination one level up", out, fixed = TRUE)))
  expect_true(any(grepl("drug A +- +- +- +- +- +no", out)))

  n_tox <- zero_matrix()
  n_tox[1, 1] <- 3
  n_pts <- zero_matrix()
  n_pts[1, 1] <- 3
  out <- capture.output(print(
    comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                   current = c(1, 1))
  ))
  expect_true(any(grepl("Decision: stop the trial: the lowest combination is eliminated",
                        out, fixed = TRUE)))
  expect_true(any(grepl("Eliminated combinations: A1B1, A1B2, A1B3, A1B4, A2B1",
                        out, fixed = TRUE)))

  n_pts <- rbind(c(3, 3, 0, 0), c(3, 12, 0, 0), c(0, 0, 0, 0))
  n_tox <- rbind(c(0, 0, 0, 0), c(0, 3, 0, 0), c(0, 0, 0, 0))
  out <- capture.output(print(
    comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
                   current = c(2, 2), n_earlystop = 12)
  ))
  expect_true(any(grepl("12 patients have been treated at A2B2", out,
                        fixed = TRUE)))
  expect_true(any(grepl("comb_boin_select_mtd()", out, fixed = TRUE)))
})
