test_that("boin_boundary reproduces the published rule for cohorts of three", {
  # With a target of 0.30 the design escalates on 0 of 3, stays on 1 of 3 and
  # de-escalates on 2 or more of 3, which is the rule given by Liu and Yuan.
  bd <- boin_boundary(target = 0.30, max_n = 18)
  at_cohort_end <- seq(3, 18, by = 3)

  expect_equal(bd$b_esc[at_cohort_end], c(0, 1, 2, 2, 3, 4))
  expect_equal(bd$b_deesc[at_cohort_end], c(2, 3, 4, 5, 6, 7))
  expect_equal(bd$b_elim[at_cohort_end], c(3, 4, 5, 7, 8, 9))
})

test_that("boin_boundary is internally consistent", {
  bd <- boin_boundary(target = 0.30, max_n = 30, extrasafe = TRUE)

  expect_s3_class(bd, "boin_boundary")
  expect_length(bd$b_esc, 30)
  expect_true(all(bd$b_esc < bd$b_deesc))
  expect_true(all(diff(bd$b_esc) >= 0))
  expect_true(all(diff(bd$b_deesc) >= 0))

  # No elimination decision is possible before three patients.
  expect_true(all(is.na(bd$b_elim[1:2])))
  expect_true(all(is.na(bd$b_stop[1:2])))

  # The de-escalation boundary never exceeds the elimination boundary, and the
  # safety boundary is never stricter than the elimination boundary.
  known <- !is.na(bd$b_elim)
  expect_true(all(bd$b_deesc[known] <= bd$b_elim[known]))
  expect_true(all(bd$b_stop[known] <= bd$b_elim[known], na.rm = TRUE))
})

test_that("boin_boundary agrees with the posterior probability it encodes", {
  target <- 0.30
  cutoff <- 0.95
  bd <- boin_boundary(target = target, max_n = 24, cutoff_eli = cutoff)

  for (n in 3:24) {
    boundary <- bd$b_elim[n]
    if (is.na(boundary)) {
      expect_true(all(1 - pbeta(target, 1:n + 1, n - 1:n + 1) <= cutoff))
    } else {
      expect_gt(1 - pbeta(target, boundary + 1, n - boundary + 1), cutoff)
      if (boundary > 1) {
        expect_lte(1 - pbeta(target, boundary, n - boundary + 2), cutoff)
      }
    }
  }
})

test_that("b_stop is only produced when extrasafe is requested", {
  expect_true(all(is.na(boin_boundary(0.30, 18)$b_stop)))
  expect_false(all(is.na(boin_boundary(0.30, 18, extrasafe = TRUE)$b_stop)))
})

test_that("as.data.frame returns one row per sample size", {
  bd <- boin_boundary(target = 0.30, max_n = 12, extrasafe = TRUE)
  out <- as.data.frame(bd)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 12)
  expect_true("stop_if_dlt_geq" %in% names(out))
  expect_false("stop_if_dlt_geq" %in% names(as.data.frame(boin_boundary(0.30, 12))))
})

test_that("boin_boundary validates its arguments", {
  expect_error(boin_boundary(0.30, max_n = 0), "at least 1")
  expect_error(boin_boundary(0.30, max_n = 2.5), "whole number")
  expect_error(boin_boundary(0.30, 18, offset = 0.6), "between 0 and 0.5")
  expect_error(boin_boundary(0.30, 18, extrasafe = NA), "TRUE or FALSE")
})
