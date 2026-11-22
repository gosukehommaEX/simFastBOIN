select_cases <- function() {
  list(
    n_pts = rbind(c(3, 6, 12, 9, 0), c(3, 3, 9, 12, 3), c(6, 9, 9, 6, 0),
                  c(3, 3, 3, 6, 15), c(9, 12, 6, 3, 0), c(6, 0, 0, 0, 0),
                  c(3, 6, 9, 6, 6)),
    n_tox = rbind(c(0, 1, 3, 4, 0), c(0, 0, 2, 5, 2), c(1, 3, 4, 4, 0),
                  c(0, 0, 0, 1, 4), c(2, 5, 3, 3, 0), c(5, 0, 0, 0, 0),
                  c(0, 1, 2, 2, 3))
  )
}

test_that("iboin_select_mtd selects the doses of the rule", {
  x <- select_cases()
  skeleton <- c(0.10, 0.19, 0.30, 0.42, 0.54)

  # Expected doses computed apart from the package (Python, scipy)
  with_prior <- iboin_select_mtd(x$n_pts, x$n_tox, 0.30, skeleton, 3)
  expect_identical(with_prior$mtd, c(3L, 3L, 2L, 5L, 2L, NA, 3L))
  expect_identical(with_prior$reason,
                   c(rep("selected", 5), "lowest_dose_eliminated",
                     "selected"))
  expect_identical(with_prior$trial, 1:7)

  robust <- iboin_select_mtd(x$n_pts, x$n_tox, 0.30, skeleton, 3,
                             robust = TRUE)
  expect_identical(robust$mtd, c(3L, 3L, 2L, 5L, 2L, NA, 4L))

  by_dose <- iboin_select_mtd(x$n_pts, x$n_tox, 0.30, skeleton,
                              c(3, 0, 6, 1, 0))
  expect_identical(by_dose$mtd, c(3L, 3L, 2L, 5L, 1L, NA, 3L))

  data_only <- iboin_select_mtd(x$n_pts, x$n_tox, 0.30, skeleton, 3,
                                with_prior = FALSE)
  expect_identical(data_only$mtd, c(3L, 3L, 2L, 5L, 1L, NA, 4L))
})

test_that("without prior information the selection is that of BOIN", {
  trials <- boin_simulate(
    target = 0.30, p_true = c(0.08, 0.15, 0.31, 0.45, 0.55), n_cohort = 10,
    cohort_size = 3, n_trials = 500, n_earlystop = 12, seed = 77
  )
  skeleton <- c(0.10, 0.19, 0.30, 0.42, 0.54)
  for (xs in c(FALSE, TRUE)) {
    ref <- boin_select_mtd(trials$n_pts, trials$n_tox, 0.30,
                           extrasafe = xs)
    expect_identical(
      iboin_select_mtd(trials$n_pts, trials$n_tox, 0.30, skeleton, 0,
                       extrasafe = xs),
      ref
    )
    expect_identical(
      iboin_select_mtd(trials$n_pts, trials$n_tox, 0.30, skeleton, 3,
                       with_prior = FALSE, extrasafe = xs),
      ref
    )
  }
  ref <- boin_select_mtd(trials$n_pts, trials$n_tox, 0.30,
                         min_mtd_sample = 6)
  expect_identical(
    iboin_select_mtd(trials$n_pts, trials$n_tox, 0.30, skeleton, 0,
                     min_mtd_sample = 6),
    ref
  )
})

test_that("only doses given to patients can be selected", {
  skeleton <- c(0.10, 0.19, 0.30, 0.42, 0.54)
  # The prior MTD is dose 3, but only doses 1 and 2 were tried
  res <- iboin_select_mtd(c(3, 27, 0, 0, 0), c(0, 2, 0, 0, 0), 0.30,
                          skeleton, 12)
  expect_identical(res$mtd, 2L)
  none <- iboin_select_mtd(c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0), 0.30,
                           skeleton, 3)
  expect_identical(none$mtd, NA_integer_)
  expect_identical(none$reason, "no_admissible_dose")
})

test_that("iboin_select_mtd rejects inconsistent data", {
  skeleton <- c(0.10, 0.19, 0.30)
  expect_error(iboin_select_mtd(c(3, 3, 3), c(4, 0, 0), 0.30, skeleton, 3),
               "must not exceed")
  expect_error(iboin_select_mtd(c(3, 3), c(1, 0), 0.30, skeleton, 3),
               "skeleton")
  expect_error(iboin_select_mtd(c(3, 3, 3), c(1, 0), 0.30, skeleton, 3),
               "same dimensions")
  expect_error(iboin_select_mtd(c(3, 3, 3), c(1, 0, 0), 0.30, skeleton, 3,
                                with_prior = "yes"), "with_prior")
  expect_error(iboin_select_mtd(c(3, 3, 3), c(1, 0, 0), 0.30, skeleton,
                                1.5), "prior_ess")
})
