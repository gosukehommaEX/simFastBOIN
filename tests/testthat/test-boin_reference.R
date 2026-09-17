# Cross-implementation checks of the decision boundaries and of MTD selection
# against the BOIN package. The trial simulation itself is checked in
# test-boin_reference_oc.R, which is kept separate to keep each test short.

test_that("boin_boundary reproduces BOIN::get.boundary exactly", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  grid <- expand.grid(target = c(0.20, 0.25, 0.30, 0.40),
                      n_cohort = c(10, 20),
                      cutoff_eli = c(0.90, 0.95))

  for (i in seq_len(nrow(grid))) {
    target <- grid$target[i]
    n_cohort <- grid$n_cohort[i]
    cutoff_eli <- grid$cutoff_eli[i]
    max_n <- n_cohort * 3

    reference <- BOIN::get.boundary(
      target = target, ncohort = n_cohort, cohortsize = 3,
      n.earlystop = max_n, cutoff.eli = cutoff_eli, extrasafe = TRUE
    )
    ours <- boin_boundary(target = target, max_n = max_n,
                          cutoff_eli = cutoff_eli, extrasafe = TRUE)

    info <- paste0("target ", target, ", max_n ", max_n,
                   ", cutoff_eli ", cutoff_eli)
    expect_equal(ours$lambda_e, reference$lambda_e, info = info)
    expect_equal(ours$lambda_d, reference$lambda_d, info = info)
    expect_equal(ours$b_esc, as.integer(reference$full_boundary_tab[2, ]),
                 ignore_attr = TRUE, info = info)
    expect_equal(ours$b_deesc, as.integer(reference$full_boundary_tab[3, ]),
                 ignore_attr = TRUE, info = info)
    expect_equal(ours$b_elim, as.integer(reference$full_boundary_tab[4, ]),
                 ignore_attr = TRUE, info = info)
    expect_equal(ours$b_stop, as.integer(reference$stop_boundary[2, ]),
                 ignore_attr = TRUE, info = info)
  }
})

test_that("boin_select_mtd reproduces BOIN::select.mtd exactly", {
  skip_if_not_installed("BOIN")
  expect_true(have_package("BOIN"))

  set.seed(314)
  n_doses <- 5
  target <- 0.30

  for (i in 1:120) {
    n <- sample(0:18, n_doses, replace = TRUE)
    y <- vapply(n, function(k) if (k == 0) 0L else sample(0:k, 1), integer(1))
    extrasafe <- i %% 2 == 0
    bound_mtd <- i %% 3 == 0

    reference <- BOIN::select.mtd(
      target = target, npts = n, ntox = y, cutoff.eli = 0.95,
      extrasafe = extrasafe, offset = 0.05, boundMTD = bound_mtd,
      p.tox = 1.4 * target
    )$MTD

    ours <- boin_select_mtd(
      n_pts = n, n_tox = y, target = target, cutoff_eli = 0.95,
      extrasafe = extrasafe, offset = 0.05, bound_mtd = bound_mtd
    )$mtd

    expect_equal(if (is.na(ours)) 99 else ours, reference,
                 info = paste("draw", i, ": n =", paste(n, collapse = "/"),
                              ", y =", paste(y, collapse = "/")))
  }
})
