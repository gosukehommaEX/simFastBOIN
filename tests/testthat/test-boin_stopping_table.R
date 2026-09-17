test_that("boin_stopping_table reproduces the safety boundary", {
  bd <- boin_boundary(target = 0.30, max_n = 18, extrasafe = TRUE)
  out <- boin_stopping_table(bd)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  # The rule is not evaluated before three patients, so those columns are gone.
  expect_equal(colnames(out), as.character(3:18))
  expect_equal(as.numeric(out[1, ]), 3:18)
  expect_equal(as.numeric(out[2, ]),
               c(2, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8))
  expect_equal(rownames(out),
               c("Number of evaluable patients treated at the lowest dose level",
                 "Stop the trial if # of DLT >="))
})

test_that("cohort_size keeps only the sample sizes reached at a cohort end", {
  bd <- boin_boundary(target = 0.30, max_n = 18, extrasafe = TRUE)
  out <- boin_stopping_table(bd, cohort_size = 3)

  expect_equal(colnames(out), as.character(seq(3, 18, by = 3)))
  expect_equal(as.numeric(out[2, ]), c(2, 4, 5, 6, 7, 8))
})

test_that("the table matches the boundary object it came from", {
  bd <- boin_boundary(target = 0.25, max_n = 24, extrasafe = TRUE)
  out <- boin_stopping_table(bd)

  defined <- !is.na(bd$b_stop)
  expect_equal(as.numeric(out[1, ]), as.numeric(bd$n[defined]))
  expect_equal(as.numeric(out[2, ]), as.numeric(bd$b_stop[defined]))
})

test_that("the safety boundary is never stricter than the elimination boundary", {
  bd <- boin_boundary(target = 0.30, max_n = 30, extrasafe = TRUE)
  out <- boin_stopping_table(bd)
  n <- as.numeric(out[1, ])

  expect_true(all(as.numeric(out[2, ]) <= bd$b_elim[n], na.rm = TRUE))
})

test_that("boin_stopping_table validates its arguments", {
  expect_error(boin_stopping_table(boin_boundary(0.30, 18)), "extrasafe")
  expect_error(boin_stopping_table(list(a = 1)), "boin_boundary")
  expect_error(
    boin_stopping_table(boin_boundary(0.30, 18, extrasafe = TRUE), cohort_size = 0),
    "at least 1"
  )
  expect_error(
    boin_stopping_table(boin_boundary(0.30, 18, extrasafe = TRUE), cohort_size = 40),
    "no sample size"
  )
})
