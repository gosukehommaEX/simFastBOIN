test_that("print.iboin_boundary shows the prior and the boundaries", {
  bd <- iboin_boundary(0.30, c(0.10, 0.19, 0.30, 0.42, 0.54), 3, 30)

  expect_output(print(bd), "iBOIN decision boundaries")
  expect_output(print(bd), "Prior information")
  expect_output(print(bd), "Skeleton")
  expect_output(print(bd), "DL1 escalate if # of DLT <=")
  expect_output(print(bd), "DL5 de-escalate if # of DLT >=")
  expect_output(print(bd), "Eliminate if # of DLT >=")
  expect_output(print(bd), "robust prior    : no")

  returned <- quiet_print(bd)
  expect_false(returned$visible)
  expect_identical(returned$value, bd)
})

test_that("cohort_size keeps the multiples of the cohort size", {
  bd <- iboin_boundary(0.30, c(0.10, 0.19, 0.30, 0.42, 0.54), 3, 30)

  out <- capture.output(print(bd, cohort_size = 3))
  header <- out[grep("^ +3 +6 +9", out)]
  expect_length(header, 1L)
  expect_false(grepl(" 4 ", header))

  # The row of the lowest dose of Table 1 of Zhou et al. (2021)
  row <- out[grep("DL1 escalate", out)]
  expect_match(row, "1 +1 +2 +3 +4 +4 +5 +6 +6 +7")

  expect_output(print(bd, cohort_size = 40), "No sample size in the table")
  expect_error(print(bd, cohort_size = 0), "at least 1")
  expect_silent(try(print(bd, cohort_size = 0), silent = TRUE))
})

test_that("boundaries that cannot be reached are printed as NA", {
  bd <- iboin_boundary(0.10, c(0.04, 0.13, 0.16, 0.25, 0.40), 8, 6)
  out <- capture.output(print(bd))
  expect_match(out[grep("DL4 escalate", out)], "NA +NA +NA +NA +NA +NA")
  expect_output(print(iboin_boundary(0.30, c(0.1, 0.3), 3, 6, robust = TRUE,
                                     extrasafe = TRUE)),
                "safety cutoff")
})
