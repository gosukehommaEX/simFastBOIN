test_that("the printed rows collapse the numbers of DLTs", {
  prev <- decision_table_3p3_rows(decision_table_3p3())

  expect_identical(names(prev), c("Stage", "Patients", "DLTs", "Decision"))
  expect_identical(prev$Stage, rep("Dose escalation", 5L))
  expect_identical(prev$Patients, c("3", "3", "3", "6", "6"))
  expect_identical(prev$DLTs, c("0", "1", ">= 2", "<= 1", ">= 2"))
  expect_identical(prev$Decision,
                   c("Escalate", "Treat 3 more patients at this dose",
                     "Stop escalation", "Escalate", "Stop escalation"))

  expand <- decision_table_3p3_rows(decision_table_3p3(mtd_rule = "expand"))
  search <- expand[expand$Stage == "Search for the MTD", ]

  expect_equal(nrow(expand), 8L)
  expect_identical(as.list(expand[1:5, ]), as.list(prev))
  expect_identical(search$Patients, c("3", "6", "6"))
  expect_identical(search$DLTs, c("0", "<= 1", ">= 2"))
  expect_identical(search$Decision,
                   c("Treat 3 more patients at this dose",
                     "Select this dose as the MTD",
                     "Move to the next lower dose"))
})

test_that("print.decision_table_3p3 shows the stages and the MTD rule", {
  prev <- decision_table_3p3()
  out <- capture.output(print(prev))

  expect_match(out[1L], "MTD rule: previous", fixed = TRUE)
  expect_true("Dose escalation" %in% out)
  expect_false("Search for the MTD" %in% out)
  # The notes are wrapped, so they are searched with the lines joined again.
  expect_true(grepl("no MTD is selected when escalation stops at the starting dose",
                    paste(out, collapse = " "), fixed = TRUE))

  expand <- decision_table_3p3(mtd_rule = "expand")
  out <- capture.output(print(expand))

  expect_match(out[1L], "MTD rule: expand", fixed = TRUE)
  expect_true("Search for the MTD" %in% out)
  expect_true(grepl("No MTD is selected when the search moves below the starting dose",
                    paste(out, collapse = " "), fixed = TRUE))
})

test_that("print.decision_table_3p3 returns its argument invisibly", {
  prev <- decision_table_3p3()
  result <- quiet_print(prev)

  expect_false(result$visible)
  expect_identical(result$value, prev)
})

test_that("a table without its columns prints as a data frame", {
  prev <- decision_table_3p3()
  out <- capture.output(print(prev[, c("n", "n_tox")]))

  expect_false(any(grepl("Dose escalation", out, fixed = TRUE)))
  expect_true(any(grepl("n_tox", out, fixed = TRUE)))
})
