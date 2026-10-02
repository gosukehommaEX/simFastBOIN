test_that("the printed rows are those of the published tables", {
  # The transcription of Table 1 and Table S1 of Yuan et al. (2018) holds one
  # row per printed row, in the published order.
  published <- read.csv(test_path("fixtures", "yuan2018-tables.csv"),
                        stringsAsFactors = FALSE)

  for (source in unique(published$source)) {
    rows <- published[published$source == source, ]
    tab <- tite_boin_decision_table(target = rows$target[1L], max_n = 15)
    shown <- tite_boin_display_rows(tab, cohort_size = 3, digits = 2L)

    tox_label <- ifelse(
      is.na(rows$n_tox_max), paste(">=", rows$n_tox_min),
      ifelse(rows$n_tox_max == rows$n_tox_min, as.character(rows$n_tox_min),
             paste0(rows$n_tox_min, ", ", rows$n_tox_max))
    )
    pending_label <- ifelse(
      is.na(rows$n_pending_max), paste(">=", rows$n_pending_min),
      ifelse(rows$n_pending_min == 0 & rows$n_pending_max > 0,
             paste("<=", rows$n_pending_max), as.character(rows$n_pending_min))
    )
    e <- formatC(rows$stft_escalate, format = "f", digits = 2)
    d <- formatC(rows$stft_deescalate, format = "f", digits = 2)
    escalate <- ifelse(rows$decision == "E", "Y",
                ifelse(rows$decision == "SUS", "Suspend",
                ifelse(rows$decision == "ES", paste(">=", e), "")))
    stay <- ifelse(rows$decision == "S", "Y",
            ifelse(rows$decision == "SUS", "Suspend",
            ifelse(rows$decision == "ES", paste("<", e),
            ifelse(rows$decision == "SD", paste(">", d), ""))))
    deescalate <- ifelse(rows$decision == "D", "Y",
                  ifelse(rows$decision == "DE", "Y & Elim",
                  ifelse(rows$decision == "SUS", "Suspend",
                  ifelse(rows$decision == "SD", paste("<=", d), ""))))

    expect_equal(nrow(shown), nrow(rows), info = source)
    expect_identical(shown$Patients, as.character(rows$n), info = source)
    expect_identical(shown$DLTs, tox_label, info = source)
    expect_identical(shown$Pending, pending_label, info = source)
    expect_identical(shown$Escalate, escalate, info = source)
    expect_identical(shown$Stay, stay, info = source)
    expect_identical(shown$`De-escalate`, deescalate, info = source)
  }
})

test_that("print.tite_boin_decision_table shows the method and the rules", {
  imp <- tite_boin_decision_table(target = 0.2, max_n = 15)
  out <- capture.output(print(imp, cohort_size = 3))

  expect_match(out[1L], "target 0.2", fixed = TRUE)
  expect_true(any(grepl("single mean imputation", out, fixed = TRUE)))
  expect_true(any(grepl("boundaries refer to STFT", out, fixed = TRUE)))
  expect_true(any(grepl("more than 50% of the patients are pending", out,
                        fixed = TRUE)))
  expect_false(any(grepl("escalation requires at least", out, fixed = TRUE)))
  expect_false(any(grepl("MF", out, fixed = TRUE)))
  expect_true(any(grepl("^ *9 +1 +4 +>= 2.15 +< 2.15", out)))

  ess <- tite_boin_decision_table(target = 0.3, max_n = 12, method = "ess")
  out <- capture.output(print(ess, cohort_size = 3))

  expect_true(any(grepl("effective sample size", out, fixed = TRUE)))
  expect_true(any(grepl("boundaries refer to ESS", out, fixed = TRUE)))
  expect_true(any(grepl("escalation requires at least 2 patients", out,
                        fixed = TRUE)))
  expect_false(any(grepl("of the patients are pending", out, fixed = TRUE)))
  expect_true(any(grepl("Suspend if >= 4.23", out, fixed = TRUE)))
})

test_that("escalations with pending patients show the condition on MF", {
  # Table A1 of the supplementary materials of Chen et al. (2026).
  tab <- tite_boin_decision_table(target = 0.25, max_n = 9,
                                  max_pending_ratio = 0.49,
                                  min_follow_up = 0.25)
  shown <- tite_boin_display_rows(tab, cohort_size = 3, digits = 2L)
  row_of <- function(n, dlts, pending) {
    shown[shown$Patients == n & shown$DLTs == dlts & shown$Pending == pending, ]
  }

  expect_identical(row_of("6", "1", "1")$Escalate, ">= 0.22 & MF >= 0.25")
  expect_identical(row_of("6", "1", "1")$Stay, "< 0.22")
  expect_identical(row_of("9", "1", "4")$Escalate, ">= 0.66 & MF >= 0.25")
  expect_identical(row_of("9", "0", "1-4")$Escalate, "MF >= 0.25")
  expect_identical(row_of("9", "0", "0")$Escalate, "Y")
  expect_identical(row_of("3", "1, 2", "<= 2")$`De-escalate`, "Y")
  expect_false(any(grepl("MF", shown$Stay, fixed = TRUE)))
  expect_false(any(grepl("MF", shown$`De-escalate`, fixed = TRUE)))

  out <- capture.output(print(tab, cohort_size = 3))
  expect_true(any(grepl("MF = shortest follow-up time", out, fixed = TRUE)))
  expect_true(any(grepl("escalation requires MF >= 0.25", out, fixed = TRUE)))
  expect_true(any(grepl("more than 49% of the patients are pending", out,
                        fixed = TRUE)))
})

test_that("print.tite_boin_decision_table returns its input invisibly", {
  tab <- tite_boin_decision_table(target = 0.3, max_n = 6)
  result <- quiet_print(tab, cohort_size = 3)

  expect_false(result$visible)
  expect_identical(result$value, tab)
})

test_that("cohort_size and digits change the output", {
  tab <- tite_boin_decision_table(target = 0.2, max_n = 12)

  all_n <- tite_boin_display_rows(tab)
  by_three <- tite_boin_display_rows(tab, cohort_size = 3)
  expect_setequal(unique(all_n$Patients), as.character(1:12))
  expect_setequal(unique(by_three$Patients), c("3", "6", "9", "12"))

  two <- capture.output(print(tab, cohort_size = 3))
  three <- capture.output(print(tab, cohort_size = 3, digits = 3))
  expect_true(any(grepl(">= 2.15", two, fixed = TRUE)))
  expect_true(any(grepl(">= 2.151", three, fixed = TRUE)))

  expect_output(print(tab, cohort_size = 25), "No sample size")
})

test_that("invalid arguments are rejected before anything is printed", {
  tab <- tite_boin_decision_table(target = 0.3, max_n = 6)

  expect_output(expect_error(print(tab, cohort_size = 0), "cohort_size"), NA)
  expect_output(expect_error(print(tab, digits = -1), "digits"), NA)
  expect_output(expect_error(print(tab, digits = 1.5), "digits"), NA)
})

test_that("a subset of the table can be printed", {
  tab <- tite_boin_decision_table(target = 0.2, max_n = 9)

  # One state keeps the layout of the table.
  one <- tab[tab$n == 9 & tab$n_tox == 1 & tab$n_pending == 4, ]
  expect_output(print(one), ">= 2.15")

  # Without the columns of a decision table it prints as a data frame.
  expect_output(print(tab[1:3, c("n", "decision")]), "decision")
})
