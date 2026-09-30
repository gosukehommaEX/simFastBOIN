# The published tables are read from a transcription of Table 1 of Yuan et al.
# (2018) and Table S1 of its supplementary appendix. Each row is one row of the
# published table; NA in n_tox_max or n_pending_max stands for a range that
# runs up to the largest attainable value, as the signs in the tables do.
read_published <- function() {
  read.csv(test_path("fixtures", "yuan2018-tables.csv"), stringsAsFactors = FALSE)
}

# Integer range that is empty when the upper end is below the lower end.
span <- function(from, to) if (to >= from) from:to else integer(0)

# One row per state covered by the published rows.
expand_published <- function(rows) {
  cells <- lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    tox_max <- if (is.na(r$n_tox_max)) r$n else r$n_tox_max
    per_tox <- lapply(span(r$n_tox_min, tox_max), function(s) {
      pending_max <- r$n - s
      if (!is.na(r$n_pending_max)) pending_max <- min(pending_max, r$n_pending_max)
      pending <- span(r$n_pending_min, pending_max)
      m <- length(pending)
      data.frame(n = rep(r$n, m), n_tox = rep(s, m), n_pending = pending,
                 decision = rep(r$decision, m),
                 esc_bound = rep(r$stft_escalate, m),
                 deesc_bound = rep(r$stft_deescalate, m),
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, per_tox)
  })
  do.call(rbind, cells)
}

# Equal as numbers, or both missing.
same_value <- function(a, b) {
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < 1e-9)
}

test_that("Table 1 and Table S1 of Yuan et al. (2018) are reproduced", {
  published <- read_published()
  code <- c(E = "E", S = "S", D = "D", DE = "DE", SUS = "SUS",
            ES = "E/S", SD = "S/D")
  # Number of cells carrying a boundary in each table, as counted in the PDFs,
  # so that the comparison below cannot pass vacuously.
  n_boundaries <- c("Table 1" = 17L, "Table S1" = 38L)

  for (source in unique(published$source)) {
    rows <- published[published$source == source, ]
    tab <- tite_boin_decision_table(target = rows$target[1L], max_n = 15)
    expected <- expand_published(rows)

    key_tab <- paste(tab$n, tab$n_tox, tab$n_pending)
    key_exp <- paste(expected$n, expected$n_tox, expected$n_pending)

    # Every attainable state at the end of a cohort of three is covered once:
    # 10 + 28 + 55 + 91 + 136 = 320 states for n = 3, 6, 9, 12 and 15.
    expect_equal(length(key_exp), 320L, info = source)
    expect_false(anyDuplicated(key_exp) > 0L, info = source)
    expect_setequal(key_exp, key_tab[tab$n %% 3L == 0L])
    expect_equal(sum(!is.na(expected$esc_bound) | !is.na(expected$deesc_bound)),
                 n_boundaries[[source]], info = source)

    # The states that disagree, listed by name; none are expected.
    k <- match(key_exp, key_tab)
    state <- paste0(source, ": n = ", expected$n, ", DLTs = ", expected$n_tox,
                    ", pending = ", expected$n_pending)
    wrong_decision <- tab$decision[k] != unname(code[expected$decision])
    wrong_esc <- !same_value(round(tab$esc_bound[k], 2), expected$esc_bound)
    wrong_deesc <- !same_value(round(tab$deesc_bound[k], 2), expected$deesc_bound)
    expect_identical(state[wrong_decision], character(0))
    expect_identical(state[wrong_esc], character(0))
    expect_identical(state[wrong_deesc], character(0))
  }
})

test_that("elimination counts the pending patients as treated", {
  # Nine patients, three DLTs among the three completed, six pending. With all
  # nine patients Pr(p > 0.2) is 0.879 and the dose is only de-escalated; with
  # the three completed patients alone it would be 0.998 and eliminated.
  tab <- tite_boin_decision_table(target = 0.2, max_n = 9)
  k <- which(tab$n == 9 & tab$n_tox == 3 & tab$n_pending == 6)

  expect_identical(tab$decision[k], "D")
})

test_that("de-escalation is possible once the observed rate reaches the target", {
  # Fifteen patients with three DLTs: 3 / 15 equals the target of 0.2, and
  # Table 1 of Yuan et al. (2018) de-escalates when STFT is at most 1.16.
  tab <- tite_boin_decision_table(target = 0.2, max_n = 15)
  k <- which(tab$n == 15 & tab$n_tox == 3 & tab$n_pending == 3)

  expect_identical(tab$decision[k], "S/D")
  expect_equal(tab$deesc_bound[k], 1.1575251082673235, tolerance = 1e-10)
})

test_that("states without pending patients follow the BOIN decision table", {
  settings <- list(
    list(target = 0.20, p_saf = NULL, p_tox = NULL),
    list(target = 0.25, p_saf = NULL, p_tox = NULL),
    list(target = 0.30, p_saf = NULL, p_tox = NULL),
    list(target = 0.35, p_saf = NULL, p_tox = NULL),
    list(target = 0.30, p_saf = 0.20, p_tox = 0.40)
  )
  for (set in settings) {
    ref <- unclass(boin_decision_table(target = set$target, max_n = 18,
                                       p_saf = set$p_saf, p_tox = set$p_tox))
    for (method in c("imputation", "ess")) {
      tab <- tite_boin_decision_table(target = set$target, max_n = 18,
                                      method = method,
                                      p_saf = set$p_saf, p_tox = set$p_tox)
      keep <- tab$n_pending == 0L
      expect_identical(tab$decision[keep],
                       ref[cbind(tab$n_tox[keep] + 1L, tab$n[keep])],
                       info = paste(set$target, method))
      expect_true(all(is.na(tab$esc_bound[keep])))
      expect_true(all(is.na(tab$deesc_bound[keep])))
    }
  }
})

test_that("the effective sample size method agrees with independent values", {
  # Reference values computed independently in Python for target 0.3:
  # lambda_e = 0.23649068523646805 and lambda_d = 0.35851946464092954, so that
  # n_tox / lambda_e = 4.228496 n_tox and n_tox / lambda_d = 2.789249 n_tox.
  tab <- tite_boin_decision_table(target = 0.30, max_n = 12, method = "ess")
  expected <- data.frame(
    n = c(3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 9, 9, 9, 12, 12),
    n_tox = c(0, 0, 0, 1, 1, 2, 0, 1, 1, 1, 2, 2, 2, 3, 4, 7),
    n_pending = c(1, 2, 3, 1, 2, 1, 5, 3, 4, 5, 2, 4, 7, 5, 8, 5),
    decision = c("E", "SUS", "SUS", "S/D", "S/D", "D", "SUS", "E/S", "E/S/D",
                 "SUS/S/D", "S/D", "E/S/D", "E/S/D", "S/D", "S/D", "DE"),
    esc_bound = c(NA, NA, NA, NA, NA, NA, NA, 4.228496, 4.228496, 4.228496,
                  NA, 8.456993, 8.456993, NA, NA, NA),
    deesc_bound = c(NA, NA, NA, 2.789249, 2.789249, NA, NA, NA, 2.789249,
                    2.789249, 5.578498, 5.578498, 5.578498, 8.367747,
                    11.156995, NA),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(expected))) {
    e <- expected[i, ]
    k <- which(tab$n == e$n & tab$n_tox == e$n_tox & tab$n_pending == e$n_pending)
    info <- paste("n =", e$n, "DLTs =", e$n_tox, "pending =", e$n_pending)
    expect_length(k, 1L)
    expect_identical(tab$decision[k], e$decision, info = info)
    expect_equal(tab$esc_bound[k], e$esc_bound, tolerance = 1e-6, info = info)
    expect_equal(tab$deesc_bound[k], e$deesc_bound, tolerance = 1e-6, info = info)
  }
})

test_that("the boundaries lie inside the range of the follow-up statistic", {
  for (target in c(0.2, 0.25, 0.3, 0.4)) {
    imp <- tite_boin_decision_table(target = target, max_n = 24)
    ess <- tite_boin_decision_table(target = target, max_n = 24, method = "ess")

    # STFT lies in [0, n_pending).
    e <- !is.na(imp$esc_bound)
    d <- !is.na(imp$deesc_bound)
    expect_true(all(imp$esc_bound[e] > 0 & imp$esc_bound[e] < imp$n_pending[e]))
    expect_true(all(imp$deesc_bound[d] > 0 & imp$deesc_bound[d] < imp$n_pending[d]))

    # ESS lies in [n - n_pending, n).
    done <- ess$n - ess$n_pending
    e <- !is.na(ess$esc_bound)
    d <- !is.na(ess$deesc_bound)
    expect_true(all(ess$esc_bound[e] > done[e] & ess$esc_bound[e] < ess$n[e]))
    expect_true(all(ess$deesc_bound[d] > done[d] & ess$deesc_bound[d] < ess$n[d]))
    expect_true(all(ess$esc_bound[e & d] > ess$deesc_bound[e & d]))

    # A boundary is present exactly when the decision code says so.
    for (tab in list(imp, ess)) {
      expect_identical(!is.na(tab$esc_bound), grepl("^(E|SUS)/", tab$decision))
      expect_identical(!is.na(tab$deesc_bound), grepl("/D$", tab$decision))
    }
  }
})

test_that("the table has the expected shape, vocabulary and attributes", {
  tab <- tite_boin_decision_table(target = 0.30, max_n = 12)
  vocabulary <- c("E", "S", "D", "DE", "SUS", "E/S", "S/D", "E/S/D",
                  "SUS/S", "SUS/S/D")

  expect_s3_class(tab, "tite_boin_decision_table")
  expect_s3_class(tab, "data.frame")
  expect_named(tab, c("n", "n_tox", "n_pending", "decision", "esc_bound",
                      "deesc_bound"))
  # One row per attainable state: (k + 1)(k + 2) / 2 states for k patients.
  expect_equal(nrow(tab), sum((1:12 + 1) * (1:12 + 2) / 2))
  expect_true(all(tab$n_tox + tab$n_pending <= tab$n))
  expect_true(all(tab$decision %in% vocabulary))
  expect_false(is.unsorted(tab$n))

  lambda <- boin_lambda(target = 0.30)
  expect_identical(attr(tab, "method"), "imputation")
  expect_identical(attr(tab, "statistic"), "STFT")
  expect_equal(attr(tab, "lambda_e"), lambda$lambda_e)
  expect_equal(attr(tab, "lambda_d"), lambda$lambda_d)
  expect_equal(attr(tab, "max_pending_ratio"), 0.5)
  expect_identical(attr(tab, "min_completed"), 0L)

  ess <- tite_boin_decision_table(target = 0.30, max_n = 12, method = "ess")
  expect_identical(attr(ess, "statistic"), "ESS")
  expect_equal(attr(ess, "max_pending_ratio"), 1)
  expect_identical(attr(ess, "min_completed"), 2L)
})

test_that("the suspension rules can be changed", {
  find <- function(tab, n, n_tox, n_pending) {
    tab$decision[tab$n == n & tab$n_tox == n_tox & tab$n_pending == n_pending]
  }

  # Two of three patients pending: suspended by default, escalated when up to
  # three quarters may be pending.
  imp <- tite_boin_decision_table(target = 0.30, max_n = 6)
  imp_75 <- tite_boin_decision_table(target = 0.30, max_n = 6,
                                     max_pending_ratio = 0.75)
  expect_identical(find(imp, 3, 0, 2), "SUS")
  expect_identical(find(imp_75, 3, 0, 2), "E")

  # One completed patient blocks escalation under the default of two.
  ess <- tite_boin_decision_table(target = 0.30, max_n = 6, method = "ess")
  ess_0 <- tite_boin_decision_table(target = 0.30, max_n = 6, method = "ess",
                                    min_completed = 0)
  expect_identical(find(ess, 3, 0, 2), "SUS")
  expect_identical(find(ess_0, 3, 0, 2), "E")

  # The pending ratio rule applies to the effective sample size as well.
  ess_50 <- tite_boin_decision_table(target = 0.30, max_n = 6, method = "ess",
                                     max_pending_ratio = 0.5)
  expect_identical(find(ess, 6, 1, 4), "E/S/D")
  expect_identical(find(ess_50, 6, 1, 4), "SUS")

  # A de-escalation that holds whatever the pending outcomes is not suspended.
  expect_identical(find(imp, 3, 1, 2), "SUS")
  imp_20 <- tite_boin_decision_table(target = 0.20, max_n = 3)
  expect_identical(find(imp_20, 3, 1, 2), "D")
})

test_that("invalid arguments are rejected", {
  expect_error(tite_boin_decision_table(0.30, 12, method = "crm"), "should be one of")
  expect_error(tite_boin_decision_table(0.30, 12, max_pending_ratio = 0),
               "max_pending_ratio")
  expect_error(tite_boin_decision_table(0.30, 12, max_pending_ratio = 1.5),
               "max_pending_ratio")
  expect_error(tite_boin_decision_table(0.30, 12, max_pending_ratio = c(0.5, 0.6)),
               "max_pending_ratio")
  expect_error(tite_boin_decision_table(0.30, 12, min_completed = -1),
               "min_completed")
  expect_error(tite_boin_decision_table(0.30, 12, min_completed = 1.5),
               "min_completed")
  expect_error(tite_boin_decision_table(0.30, 0), "max_n")
  expect_error(tite_boin_decision_table(0.02, 12), "target")
})
