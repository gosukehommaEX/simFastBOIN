scen <- list(
  "Scenario 1" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                      p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
  "Scenario 2" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                      p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
)

test_that("the values are those of the design results", {
  b12 <- sim_boin12_multi(target = 0.35, target_resp = 0.25,
                          scenarios = scen, n_cohort = 12, cohort_size = 3,
                          n_trials = 100, seed = 1)
  boin <- sim_boin_multi(target = 0.35, scenarios = scen, n_cohort = 12,
                         cohort_size = 3, n_trials = 100, seed = 1)
  res <- oc_results(BOIN12 = b12, BOIN = boin)
  expect_s3_class(res, "oc_results")
  expect_identical(levels(res$design), c("BOIN12", "BOIN"))
  expect_identical(levels(res$scenario), c("Scenario 1", "Scenario 2"))
  expect_identical(attr(res, "type"), "single")
  expect_identical(attr(res, "dims"), 5L)
  pick <- function(d, s, m) {
    res[res$design == d & res$scenario == s & res$metric == m, ]
  }
  for (s in names(scen)) {
    v <- pick("BOIN12", s, "sel_percent")
    expect_identical(v$dose, 1:5)
    expect_identical(v$value, unname(b12$results[[s]]$sel_percent))
    expect_identical(pick("BOIN12", s, "n_resp")$value,
                     unname(b12$results[[s]]$n_resp_dose))
    expect_identical(pick("BOIN", s, "n_pts")$value,
                     unname(boin$results[[s]]$n_pts_dose))
    expect_identical(pick("BOIN12", s, "percent_no_selection")$value,
                     b12$results[[s]]$percent_no_obd)
    expect_identical(pick("BOIN", s, "percent_no_selection")$value,
                     boin$results[[s]]$percent_no_mtd)
    expect_identical(pick("BOIN", s, "total_n_pts")$value,
                     boin$results[[s]]$total_n_pts)
    expect_true(is.na(pick("BOIN", s, "total_n_pts")$dose))
  }
  # BOIN records no responses and no design here has time.
  expect_identical(nrow(pick("BOIN", "Scenario 1", "n_resp")), 0L)
  expect_false("duration" %in% levels(res$metric))
  expect_identical(unique(res$target_type[res$design == "BOIN"]), "MTD")
  expect_identical(unique(res$target_type[res$design == "BOIN12"]), "OBD")
})

test_that("single results, 3+3 lists and durations are read", {
  p <- c(0.10, 0.25, 0.40, 0.55)
  tite <- sim_tite_boin(target = 0.30, p_true = p, n_cohort = 6,
                        cohort_size = 3, window = 1, accrual_rate = 3,
                        n_trials = 100, seed = 2)
  res <- oc_results(tite, "3+3" = list(A = oc_3p3(p), B = oc_3p3(p / 2)))
  expect_identical(levels(res$design), c("TITE-BOIN", "3+3"))
  expect_identical(levels(res$scenario), c("Scenario 1", "A", "B"))
  dur <- res$value[res$design == "TITE-BOIN" & res$metric == "duration"]
  expect_identical(dur, tite$duration_mean)
  expect_identical(
    res$value[res$design == "3+3" & res$scenario == "B" &
                res$metric == "sel_percent"],
    oc_3p3(p / 2)$sel_percent
  )
  expect_identical(oc_results(tite)$design[1], factor("TITE-BOIN"))
})

test_that("combinations are read cell by cell", {
  p <- matrix(c(0.05, 0.10, 0.20, 0.10, 0.20, 0.30), nrow = 3)
  comb <- sim_comb_boin(target = 0.30, p_true = p, n_cohort = 6,
                        cohort_size = 3, n_trials = 100, seed = 3)
  res <- oc_results(comb)
  expect_identical(attr(res, "type"), "combination")
  expect_identical(attr(res, "dims"), c(3L, 2L))
  sel <- res[res$metric == "sel_percent", ]
  expect_identical(sel$value, as.vector(unname(comb$sel_percent)))
  expect_identical(sel$value[sel$dose_a == 2 & sel$dose_b == 2],
                   unname(comb$sel_percent[2, 2]))
  expect_identical(levels(res$design), "BOIN combination")
})

test_that("invalid inputs are rejected", {
  p <- c(0.10, 0.25, 0.40, 0.55)
  m <- matrix(c(0.05, 0.10, 0.20, 0.10, 0.20, 0.30), nrow = 3)
  boin <- sim_boin(target = 0.30, p_true = p, n_cohort = 4, cohort_size = 3,
                   n_trials = 20, seed = 1)
  comb <- sim_comb_boin(target = 0.30, p_true = m, n_cohort = 4,
                        cohort_size = 3, n_trials = 20, seed = 1)
  expect_error(oc_results(), "at least one design")
  expect_error(oc_results(boin, boin), "must be unique")
  expect_error(oc_results(boin, comb), "single agent")
  expect_error(oc_results(boin, "3+3" = oc_3p3(c(p, 0.7))), "same number")
  expect_error(oc_results(list(1, 2)), "result of a simulation function")
  expect_error(oc_results(data.frame(a = 1)), "result of a simulation")
})
