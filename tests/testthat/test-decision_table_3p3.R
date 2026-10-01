# Operating characteristics of the 3+3 design obtained by following a decision
# table through every possible course of a trial, with binomial probabilities
# for the DLTs of each group of three patients. It uses nothing of oc_3p3(), so
# agreement with it shows that the table describes the same design.
oc_from_table_3p3 <- function(table, p_true, start_dose) {
  n_doses <- length(p_true)
  rule <- attr(table, "mtd_rule")
  decide <- function(stage, n, n_tox) {
    table$decision[table$stage == stage & table$n == n & table$n_tox == n_tox]
  }

  sel <- numeric(n_doses)
  no_mtd <- 0
  pts <- numeric(n_doses)
  tox <- numeric(n_doses)

  finish <- function(prob, n_pts, n_tox, mtd) {
    pts <<- pts + prob * n_pts
    tox <<- tox + prob * n_tox
    if (is.na(mtd)) {
      no_mtd <<- no_mtd + prob
    } else {
      sel[mtd] <<- sel[mtd] + prob
    }
  }

  search <- function(k, prob, n_pts, n_tox) {
    if (k < start_dose) return(finish(prob, n_pts, n_tox, NA))
    action <- decide("search", n_pts[k], n_tox[k])
    if (action == "S") {
      for (z in 0:3) {
        n_new <- n_pts
        t_new <- n_tox
        n_new[k] <- 6
        t_new[k] <- t_new[k] + z
        p_new <- prob * stats::dbinom(z, 3, p_true[k])
        if (decide("search", 6, t_new[k]) == "MTD") {
          finish(p_new, n_new, t_new, k)
        } else {
          search(k - 1L, p_new, n_new, t_new)
        }
      }
    } else if (action == "MTD") {
      finish(prob, n_pts, n_tox, k)
    } else {
      search(k - 1L, prob, n_pts, n_tox)
    }
  }

  # Escalation stopped at dose d, or passed the highest dose when d exceeds it.
  stopped <- function(d, prob, n_pts, n_tox) {
    if (rule == "previous") {
      finish(prob, n_pts, n_tox, if (d - 1L >= start_dose) d - 1L else NA)
    } else {
      search(d - 1L, prob, n_pts, n_tox)
    }
  }

  escalate <- function(d, prob, n_pts, n_tox) {
    if (d > n_doses) return(stopped(d, prob, n_pts, n_tox))
    for (y in 0:3) {
      n_new <- n_pts
      t_new <- n_tox
      n_new[d] <- 3
      t_new[d] <- y
      p_new <- prob * stats::dbinom(y, 3, p_true[d])
      action <- decide("escalation", 3, y)
      if (action == "S") {
        for (z in 0:3) {
          n_six <- n_new
          t_six <- t_new
          n_six[d] <- 6
          t_six[d] <- y + z
          p_six <- p_new * stats::dbinom(z, 3, p_true[d])
          if (decide("escalation", 6, y + z) == "E") {
            escalate(d + 1L, p_six, n_six, t_six)
          } else {
            stopped(d, p_six, n_six, t_six)
          }
        }
      } else if (action == "E") {
        escalate(d + 1L, p_new, n_new, t_new)
      } else {
        stopped(d, p_new, n_new, t_new)
      }
    }
  }

  escalate(start_dose, 1, numeric(n_doses), numeric(n_doses))
  list(sel_percent = sel * 100, percent_no_mtd = no_mtd * 100,
       n_pts_dose = pts, n_tox_dose = tox)
}

test_that("decision_table_3p3 lists the rules of the 3+3 design", {
  prev <- decision_table_3p3()

  expect_s3_class(prev, "decision_table_3p3")
  expect_identical(attr(prev, "mtd_rule"), "previous")
  expect_identical(unique(prev$stage), "escalation")
  expect_identical(prev$n, c(rep(3L, 4L), rep(6L, 7L)))
  expect_identical(prev$n_tox, c(0:3, 0:6))
  expect_identical(prev$decision,
                   c("E", "S", "STOP", "STOP", "E", "E", rep("STOP", 5L)))

  expand <- decision_table_3p3(mtd_rule = "expand")
  expect_identical(attr(expand, "mtd_rule"), "expand")
  expect_identical(expand$decision[expand$stage == "escalation"], prev$decision)

  search <- expand$stage == "search"
  expect_identical(expand$n[search], c(3L, rep(6L, 7L)))
  expect_identical(expand$n_tox[search], c(0L, 0:6))
  expect_identical(expand$decision[search],
                   c("S", "MTD", "MTD", rep("D", 5L)))

  expect_error(decision_table_3p3(mtd_rule = "other"))
})

test_that("the decision table gives the operating characteristics of oc_3p3", {
  scenarios <- list(
    c(0.30, 0.48, 0.67),
    c(0.16, 0.30, 0.44),
    c(0.05, 0.12, 0.20, 0.30, 0.45),
    c(0.02, 0.05, 0.08, 0.10),
    0.25
  )

  for (rule in c("previous", "expand")) {
    table <- decision_table_3p3(mtd_rule = rule)
    for (p in scenarios) {
      for (start in unique(c(1L, min(2L, length(p))))) {
        walked <- oc_from_table_3p3(table, p, start)
        exact <- oc_3p3(p, mtd_rule = rule, start_dose = start)
        info <- paste(rule, paste(p, collapse = ", "), "start", start)

        expect_equal(walked$sel_percent, as.numeric(exact$sel_percent),
                     tolerance = 1e-10, info = info)
        expect_equal(walked$percent_no_mtd, exact$percent_no_mtd,
                     tolerance = 1e-10, info = info)
        expect_equal(walked$n_pts_dose, as.numeric(exact$n_pts_dose),
                     tolerance = 1e-10, info = info)
        expect_equal(walked$n_tox_dose, as.numeric(exact$n_tox_dose),
                     tolerance = 1e-10, info = info)
      }
    }
  }
})
