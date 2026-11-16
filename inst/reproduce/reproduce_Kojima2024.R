# Reproduction of Kojima, Wende and Zhao (2024) with simFastBOIN
#
# Kojima M, Wende W, Zhao H. Tips for accelerating BOIN design. Therapeutic
# Innovation and Regulatory Science 2024; 58: 1129-1137.
# doi:10.1007/s43441-024-00692-9
#
# Tables 1 (target 0.30) and 2 (target 0.25) give, for six scenarios each, the
# percentage of trials selecting each dose as the MTD (%MTD), the average
# sample size (Ave SS) and the percentage of trials stopped early, including
# the safety stopping (%ES), of the accelerated BOIN design (ABOIN), the BOIN
# design and the 3+3 design. The settings of the simulation study are five
# doses, at most 24 patients, cohorts of 3, the elimination rule of Liu and
# Yuan (2015) and 10,000 trials. BOIN stops when 12 patients are at the
# current dose and the decision is to stay; ABOIN stops at 6, and treats one
# patient at doses 1 and 2 until the first DLT (titration = TRUE,
# titration_doses = 2). The DLT probabilities are those of Supplemental
# Table 1. The values were read from the text of the article and the
# Supplemental Material.
#
# ABOIN and BOIN are simulated with 10,000 trials and a cell is judged PASS
# when it lies within four standard errors of the difference between two
# independent simulations of 10,000 trials plus 0.005 for the rounding, and
# FAIL otherwise. The %MTD and %ES of 3+3 are multiples of 0.1, which points
# to 1,000 simulated trials. They are compared with the exact values of
# oc_3p3() with mtd_rule = "expand", which gives the average sample size of
# 18 in Scenario 1 (three patients at each dose and three more at the highest),
# within four standard errors of 1,000 trials plus 0.05; these rows are INFO,
# because the article does not state the rule. The average study duration
# (Ave SD) is not reproduced: the article states an accrual of one patient
# per 30 days and a DLT window of 84 days but not how the duration was
# computed.
#
# Usage, from the package root:
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Kojima2024.R")
#   res <- run_reproduce_kojima2024()
#
# The judgements are written to inst/reproduce/results/reproduce_Kojima2024.csv.

kojima2024_scenarios <- function() {
  list(
    "0.30" = list(
      c(0.002, 0.005, 0.01, 0.02, 0.03),
      c(0.06, 0.08, 0.11, 0.19, 0.32),
      c(0.08, 0.12, 0.18, 0.30, 0.41),
      c(0.10, 0.20, 0.31, 0.50, 0.70),
      c(0.15, 0.29, 0.43, 0.50, 0.57),
      c(0.28, 0.41, 0.50, 0.60, 0.70)
    ),
    "0.25" = list(
      c(0.002, 0.005, 0.01, 0.02, 0.03),
      c(0.05, 0.07, 0.09, 0.16, 0.27),
      c(0.07, 0.10, 0.15, 0.25, 0.34),
      c(0.08, 0.17, 0.26, 0.42, 0.58),
      c(0.13, 0.24, 0.36, 0.42, 0.48),
      c(0.23, 0.34, 0.42, 0.50, 0.58)
    )
  )
}

# Rows: %MTD at doses 1 to 5, Ave SS and %ES of each scenario.
kojima2024_tables <- function() {
  mk <- function(...) {
    m <- rbind(...)
    colnames(m) <- c(paste0("MTD_DL", 1:5), "SS", "ES")
    m
  }
  list(
    "0.30" = list(
      ABOIN = mk(
        c(0.01, 0.01, 0.12, 0.49, 99.37, 14.28, 99.94),
        c(1.63, 3.01, 12.45, 32.53, 50.33, 16.07, 88.40),
        c(3.01, 8.20, 30.89, 37.77, 20.05, 15.53, 87.58),
        c(7.47, 27.71, 49.83, 14.39, 0.44, 14.49, 89.31),
        c(21.16, 44.08, 27.53, 5.91, 0.73, 13.95, 89.07),
        c(51.91, 27.70, 13.70, 1.24, 0.04, 12.22, 92.70)
      ),
      BOIN = mk(
        c(0.00, 0.00, 0.00, 0.05, 99.95, 24.00, 0.00),
        c(0.32, 1.66, 8.89, 36.57, 52.53, 23.99, 0.03),
        c(0.89, 7.12, 29.32, 40.90, 21.64, 23.98, 0.13),
        c(5.28, 32.84, 46.49, 14.18, 0.89, 23.94, 0.32),
        c(20.57, 52.57, 20.74, 4.39, 0.65, 23.81, 1.08),
        c(55.45, 25.61, 5.42, 0.55, 0.03, 22.11, 12.94)
      ),
      "3+3" = mk(
        c(0.00, 0.10, 0.80, 1.20, 97.90, 18.30, 0.00),
        c(6.10, 9.00, 23.80, 37.30, 20.20, 18.70, 3.60),
        c(13.30, 21.60, 32.00, 20.40, 5.90, 16.70, 6.80),
        c(31.20, 35.50, 20.70, 3.10, 0.10, 13.70, 9.40),
        c(41.80, 29.90, 7.50, 0.90, 0.20, 11.80, 19.70),
        c(36.90, 9.40, 1.40, 0.30, 0.00, 8.60, 52.00)
      )
    ),
    "0.25" = list(
      ABOIN = mk(
        c(0.01, 0.00, 0.08, 0.31, 99.60, 14.75, 99.52),
        c(1.34, 3.59, 12.56, 33.79, 48.68, 20.27, 61.59),
        c(3.15, 10.05, 31.44, 35.52, 19.80, 20.48, 56.74),
        c(8.58, 34.77, 42.80, 12.75, 0.97, 19.19, 65.41),
        c(23.10, 49.53, 21.02, 4.88, 1.06, 17.91, 70.04),
        c(54.04, 32.53, 8.08, 1.19, 0.09, 14.79, 84.13)
      ),
      BOIN = mk(
        c(0.00, 0.00, 0.12, 0.49, 99.39, 24.00, 0.00),
        c(0.43, 4.16, 16.40, 38.22, 40.76, 23.99, 0.03),
        c(1.74, 13.35, 33.53, 34.91, 16.43, 23.99, 0.04),
        c(6.71, 39.74, 39.91, 12.49, 1.00, 23.97, 0.15),
        c(24.22, 51.55, 18.39, 4.41, 0.83, 23.90, 0.60),
        c(60.65, 23.30, 5.87, 0.99, 0.06, 22.79, 9.13)
      ),
      "3+3" = mk(
        c(0.00, 0.10, 0.80, 1.20, 97.90, 18.30, 0.00),
        c(4.90, 6.20, 19.20, 32.60, 35.30, 19.10, 1.80),
        c(10.10, 16.00, 30.20, 25.70, 13.10, 17.60, 4.90),
        c(24.50, 31.60, 28.70, 7.10, 0.50, 14.80, 7.60),
        c(33.50, 33.90, 11.30, 3.90, 0.40, 12.90, 17.00),
        c(37.90, 17.40, 5.20, 0.90, 0.00, 10.20, 38.60)
      )
    )
  )
}

# %MTD by dose, Ave SS and %ES of one BOIN or ABOIN run, with the standard
# errors of the difference from a published simulation of n_pub trials.
kojima2024_summary <- function(oc, n_pub = 10000) {
  tr <- oc$trials
  n_trials <- nrow(tr$n_pts)
  factor_se <- sqrt(1 / n_pub + 1 / n_trials)
  pct_se <- function(p) 100 * sqrt(p * (1 - p)) * factor_se
  sel <- vapply(1:5, function(j) mean(!is.na(tr$mtd) & tr$mtd == j),
                numeric(1))
  total <- rowSums(tr$n_pts)
  early <- mean(!tr$stop_reason %in% c("max_sample_size", "max_cohorts"))
  value <- c(100 * sel, mean(total), 100 * early)
  se <- c(pct_se(sel), stats::sd(total) * factor_se, pct_se(early))
  cbind(value = value, se = se)
}

reproduce_kojima2024_tables <- function(n_trials = 10000, seed = 2024) {
  scenarios <- kojima2024_scenarios()
  pub <- kojima2024_tables()
  quantity <- c(paste0("MTD_DL", 1:5), "SS", "ES")
  rows <- list()
  for (tt in names(scenarios)) {
    target <- as.numeric(tt)
    for (s in seq_along(scenarios[[tt]])) {
      p <- scenarios[[tt]][[s]]
      for (design in c("ABOIN", "BOIN")) {
        accelerated <- design == "ABOIN"
        oc <- sim_boin(
          target = target, p_true = p, n_cohort = 8, cohort_size = 3,
          n_trials = n_trials, n_earlystop = if (accelerated) 6 else 12,
          titration = accelerated,
          titration_doses = if (accelerated) 2 else NULL,
          keep_trials = TRUE, seed = seed + s
        )
        sm <- kojima2024_summary(oc)
        rows[[length(rows) + 1L]] <- data.frame(
          item = paste("Scenario", s), table = paste("target", tt),
          design = design, quantity = quantity,
          published = unname(pub[[tt]][[design]][s, ]),
          reproduced = sm[, "value"], tolerance = 4 * sm[, "se"] + 0.005,
          judged = TRUE, stringsAsFactors = FALSE
        )
      }
      ex <- oc_3p3(p, mtd_rule = "expand")
      exact <- c(unname(ex$sel_percent), ex$total_n_pts, ex$percent_no_mtd)
      prob <- c(exact[1:5], exact[7]) / 100
      se <- 100 * sqrt(prob * (1 - prob) / 1000)
      rows[[length(rows) + 1L]] <- data.frame(
        item = paste("Scenario", s), table = paste("target", tt),
        design = "3+3", quantity = quantity,
        published = unname(pub[[tt]][["3+3"]][s, ]),
        reproduced = exact, tolerance = c(4 * se[1:5], NA, 4 * se[6]) + 0.05,
        judged = FALSE, stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

run_reproduce_kojima2024 <- function(n_trials = 10000, seed = 2024,
                                     out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_boin", mode = "function") ||
      !"titration_doses" %in% names(formals(sim_boin))) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- reproduce_kojima2024_tables(n_trials = n_trials, seed = seed)
  within <- abs(res$published - res$reproduced) <= res$tolerance + 1e-12
  res$judgement <- ifelse(!res$judged, "INFO",
                          ifelse(within, "PASS", "FAIL"))
  res$judged <- NULL
  res$reproduced <- round(res$reproduced, 2)
  res$tolerance <- round(res$tolerance, 2)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Kojima2024.csv"),
                   row.names = FALSE)
  res
}
