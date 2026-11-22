# Reproduction of Zhou, Lee, Wang, Bailey and Yuan (2021) with simFastBOIN
#
# Zhou Y, Lee JJ, Wang S, Bailey S, Yuan Y. Incorporating historical
# information to improve phase I clinical trials. Pharmaceutical Statistics
# 2021; 20: 1017-1034. doi:10.1002/pst.2121
#
# Table 1 gives the iBOIN decision boundaries for the skeleton (0.10, 0.19,
# 0.30, 0.42, 0.54), a prior effective sample size (PESS) of 3 at every dose
# and a target of 0.30, for 3, 6, ..., 30 patients at each of the five doses.
# iboin_boundary() must give these 100 values exactly.
#
# Table S.1 of the supporting information gives, for the ten scenarios of
# Table 2 of the article, five operating characteristics of eight designs.
# The rows BOIN, iBOIN and iBOINR (iBOIN with the robust prior) are
# reproduced here: the percentage of correct selection (PCS), the percentage
# of patients treated at the MTD and above the MTD, the risk of overdosing
# and the risk of poor allocation. The settings are five doses, a target of
# 0.30, 30 patients in cohorts of 3 and a PESS of 3; the trials are not
# stopped early (n_earlystop = 100).
#
# What the article does not state, and how it is read here:
#
#   - The selection of the MTD. The article gives two estimates, with and
#     without the prior information, without saying which one the simulation
#     study uses. The iBOIN and iBOINR rows are judged with the prior
#     information (select_with_prior = TRUE). The runs without it are added
#     as INFO: their PCS is lower by 4 to 7 points in most scenarios and
#     higher by 10 in Scenario 9 of iBOIN, where the prior MTD is wrong.
#   - The denominator of the percentages of patients. The 30 planned
#     patients are used, also for trials that stop early because the lowest
#     dose is eliminated; this agrees with the BOIN row of Scenario 1.
#   - The risk of overdosing. The text defines it as 50% or more of the
#     patients above the MTD; the published values agree with more than 50%
#     (more than 15 of 30 patients), which is used here.
#   - The number of simulated trials of the fixed scenarios. The 2,000
#     trials stated for the random scenarios (Section 4.2) are assumed.
#
# A cell is judged PASS when it lies within four standard errors of the
# difference between a simulation of 2,000 trials and one of n_trials trials,
# plus 0.05 for the rounding, and FAIL otherwise. The standard error of a
# percentage of trials uses the mean of the published and the reproduced
# value; that of a percentage of patients uses the standard deviation
# between the reproduced trials.
#
# In Scenarios 2 and 10 the published percentage of patients above the MTD
# (6.0 and 5.4 for BOIN) is far below the reproduced one (about 28 and 24),
# while the risk of overdosing of the same rows, which also counts the
# patients above the MTD, is reproduced. The published values agree with the
# percentage of patients treated at doses whose true DLT probability exceeds
# 0.40, which the column "alternative" reports; in the other scenarios the
# two definitions select the same doses. These cells are registered in
# inst/validation/explained.csv, and their judgements of FAIL become
# EXPLAINED.
#
# Usage, from the package root:
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Zhou2021.R")
#   res <- run_reproduce_zhou2021()
#
# The judgements are written to inst/reproduce/results/reproduce_Zhou2021.csv.

# True DLT probabilities, skeletons and true MTDs of Table 2.
zhou2021_scenarios <- function() {
  list(
    p_true = rbind(
      c(0.30, 0.42, 0.50, 0.60, 0.65), c(0.15, 0.27, 0.40, 0.50, 0.65),
      c(0.08, 0.15, 0.31, 0.45, 0.55), c(0.09, 0.12, 0.15, 0.30, 0.45),
      c(0.05, 0.08, 0.10, 0.14, 0.30), c(0.09, 0.12, 0.15, 0.30, 0.45),
      c(0.08, 0.15, 0.31, 0.45, 0.55), c(0.08, 0.15, 0.31, 0.45, 0.55),
      c(0.04, 0.08, 0.10, 0.18, 0.27), c(0.08, 0.10, 0.28, 0.40, 0.45)
    ),
    skeleton = rbind(
      c(0.30, 0.42, 0.54, 0.64, 0.73), c(0.19, 0.30, 0.42, 0.54, 0.64),
      c(0.10, 0.19, 0.30, 0.42, 0.54), c(0.04, 0.10, 0.19, 0.30, 0.42),
      c(0.01, 0.04, 0.10, 0.19, 0.30), c(0.01, 0.04, 0.10, 0.19, 0.30),
      c(0.19, 0.30, 0.42, 0.54, 0.64), c(0.01, 0.04, 0.10, 0.19, 0.30),
      c(0.04, 0.09, 0.30, 0.40, 0.45), c(0.30, 0.42, 0.54, 0.64, 0.73)
    ),
    mtd = c(1, 2, 3, 4, 5, 4, 3, 3, 5, 3)
  )
}

# Table S.1: PCS, % patients at MTD, % patients above MTD, risk of
# overdosing and risk of poor allocation, one row per scenario.
zhou2021_table_s1 <- function() {
  mk <- function(...) {
    m <- rbind(...)
    colnames(m) <- zhou2021_quantities()
    m
  }
  list(
    BOIN = mk(
      c(60.0, 59.5, 29.5, 23.8, 10.6), c(51.2, 41.6, 6.0, 23.4, 16.0),
      c(52.4, 35.4, 17.6, 8.3, 19.8), c(50.0, 28.7, 13.8, 1.3, 24.5),
      c(71.1, 34.8, 0.0, 0.0, 17.6), c(50.0, 28.7, 13.8, 1.3, 24.5),
      c(52.4, 35.4, 17.6, 8.3, 19.8), c(52.4, 35.4, 17.6, 8.3, 19.8),
      c(70.3, 33.7, 0.0, 0.0, 22.1), c(53.5, 37.5, 5.4, 14.7, 17.6)
    ),
    iBOIN = mk(
      c(65.4, 66.5, 22.8, 12.6, 3.8), c(58.7, 47.9, 3.6, 10.0, 8.7),
      c(59.7, 41.3, 14.4, 3.4, 11.7), c(60.8, 36.4, 11.7, 0.4, 13.1),
      c(77.8, 42.8, 0.0, 0.0, 8.2), c(59.2, 35.7, 18.3, 3.4, 11.7),
      c(62.7, 33.0, 10.9, 2.2, 15.0), c(56.1, 36.7, 28.6, 19.4, 20.0),
      c(51.7, 25.6, 0.0, 0.0, 34.4), c(66.2, 36.1, 0.6, 3.4, 13.7)
    ),
    iBOINR = mk(
      c(65.4, 66.5, 22.8, 12.6, 3.8), c(58.7, 47.9, 3.6, 10.0, 8.7),
      c(59.5, 37.9, 17.3, 8.9, 15.4), c(58.8, 32.8, 16.0, 3.4, 17.8),
      c(77.8, 42.8, 0.0, 0.0, 8.2), c(59.2, 35.7, 18.3, 3.4, 11.7),
      c(62.7, 33.0, 10.9, 2.2, 15.0), c(56.1, 36.7, 28.6, 19.4, 20.0),
      c(67.9, 36.3, 0.0, 0.0, 22.0), c(66.2, 36.1, 0.6, 3.4, 13.7)
    )
  )
}

zhou2021_quantities <- function() {
  c("PCS", "% patients at MTD", "% patients above MTD",
    "Risk of overdosing", "Risk of poor allocation")
}

# Table 1: escalate if the number of DLTs is at most, de-escalate if it is
# at least, for 3, 6, ..., 30 patients (columns) at doses 1 to 5 (rows).
zhou2021_table1 <- function() {
  list(
    escalate = rbind(
      c(1, 1, 2, 3, 4, 4, 5, 6, 6, 7), c(0, 1, 2, 3, 3, 4, 5, 5, 6, 7),
      c(0, 1, 2, 2, 3, 4, 4, 5, 6, 7), c(0, 1, 1, 2, 3, 3, 4, 5, 6, 6),
      c(0, 0, 1, 2, 2, 3, 4, 5, 5, 6)
    ),
    deescalate = rbind(
      c(2, 3, 4, 5, 7, 8, 9, 10, 11, 12), c(2, 3, 4, 5, 6, 7, 8, 9, 11, 12),
      c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11), c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11),
      c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11)
    )
  )
}

reproduce_zhou2021_table1 <- function() {
  pub <- zhou2021_table1()
  bd <- iboin_boundary(target = 0.30,
                       skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
                       prior_ess = 3, max_n = 30)
  at <- seq(3, 30, by = 3)
  rows <- list()
  for (j in 1:5) {
    for (what in c("escalate", "deescalate")) {
      ours <- if (what == "escalate") bd$b_esc[j, at] else bd$b_deesc[j, at]
      rows[[length(rows) + 1L]] <- data.frame(
        item = paste0("Table 1, dose ", j),
        design = "iBOIN",
        quantity = paste0(if (what == "escalate") "Escalate if DLTs <="
                          else "De-escalate if DLTs >=", ", n = ", at),
        published = pub[[what]][j, ], reproduced = as.numeric(ours),
        tolerance = 0, alternative = NA_real_, judged = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# The five quantities of one run, with the standard errors of their
# differences from published values pub of a simulation of n_pub trials, and
# the percentage of patients at doses with a true DLT probability above 0.40.
zhou2021_summary <- function(oc, mtd, pub, n_pub = 2000, n_planned = 30) {
  tr <- oc$trials
  n_trials <- nrow(tr$n_pts)
  factor_se <- sqrt(1 / n_pub + 1 / n_trials)
  above_doses <- seq_len(ncol(tr$n_pts)) > mtd
  n_above <- rowSums(tr$n_pts[, above_doses, drop = FALSE])
  at <- tr$n_pts[, mtd] / n_planned
  above <- n_above / n_planned
  correct <- !is.na(tr$mtd) & tr$mtd == mtd
  overdosing <- n_above > 0.5 * n_planned
  poor <- tr$n_pts[, mtd] < 6
  value <- 100 * c(mean(correct), mean(at), mean(above), mean(overdosing),
                   mean(poor))
  trial_se <- function(reproduced, published) {
    p <- (reproduced + published / 100) / 2
    100 * sqrt(p * (1 - p)) * factor_se
  }
  se <- c(trial_se(mean(correct), pub[1L]),
          100 * stats::sd(at) * factor_se,
          100 * stats::sd(above) * factor_se,
          trial_se(mean(overdosing), pub[4L]),
          trial_se(mean(poor), pub[5L]))
  high <- oc$p_true > 0.40
  alternative <- 100 * mean(rowSums(tr$n_pts[, high, drop = FALSE]) /
                              n_planned)
  list(value = value, se = se, alternative = alternative)
}

reproduce_zhou2021_table_s1 <- function(n_trials = 10000, seed = 2021) {
  sc <- zhou2021_scenarios()
  pub <- zhou2021_table_s1()
  quantity <- zhou2021_quantities()
  runs <- list(
    list(label = "BOIN", pub = "BOIN", ess = 0, robust = FALSE,
         prior = TRUE, judged = TRUE),
    list(label = "iBOIN", pub = "iBOIN", ess = 3, robust = FALSE,
         prior = TRUE, judged = TRUE),
    list(label = "iBOINR", pub = "iBOINR", ess = 3, robust = TRUE,
         prior = TRUE, judged = TRUE),
    list(label = "iBOIN (MTD from the data alone)", pub = "iBOIN", ess = 3,
         robust = FALSE, prior = FALSE, judged = FALSE),
    list(label = "iBOINR (MTD from the data alone)", pub = "iBOINR",
         ess = 3, robust = TRUE, prior = FALSE, judged = FALSE)
  )
  rows <- list()
  for (s in seq_along(sc$mtd)) {
    for (run in runs) {
      published <- unname(pub[[run$pub]][s, ])
      oc <- sim_iboin(
        target = 0.30, p_true = sc$p_true[s, ], skeleton = sc$skeleton[s, ],
        prior_ess = run$ess, n_cohort = 10, cohort_size = 3,
        n_trials = n_trials, n_earlystop = 100, robust = run$robust,
        select_with_prior = run$prior, keep_trials = TRUE, seed = seed
      )
      sm <- zhou2021_summary(oc, sc$mtd[s], published)
      rows[[length(rows) + 1L]] <- data.frame(
        item = paste0("Scenario ", s), design = run$label,
        quantity = quantity, published = published, reproduced = sm$value,
        tolerance = 4 * sm$se + 0.05,
        alternative = c(NA, NA, sm$alternative, NA, NA),
        judged = run$judged, stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
zhou2021_apply_explained <- function(res, file,
                                     path = file.path("inst", "validation",
                                                      "explained.csv")) {
  res$reason <- ""
  if (!file.exists(path)) return(res)
  ex <- utils::read.csv(path, stringsAsFactors = FALSE)
  ex <- ex[ex$file == file, , drop = FALSE]
  key <- paste(paste(res$item, res$design, sep = ", "), res$quantity,
               sep = "\r")
  hit <- match(key, paste(ex$item, ex$quantity, sep = "\r"))
  use <- !is.na(hit) & res$judgement == "FAIL"
  res$judgement[use] <- "EXPLAINED"
  res$reason[use] <- ex$reason[hit[use]]
  res
}

run_reproduce_zhou2021 <- function(n_trials = 10000, seed = 2021,
                                   out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_iboin", mode = "function") ||
      !exists("iboin_boundary", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- rbind(reproduce_zhou2021_table1(),
               reproduce_zhou2021_table_s1(n_trials = n_trials, seed = seed))
  rownames(res) <- NULL
  within <- abs(res$published - res$reproduced) <= res$tolerance + 1e-12
  res$judgement <- ifelse(!res$judged, "INFO",
                          ifelse(within, "PASS", "FAIL"))
  res$judged <- NULL
  res <- zhou2021_apply_explained(res, "reproduce_Zhou2021.R")
  res$reproduced <- round(res$reproduced, 2)
  res$tolerance <- round(res$tolerance, 2)
  res$alternative <- round(res$alternative, 2)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Zhou2021.csv"),
                   row.names = FALSE)
  res
}
