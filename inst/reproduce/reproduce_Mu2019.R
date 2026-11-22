# Reproduction of Mu, Yuan, Xu, Mandrekar and Yin (2019) with simFastBOIN
#
# Mu R, Yuan Y, Xu J, Mandrekar SJ, Yin J. gBOIN: a unified model-assisted
# phase I trial design accounting for toxicity grades, and binary or
# continuous end points. Journal of the Royal Statistical Society: Series C
# 2019; 68: 289-308. doi:10.1111/rssc.12263
#
# Table 2 gives the escalation and de-escalation boundaries of the Bernoulli
# and quasi-Bernoulli end point for seven targets, which boin_lambda() must
# give to the three printed decimals. Some of the printed values are cut
# rather than rounded (0.358 for 0.3585), so that a cell is judged PASS
# within 0.001.
#
# Table 5 gives, for the ten scenarios of Table 4, the selection proportion
# and the average number of patients of each of six doses under gBOIN with
# the equivalent toxicity score (ETS: 0 for grades 0 and 1, 0.5, 1 and 1.5
# for grades 2, 3 and 4), a target ETS of 0.47, 30 patients in cohorts of 3
# and 4,000 simulated trials. The trials are not stopped early
# (n_earlystop = 100).
#
# What differs from the article, and what it does not state:
#
#   - The scale. sim_gboin() divides the scores by their maximum of 1.5 and
#     takes the target on that scale, 0.47 / 1.5. The article computes its
#     boundaries (0.373 and 0.566) from the unnormalized target of 0.47; the
#     boundaries used here, multiplied by 1.5, are 0.371 and 0.562. For the
#     ETS, whose totals are multiples of 0.5, the two pairs give the same
#     decision except for a mean ETS of 0.5625 at 24 patients.
#   - The elimination rule. The article describes the posterior probability
#     of a mean score above the target; how it was evaluated in the
#     simulation (scale, whole numbers) is not stated. The judged rows use
#     elim_rule = "posterior"; the rows with elim_rule = "count", the rule of
#     the 'UnifiedDoseFinding' package on the normalized scale, are INFO.
#   - Three columns of Table 4 do not sum to 1 as printed (Scenario 2 dose
#     6: 1.06, Scenario 3 dose 2: 0.99, Scenario 4 dose 6: 1.14). Every
#     column is divided by its sum, as rmultinom() of the reference
#     implementation does.
#
# A cell is judged PASS when it lies within four standard errors of the
# difference between a simulation of 4,000 trials and one of n_trials
# trials, plus half a unit of the last printed digit (0.5 for the selection
# percentages, 0.05 for the numbers of patients), and FAIL otherwise.
#
# Two cells hold misprints and are registered in
# inst/validation/explained.csv: the selection proportion of dose 3 in
# Scenario 9 is printed as 0.00 although the row sums to 0.91 and 4.1
# patients are treated there; the average number of patients of dose 1 in
# Scenario 8 is printed as 2.7 although every trial treats three patients at
# the starting dose. Agreement with Table 5 is otherwise approximate: a
# number of cells differ by more than the allowance, see the results file.
# The rules of the package were checked trial by trial against an
# implementation in Python and, in distribution, against the
# 'UnifiedDoseFinding' package (inst/validation/
# compare-with-UnifiedDoseFinding.R).
#
# Usage, from the package root:
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Mu2019.R")
#   res <- run_reproduce_mu2019()
#
# The judgements are written to inst/reproduce/results/reproduce_Mu2019.csv.

# Table 4: true probabilities of grades 0 or 1, 2, 3 and 4 (rows) at the six
# doses (columns) of each scenario, columns divided by their sums.
mu2019_scenarios <- function() {
  mk <- function(...) {
    m <- rbind(...)
    sweep(m, 2L, colSums(m), "/")
  }
  list(
    mk(c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19),
       c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11),
       c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11),
       c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)),
    mk(c(0.92, 0.85, 0.70, 0.55, 0.24, 0.00),
       c(0.03, 0.05, 0.10, 0.15, 0.26, 0.36),
       c(0.03, 0.07, 0.14, 0.21, 0.35, 0.49),
       c(0.02, 0.03, 0.06, 0.09, 0.15, 0.21)),
    mk(c(0.78, 0.56, 0.50, 0.40, 0.30, 0.16),
       c(0.14, 0.19, 0.18, 0.17, 0.15, 0.09),
       c(0.06, 0.12, 0.14, 0.15, 0.14, 0.10),
       c(0.02, 0.12, 0.18, 0.28, 0.41, 0.65)),
    mk(c(0.88, 0.64, 0.52, 0.35, 0.17, 0.00),
       c(0.04, 0.12, 0.16, 0.22, 0.28, 0.39),
       c(0.06, 0.17, 0.22, 0.30, 0.38, 0.52),
       c(0.02, 0.07, 0.10, 0.13, 0.17, 0.23)),
    mk(c(1.00, 0.91, 0.88, 0.86, 0.80, 0.65),
       c(0.00, 0.06, 0.07, 0.08, 0.10, 0.13),
       c(0.00, 0.03, 0.04, 0.05, 0.08, 0.14),
       c(0.00, 0.00, 0.01, 0.01, 0.02, 0.08)),
    mk(c(0.50, 0.38, 0.29, 0.19, 0.13, 0.08),
       c(0.25, 0.24, 0.21, 0.16, 0.11, 0.07),
       c(0.11, 0.12, 0.12, 0.10, 0.08, 0.05),
       c(0.14, 0.26, 0.38, 0.55, 0.68, 0.80)),
    mk(c(0.78, 0.58, 0.50, 0.40, 0.30, 0.16),
       c(0.14, 0.18, 0.18, 0.17, 0.15, 0.09),
       c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
       c(0.08, 0.24, 0.32, 0.43, 0.55, 0.75)),
    mk(c(0.92, 0.76, 0.68, 0.57, 0.45, 0.25),
       c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
       c(0.08, 0.24, 0.32, 0.43, 0.55, 0.75),
       c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00)),
    mk(c(0.66, 0.10, 0.00, 0.00, 0.00, 0.00),
       c(0.34, 0.90, 0.86, 0.54, 0.20, 0.33),
       c(0.00, 0.00, 0.14, 0.46, 0.80, 0.00),
       c(0.00, 0.00, 0.00, 0.00, 0.00, 0.67)),
    mk(c(0.84, 0.52, 0.36, 0.14, 0.45, 0.25),
       c(0.16, 0.48, 0.64, 0.86, 0.00, 0.00),
       c(0.00, 0.00, 0.00, 0.00, 0.55, 0.75),
       c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00))
  )
}

# gBOIN rows of Table 5: selection percentage of doses 1 to 6 and average
# number of patients at doses 1 to 6, one row per scenario.
mu2019_table5 <- function() {
  list(
    sel = 100 * rbind(
      c(0.00, 0.02, 0.29, 0.56, 0.12, 0.01),
      c(0.00, 0.00, 0.16, 0.66, 0.18, 0.00),
      c(0.02, 0.47, 0.36, 0.13, 0.00, 0.00),
      c(0.01, 0.28, 0.47, 0.23, 0.01, 0.00),
      c(0.00, 0.00, 0.00, 0.00, 0.04, 0.96),
      c(0.69, 0.29, 0.01, 0.00, 0.00, 0.00),
      c(0.11, 0.55, 0.28, 0.06, 0.00, 0.00),
      c(0.00, 0.06, 0.22, 0.41, 0.28, 0.04),
      c(0.00, 0.91, 0.00, 0.00, 0.00, 0.00),
      c(0.00, 0.00, 0.03, 0.71, 0.23, 0.03)
    ),
    pts = rbind(
      c(3.5, 5.1, 8.4, 8.9, 3.7, 0.4),
      c(3.5, 4.4, 7.5, 9.9, 4.3, 0.4),
      c(6.3, 12.1, 7.7, 3.2, 0.5, 0.0),
      c(4.8, 9.0, 10.8, 4.4, 0.9, 0.1),
      c(3.0, 3.1, 3.3, 3.4, 4.6, 12.6),
      c(19.4, 8.0, 1.7, 0.2, 0.0, 0.0),
      c(8.9, 13.0, 6.2, 1.7, 0.2, 0.0),
      c(2.7, 5.6, 7.2, 8.5, 4.4, 1.6),
      c(3.2, 22.7, 4.1, 0.0, 0.0, 0.0),
      c(3.0, 3.5, 5.0, 13.3, 3.8, 1.3)
    )
  )
}

# Table 2, row of the Bernoulli and quasi-Bernoulli end point.
reproduce_mu2019_table2 <- function() {
  target <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
  pub_e <- c(0.078, 0.118, 0.157, 0.197, 0.236, 0.276, 0.316)
  pub_d <- c(0.119, 0.179, 0.238, 0.298, 0.358, 0.419, 0.479)
  ours <- lapply(target, boin_lambda)
  data.frame(
    item = "Table 2", design = "gBOIN",
    quantity = c(paste("lambda_e, target", format(target, nsmall = 2)),
                 paste("lambda_d, target", format(target, nsmall = 2))),
    published = c(pub_e, pub_d),
    reproduced = c(vapply(ours, function(z) z$lambda_e, numeric(1)),
                   vapply(ours, function(z) z$lambda_d, numeric(1))),
    tolerance = 0.001, judged = TRUE, stringsAsFactors = FALSE
  )
}

reproduce_mu2019_table5 <- function(n_trials = 10000, seed = 2019,
                                    n_pub = 4000) {
  sc <- mu2019_scenarios()
  pub <- mu2019_table5()
  factor_se <- sqrt(1 / n_pub + 1 / n_trials)
  rows <- list()
  for (s in seq_along(sc)) {
    for (rule in c("posterior", "count")) {
      oc <- sim_gboin(
        target = 0.47 / 1.5, p_grade = sc[[s]], score = c(0, 0.5, 1, 1.5),
        n_cohort = 10, cohort_size = 3, n_trials = n_trials,
        n_earlystop = 100, elim_rule = rule, keep_trials = TRUE, seed = seed
      )
      label <- if (rule == "posterior") "gBOIN" else
        "gBOIN (elimination by count)"
      p <- (oc$sel_percent + pub$sel[s, ]) / 200
      rows[[length(rows) + 1L]] <- data.frame(
        item = paste0("Scenario ", s), design = label,
        quantity = paste0("Selection (%) DL", 1:6),
        published = pub$sel[s, ], reproduced = unname(oc$sel_percent),
        tolerance = 400 * sqrt(p * (1 - p)) * factor_se + 0.5,
        judged = rule == "posterior", stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- data.frame(
        item = paste0("Scenario ", s), design = label,
        quantity = paste0("Patients DL", 1:6),
        published = pub$pts[s, ], reproduced = unname(oc$n_pts_dose),
        tolerance = 4 * apply(oc$trials$n_pts, 2L, stats::sd) * factor_se +
          0.05,
        judged = rule == "posterior", stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
mu2019_apply_explained <- function(res, file,
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

run_reproduce_mu2019 <- function(n_trials = 10000, seed = 2019,
                                 out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_gboin", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- rbind(reproduce_mu2019_table2(),
               reproduce_mu2019_table5(n_trials = n_trials, seed = seed))
  rownames(res) <- NULL
  within <- abs(res$published - res$reproduced) <= res$tolerance + 1e-12
  res$judgement <- ifelse(!res$judged, "INFO",
                          ifelse(within, "PASS", "FAIL"))
  res$judged <- NULL
  res <- mu2019_apply_explained(res, "reproduce_Mu2019.R")
  res$reproduced <- round(res$reproduced, 3)
  res$tolerance <- round(res$tolerance, 3)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Mu2019.csv"),
                   row.names = FALSE)
  res
}
