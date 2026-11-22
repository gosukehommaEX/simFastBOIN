# Comparison of sim_gboin() with the 'UnifiedDoseFinding' package
#
# UnifiedDoseFinding::get_oc_QuasiBOIN() (version 0.1.10) is the reference
# implementation of the gBOIN design of Mu et al. (2019) for a
# quasi-Bernoulli toxicity outcome. It draws the grades of a cohort with
# rmultinom(), sim_gboin() with one uniform variate per patient, so that the
# two cannot be compared trial by trial; they are compared in distribution,
# for the ten scenarios of Table 4 of the article with the normalized
# equivalent toxicity score (scores 0, 1/3, 2/3 and 1, target 0.47 / 1.5), 30
# patients in cohorts of 3 and no early stop.
#
# The reference applies the elimination boundary of the number of DLTs to
# the total score. sim_gboin() does so with elim_rule = "count"; these rows
# are judged. Its default, the posterior probability at the total score
# (elim_rule = "posterior"), eliminates a little earlier; these rows are
# INFO. Every other rule is the same: boundaries lambda_e * n and
# lambda_d * n, de-escalation from an eliminated dose, selection of the MTD
# by the isotonic regression of the BOIN design with the total scores.
#
# A cell (selection percentage or average number of patients of a dose, or
# the percentage of trials without MTD) is judged PASS when the two values
# lie within four standard errors of the difference between two independent
# simulations of n_trials trials, plus 0.005, and FAIL otherwise.
#
# get_oc_QuasiBOIN() prints the dose of every cohort; the output is
# discarded.
#
# Usage, from the package root:
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Mu2019.R")
#   source("inst/validation/compare-with-UnifiedDoseFinding.R")
#   res <- run_compare_with_unifieddosefinding()
#
# The judgements are written to
# inst/validation/results/compare-with-UnifiedDoseFinding.csv.

compare_with_unifieddosefinding <- function(n_trials = 4000, seed = 2019) {
  sc <- mu2019_scenarios()
  score <- c(0, 0.5, 1, 1.5) / 1.5
  target <- 0.47 / 1.5
  factor_se <- sqrt(2 / n_trials)
  rows <- list()
  for (s in seq_along(sc)) {
    utils::capture.output(
      ref <- UnifiedDoseFinding::get_oc_QuasiBOIN(
        target = target, p.true = t(sc[[s]]), score = score, ncohort = 10,
        cohortsize = 3, n.earlystop = 100, startdose = 1, ntrial = n_trials,
        seed = seed + s
      )
    )
    ref_value <- c(ref$selpercent, 100 - sum(ref$selpercent), ref$nptsdose)
    for (rule in c("count", "posterior")) {
      oc <- sim_gboin(
        target = target, p_grade = sc[[s]], score = score, n_cohort = 10,
        cohort_size = 3, n_trials = n_trials, n_earlystop = 100,
        elim_rule = rule, keep_trials = TRUE, seed = seed
      )
      ours <- c(unname(oc$sel_percent), oc$percent_no_mtd,
                unname(oc$n_pts_dose))
      p <- (ours[1:7] + ref_value[1:7]) / 200
      se <- c(100 * sqrt(p * (1 - p)) * factor_se,
              apply(oc$trials$n_pts, 2L, stats::sd) * factor_se)
      rows[[length(rows) + 1L]] <- data.frame(
        item = paste0("Scenario ", s),
        design = paste0("elim_rule = ", rule),
        quantity = c(paste0("Selection (%) DL", 1:6), "No MTD (%)",
                     paste0("Patients DL", 1:6)),
        UnifiedDoseFinding = ref_value, simFastBOIN = ours,
        tolerance = 4 * se + 0.005, judged = rule == "count",
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

run_compare_with_unifieddosefinding <- function(
    n_trials = 4000, seed = 2019, out_dir = "inst/validation/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_gboin", mode = "function") ||
      !exists("mu2019_scenarios", mode = "function")) {
    stop("run devtools::load_all() and source ",
         "inst/reproduce/reproduce_Mu2019.R first", call. = FALSE)
  }
  if (!requireNamespace("UnifiedDoseFinding", quietly = TRUE)) {
    stop("install the 'UnifiedDoseFinding' package first", call. = FALSE)
  }
  res <- compare_with_unifieddosefinding(n_trials = n_trials, seed = seed)
  rownames(res) <- NULL
  within <- abs(res$UnifiedDoseFinding - res$simFastBOIN) <=
    res$tolerance + 1e-12
  res$judgement <- ifelse(!res$judged, "INFO",
                          ifelse(within, "PASS", "FAIL"))
  res$judged <- NULL
  res$UnifiedDoseFinding <- round(res$UnifiedDoseFinding, 3)
  res$simFastBOIN <- round(res$simFastBOIN, 3)
  res$tolerance <- round(res$tolerance, 3)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir,
                                  "compare-with-UnifiedDoseFinding.csv"),
                   row.names = FALSE)
  res
}
