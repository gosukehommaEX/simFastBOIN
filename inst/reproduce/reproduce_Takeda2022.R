# Reproduction of Takeda, Morita and Taguri (2022) with simFastBOIN
#
# Takeda K, Morita S, Taguri M. gBOIN-ET: The generalized Bayesian optimal
# interval design for optimal dose-finding accounting for ordinal graded
# efficacy and toxicity in early clinical trials. Biometrical Journal 2022;
# 64(7): 1178-1191. doi:10.1002/bimj.202100263
#
# The article prints the boundaries of gBOIN-ET and of BOIN-ET for its
# simulation study (Section 3.1) and for its sensitivity analysis (Section
# 3.3), which boin_et_boundary() must give to the two printed decimals.
#
# Table 1 gives, for the ten scenarios of Table S1 of the supporting
# information, the percentage of trials recommending each of six doses as
# the optimal dose, the percentage of early terminations (ET) and the
# average number of patients at each dose under gBOIN-ET with the equivalent
# toxicity score (0 for grades 0 and 1, 0.5, 1 and 1.5 for grades 2, 3 and
# 4), the equivalent efficacy score (0 for stable or progressive disease, 1
# for a partial and 1.5 for a complete response), targets of 0.47 and 0.95,
# 30 patients in cohorts of 3 and 4,000 simulated trials.
#
# The SAS program of the authors (supporting information, folder Macro) was
# read, and it differs from the article and from the package as follows.
#
#   - The categories. The program does not draw the category of a patient
#     from the probabilities of Table S1. It draws one Bernoulli variate for
#     each category with the tabulated probability and takes the highest
#     category that comes up (macro gendata). With t2, t3 and t4 the
#     tabulated probabilities of grades 2, 3 and 4, the probabilities
#     actually simulated are t4 for grade 4, (1 - t4) t3 for grade 3 and
#     (1 - t4) (1 - t3) t2 for grade 2; with cr and pr those of a complete
#     and a partial response, cr and (1 - cr) pr. The mean scores are lower
#     than those of Table S1 (Scenario 1, dose 4: 0.43 and 0.85 in place of
#     0.47 and 0.95). takeda2022_scenarios("simulated") gives these
#     probabilities, and the judged rows use them.
#   - The boundaries. The program uses 0.14133, 0.37333 and 0.50, which the
#     judged rows take through lambda1, lambda2 and eta1.
#   - The selection of the optimal dose. The program fits the linear
#     fractional polynomial of Section 2.3 (est_method = "fp_linear") and
#     selects among the doses up to the MTD with at least six patients
#     (min_obd_sample = 6); the six patients are not in the article.
#   - The exclusion of doses. The program has none of the posterior
#     probability rules of Section 2.2. A dose at which the design would
#     stay is left for the dose below when its mean efficacy score is at
#     most half the target or its toxicity score reaches a count boundary,
#     the trial ends when this happens at the lowest dose, and no dose is
#     excluded when the optimal dose is selected. The package keeps the
#     posterior probability rules; the judged rows use
#     boinet_compatible = TRUE, under which a dose excluded for toxicity
#     does not exclude the doses above it, which is closer to the program.
#   - The toxicity and the efficacy of a patient are independent
#     (resp_cor = 0).
#
# A cell is judged PASS when it lies within four standard errors of the
# difference between a simulation of 4,000 trials and one of n_trials
# trials, plus half a unit of the last printed digit (0.5 for the
# percentages, 0.05 for the numbers of patients), and FAIL otherwise. The
# rows with the probabilities of Table S1 as tabulated, and those with the
# default rules of the package, are INFO.
#
# The program was written again in Python (02_Rpackage/Validation/gboin_et
# of the development files, not part of the package). That transcription
# gives all percentages of Table 1 but one within the allowance, and ten of
# the sixty numbers of patients outside it: the published trials treat
# somewhat more patients at the higher doses than the program as read, and
# the published numbers of patients of Scenario 6 sum to 26.3 although 3% of
# the trials are reported as terminated early. The cells that the package
# misses and the transcription reproduces are registered in
# inst/validation/explained.csv; the cells that neither reproduces stay
# FAIL.
#
# The rules of the package were checked trial by trial against an
# implementation in Python and, in distribution, against boinet::gboinet()
# (inst/validation/compare-with-gboinet.R).
#
# Usage, from the package root:
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Takeda2022.R")
#   res <- run_reproduce_takeda2022()
#
# The judgements are written to
# inst/reproduce/results/reproduce_Takeda2022.csv.

# Table S1, as entered in the program of the authors: probabilities of
# grades 0 or 1, 2, 3 and 4 and of a complete and a partial response at the
# six doses of each scenario. With type "table" the categories have the
# tabulated probabilities; with "simulated" those that the program draws.
# Three columns of the toxicity grades do not sum to 1 as printed (Scenario
# 2 dose 6: 1.06, Scenario 3 dose 2: 0.99, Scenario 4 dose 6: 1.14); with
# type "table" every column is divided by its sum, as in
# inst/reproduce/reproduce_Mu2019.R. The probabilities that the program
# draws need no such step, because each grade is drawn on its own.
takeda2022_scenarios <- function(type = c("simulated", "table")) {
  type <- match.arg(type)
  tox <- list(
    rbind(c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19),
          c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11),
          c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11),
          c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)),
    rbind(c(0.92, 0.85, 0.70, 0.55, 0.24, 0.00),
          c(0.03, 0.05, 0.10, 0.15, 0.26, 0.36),
          c(0.03, 0.07, 0.14, 0.21, 0.35, 0.49),
          c(0.02, 0.03, 0.06, 0.09, 0.15, 0.21)),
    rbind(c(0.78, 0.56, 0.50, 0.40, 0.30, 0.16),
          c(0.14, 0.19, 0.18, 0.17, 0.15, 0.09),
          c(0.06, 0.12, 0.14, 0.15, 0.14, 0.10),
          c(0.02, 0.12, 0.18, 0.28, 0.41, 0.65)),
    rbind(c(0.88, 0.64, 0.52, 0.35, 0.17, 0.00),
          c(0.04, 0.12, 0.16, 0.22, 0.28, 0.39),
          c(0.06, 0.17, 0.22, 0.30, 0.38, 0.52),
          c(0.02, 0.07, 0.10, 0.13, 0.17, 0.23)),
    rbind(c(1.00, 0.91, 0.88, 0.86, 0.80, 0.65),
          c(0.00, 0.06, 0.07, 0.08, 0.10, 0.13),
          c(0.00, 0.03, 0.04, 0.05, 0.08, 0.14),
          c(0.00, 0.00, 0.01, 0.01, 0.02, 0.08)),
    rbind(c(0.50, 0.38, 0.29, 0.19, 0.13, 0.08),
          c(0.25, 0.24, 0.21, 0.16, 0.11, 0.07),
          c(0.11, 0.12, 0.12, 0.10, 0.08, 0.05),
          c(0.14, 0.26, 0.38, 0.55, 0.68, 0.80)),
    rbind(c(0.78, 0.58, 0.50, 0.40, 0.30, 0.16),
          c(0.14, 0.18, 0.18, 0.17, 0.15, 0.09),
          c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
          c(0.08, 0.24, 0.32, 0.43, 0.55, 0.75)),
    rbind(c(0.92, 0.76, 0.68, 0.57, 0.45, 0.25),
          c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
          c(0.08, 0.24, 0.32, 0.43, 0.55, 0.75),
          c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00)),
    rbind(c(0.66, 0.10, 0.00, 0.00, 0.00, 0.00),
          c(0.34, 0.90, 0.86, 0.54, 0.20, 0.33),
          c(0.00, 0.00, 0.14, 0.46, 0.80, 0.00),
          c(0.00, 0.00, 0.00, 0.00, 0.00, 0.67)),
    rbind(c(0.84, 0.52, 0.36, 0.14, 0.45, 0.25),
          c(0.16, 0.48, 0.64, 0.86, 0.00, 0.00),
          c(0.00, 0.00, 0.00, 0.00, 0.55, 0.75),
          c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00))
  )
  cr <- rbind(
    c(0.05, 0.10, 0.15, 0.50, 0.55, 0.60),
    c(0.10, 0.20, 0.55, 0.45, 0.35, 0.25),
    c(0.20, 0.50, 0.55, 0.60, 0.65, 0.70),
    c(0.10, 0.20, 0.50, 0.51, 0.52, 0.53),
    c(0.10, 0.20, 0.30, 0.55, 0.40, 0.30),
    c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75),
    c(0.55, 0.50, 0.45, 0.40, 0.35, 0.30),
    c(0.10, 0.20, 0.50, 0.60, 0.65, 0.70),
    c(0.30, 0.45, 0.46, 0.47, 0.48, 0.49),
    c(0.10, 0.15, 0.20, 0.50, 0.51, 0.52)
  )
  pr <- rbind(
    c(0.05, 0.20, 0.35, 0.20, 0.20, 0.20),
    c(0.30, 0.40, 0.20, 0.10, 0.10, 0.10),
    c(0.40, 0.20, 0.20, 0.20, 0.20, 0.20),
    c(0.30, 0.40, 0.20, 0.21, 0.22, 0.23),
    c(0.30, 0.30, 0.30, 0.20, 0.20, 0.20),
    c(0.20, 0.20, 0.20, 0.20, 0.20, 0.20),
    c(0.20, 0.20, 0.20, 0.20, 0.20, 0.20),
    c(0.30, 0.40, 0.20, 0.20, 0.20, 0.20),
    c(0.30, 0.30, 0.30, 0.30, 0.30, 0.30),
    c(0.20, 0.25, 0.45, 0.20, 0.20, 0.20)
  )
  lapply(seq_along(tox), function(s) {
    t2 <- tox[[s]][2, ]
    t3 <- tox[[s]][3, ]
    t4 <- tox[[s]][4, ]
    if (type == "table") {
      list(p_grade = sweep(tox[[s]], 2L, colSums(tox[[s]]), "/"),
           p_resp_grade = rbind(1 - cr[s, ] - pr[s, ], pr[s, ], cr[s, ]))
    } else {
      g4 <- t4
      g3 <- (1 - t4) * t3
      g2 <- (1 - t4) * (1 - t3) * t2
      list(p_grade = rbind(1 - g4 - g3 - g2, g2, g3, g4),
           p_resp_grade = rbind((1 - cr[s, ]) * (1 - pr[s, ]),
                                (1 - cr[s, ]) * pr[s, ], cr[s, ]))
    }
  })
}

# gBOIN-ET rows of Table 1: percentage of early terminations and of
# recommendations of doses 1 to 6, and average number of patients at doses
# 1 to 6, one row per scenario.
takeda2022_table1 <- function() {
  list(
    rec = rbind(
      c(2, 1, 3, 22, 61, 11, 0),
      c(1, 4, 9, 57, 22, 6, 0),
      c(1, 11, 53, 28, 6, 1, 0),
      c(2, 5, 22, 52, 18, 3, 0),
      c(0, 4, 4, 15, 47, 16, 14),
      c(3, 72, 22, 3, 0, 0, 0),
      c(1, 64, 25, 7, 2, 0, 0),
      c(1, 3, 11, 51, 30, 4, 0),
      c(2, 15, 56, 24, 3, 0, 0),
      c(1, 2, 1, 13, 67, 16, 1)
    ),
    pts = rbind(
      c(3.4, 4.7, 7.4, 10.5, 3.2, 0.4),
      c(3.9, 6.6, 12.8, 4.1, 1.8, 0.5),
      c(7.9, 13.3, 6.6, 1.6, 0.2, 0.0),
      c(4.9, 8.9, 10.7, 4.0, 0.9, 0.1),
      c(3.8, 4.9, 6.4, 9.4, 3.1, 2.4),
      c(17.7, 7.1, 1.4, 0.1, 0.0, 0.0),
      c(18.7, 7.0, 2.5, 0.8, 0.2, 0.0),
      c(4.4, 7.7, 11.1, 4.9, 1.3, 0.2),
      c(7.4, 13.3, 7.1, 1.9, 0.2, 0.0),
      c(3.5, 4.0, 6.4, 11.3, 4.2, 0.6)
    )
  )
}

# Boundaries printed in Sections 3.1 and 3.3.
reproduce_takeda2022_boundaries <- function() {
  sets <- list(
    list(item = "Section 3.1", design = "gBOIN-ET", phi = 0.47 / 1.5,
         delta = 0.95 / 1.5, pub = c(0.14, 0.37, 0.50)),
    list(item = "Section 3.1", design = "BOIN-ET", phi = 0.33, delta = 0.70,
         pub = c(0.15, 0.39, 0.56)),
    list(item = "Section 3.3", design = "gBOIN-ET", phi = 0.35 / 1.5,
         delta = 0.65 / 1.5, pub = c(0.10, 0.27, 0.34)),
    list(item = "Section 3.3", design = "BOIN-ET", phi = 0.20, delta = 0.50,
         pub = c(0.09, 0.23, 0.40))
  )
  do.call(rbind, lapply(sets, function(z) {
    bd <- boin_et_boundary(z$phi, z$delta)
    data.frame(
      item = z$item, design = z$design,
      quantity = c("lambda1", "lambda2", "eta1"), published = z$pub,
      reproduced = c(bd$lambda1, bd$lambda2, bd$eta1), tolerance = 0.005,
      judged = TRUE, stringsAsFactors = FALSE
    )
  }))
}

reproduce_takeda2022_table1 <- function(n_trials = 10000, seed = 2022,
                                        n_pub = 4000) {
  pub <- takeda2022_table1()
  factor_se <- sqrt(1 / n_pub + 1 / n_trials)
  runs <- list(
    list(label = "gBOIN-ET", type = "simulated", compat = TRUE,
         judged = TRUE),
    list(label = "gBOIN-ET (rules of the package)", type = "simulated",
         compat = FALSE, judged = FALSE),
    list(label = "gBOIN-ET (probabilities of Table S1)", type = "table",
         compat = TRUE, judged = FALSE)
  )
  rows <- list()
  for (run in runs) {
    sc <- takeda2022_scenarios(run$type)
    for (s in seq_along(sc)) {
      oc <- sim_gboin_et(
        target = 0.47 / 1.5, target_resp = 0.95 / 1.5,
        p_grade = sc[[s]]$p_grade, p_resp_grade = sc[[s]]$p_resp_grade,
        score = c(0, 0.5, 1, 1.5), score_resp = c(0, 1, 1.5), n_cohort = 10,
        cohort_size = 3, n_trials = n_trials, lambda1 = 0.14133,
        lambda2 = 0.37333, eta1 = 0.50, resp_cor = 0,
        boinet_compatible = run$compat, est_method = "fp_linear",
        min_obd_sample = 6, keep_trials = TRUE, seed = seed
      )
      rec <- c(oc$percent_no_obd, unname(oc$sel_percent))
      p <- (rec + pub$rec[s, ]) / 200
      rows[[length(rows) + 1L]] <- data.frame(
        item = paste0("Scenario ", s), design = run$label,
        quantity = c("No OD (%)", paste0("Selection (%) DL", 1:6)),
        published = pub$rec[s, ], reproduced = rec,
        tolerance = 400 * sqrt(p * (1 - p)) * factor_se + 0.5,
        judged = run$judged, stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- data.frame(
        item = paste0("Scenario ", s), design = run$label,
        quantity = paste0("Patients DL", 1:6),
        published = pub$pts[s, ], reproduced = unname(oc$n_pts_dose),
        tolerance = 4 * apply(oc$trials$n_pts, 2L, stats::sd) * factor_se +
          0.05,
        judged = run$judged, stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
takeda2022_apply_explained <- function(res, file,
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

run_reproduce_takeda2022 <- function(n_trials = 10000, seed = 2022,
                                     out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_gboin_et", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- rbind(reproduce_takeda2022_boundaries(),
               reproduce_takeda2022_table1(n_trials = n_trials, seed = seed))
  rownames(res) <- NULL
  within <- abs(res$published - res$reproduced) <= res$tolerance + 1e-12
  res$judgement <- ifelse(!res$judged, "INFO",
                          ifelse(within, "PASS", "FAIL"))
  res$judged <- NULL
  res <- takeda2022_apply_explained(res, "reproduce_Takeda2022.R")
  res$reproduced <- round(res$reproduced, 3)
  res$tolerance <- round(res$tolerance, 3)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Takeda2022.csv"),
                   row.names = FALSE)
  res
}
