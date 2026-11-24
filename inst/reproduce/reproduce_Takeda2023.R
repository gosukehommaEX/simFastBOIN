# Reproduction of Takeda, Yamaguchi, Taguri and Morita (2023) with simFastBOIN
#
# Takeda K, Yamaguchi Y, Taguri M, Morita S. TITE-gBOIN-ET: Time-to-event
# generalized Bayesian optimal interval design to accelerate dose-finding
# accounting for ordinal graded efficacy and toxicity outcomes. Biometrical
# Journal 2023; 65(7): 2200265. doi:10.1002/bimj.202200265
#
# The article prints the boundaries of its simulation study (Section 3.1),
# 0.14, 0.37 and 0.46, which boin_et_boundary() must give to the two printed
# decimals for the targets 0.47 / 1.5 and 1.75 / 3. (The article rounds the
# targets to 0.313 and 0.583, and its program searches a grid rounded to two
# decimals; with the rounded targets the grid of boin_et_boundary(), which
# is not rounded, gives 0.4698 for the third boundary. The simulations below
# take the three printed boundaries through lambda1, lambda2 and eta1.)
#
# Table 1 gives, for the ten scenarios of Table S1 of the supporting
# information, the percentage of trials selecting the optimal dose (PCS),
# the percentage of patients treated at it (PCA), the percentages of trials
# selecting, and of patients treated at, a higher dose (POS, POA), the
# percentage of trials without a selected dose (PET) and the average
# duration in days, under TITE-gBOIN-ET with the equivalent toxicity score
# (0 for grades 0 and 1, 0.5, 1 and 1.5 for grades 2, 3 and 4), the
# equivalent efficacy score (0, 0.25, 1 and 3 for progressive disease,
# stable disease, partial response and complete response), targets of 0.313
# and 0.583, 36 patients in cohorts of 3, a toxicity window of 30 days, an
# efficacy window of 60 days, one patient per 10 days on average and 10,000
# simulated trials. Table S2 gives the percentage of selections of each dose
# and the average number of patients at each dose. Tables S3 and S4 give the
# same for one patient per 5 days, and Tables S5 and S6 for an efficacy
# window of 90 days.
#
# The R program of the authors (supporting information, folder prog) was
# read and written again in C++ with R's random numbers
# (02_Rpackage/Validation/tite_gboin_et of the development files, not part
# of the package). That transcription reproduces all 4,560 summaries that
# the supporting information holds for the program (data/simulation1_*.csv,
# both designs, ten scenarios, three efficacy windows, two accrual rates)
# exactly. The program differs from the article, from
# boinet::tite.gboinet() and from the package in two points.
#
#   - sample(). When the current dose is neither left nor kept by the first
#     three rules, the program takes the dose with the highest mean efficacy
#     score among the current dose and its neighbors through
#     sample(doses, 1). When one dose d holds the highest score, sample(d, 1)
#     draws from 1, ..., d, so that the next cohort can be treated at any
#     dose up to d. 'boinet' takes a single candidate as it is, and so does
#     the package. This point moves the published tables substantially: in
#     Scenario 1, for example, the percentage of correct selections is 64.3
#     as published and about 75 without it.
#   - The elimination of doses. During the trial the program computes the
#     posterior probabilities of the elimination rule with the effective
#     numbers of patients (the patients assessed plus the fractions of the
#     window followed by the others); 'boinet' and the package use the
#     numbers of patients treated. This point moves mainly the percentage of
#     trials without a selected dose.
#
# The rows judged against the published tables use the rules of 'boinet'
# (boinet_compatible = TRUE), independent toxicity and efficacy
# (resp_cor = 0), the efficacy cutoff 0.95 of the article and the printed
# boundaries. A cell is judged PASS when it lies within four standard errors
# of the difference between two simulations of 10,000 and of n_trials
# trials, plus half a unit of the last printed digit, and FAIL otherwise.
# The cells that the two points above move out of this allowance are
# registered in inst/validation/explained.csv and become EXPLAINED.
#
# The same simulations are also judged against the transcription of the
# program with the two points removed (a single candidate dose is taken as
# it is, and the elimination uses the numbers of patients), run for 10,000
# trials: rows "TITE-gBOIN-ET against the corrected program", whose column
# published holds the value of that run. These values and the published
# cells are in inst/reproduce/takeda2023-tables.csv.
#
# The rows with the default rules of the package (the changes from 'boinet'
# described in boin_et_simulate(), an efficacy cutoff of 0.99 and a copula
# correlation of 0.2) are INFO.
#
# The rules of the package were checked trial by trial against an
# implementation in Python and, in distribution, against
# boinet::tite.gboinet() (inst/validation/compare-with-tite-gboinet.R).
#
# Usage, from the package root:
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Takeda2023.R")
#   res <- run_reproduce_takeda2023()
#
# The judgements are written to
# inst/reproduce/results/reproduce_Takeda2023.csv.

# Table S1 (data/true_probbaility.csv of the supporting information):
# probabilities of grades 0 or 1, 2, 3 and 4 and of progressive disease,
# stable disease, partial response and complete response at the six doses of
# each scenario, and the optimal dose of each scenario, the dose with the
# highest utility in that file.
takeda2023_scenarios <- function() {
  tox <- list(
    rbind(c(0.94, 0.87, 0.79, 0.68, 0.62, 0.50),
          c(0.05, 0.10, 0.15, 0.20, 0.20, 0.20),
          c(0.01, 0.03, 0.05, 0.10, 0.15, 0.25),
          c(0.00, 0.00, 0.01, 0.02, 0.03, 0.05)),
    rbind(c(0.93, 0.88, 0.77, 0.65, 0.50, 0.30),
          c(0.03, 0.05, 0.10, 0.15, 0.20, 0.30),
          c(0.03, 0.05, 0.10, 0.15, 0.20, 0.20),
          c(0.01, 0.02, 0.03, 0.05, 0.10, 0.20)),
    rbind(c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19),
          c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11),
          c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11),
          c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)),
    rbind(c(0.82, 0.65, 0.41, 0.42, 0.34, 0.26),
          c(0.10, 0.20, 0.34, 0.28, 0.31, 0.34),
          c(0.05, 0.10, 0.15, 0.18, 0.21, 0.24),
          c(0.03, 0.05, 0.10, 0.12, 0.14, 0.16)),
    rbind(c(0.88, 0.64, 0.52, 0.35, 0.17, 0.01),
          c(0.04, 0.12, 0.16, 0.22, 0.28, 0.24),
          c(0.06, 0.17, 0.22, 0.23, 0.25, 0.35),
          c(0.02, 0.07, 0.10, 0.20, 0.30, 0.40)),
    rbind(c(0.50, 0.35, 0.26, 0.17, 0.08, 0.01),
          c(0.40, 0.45, 0.50, 0.55, 0.60, 0.65),
          c(0.05, 0.10, 0.12, 0.14, 0.16, 0.18),
          c(0.05, 0.10, 0.12, 0.14, 0.16, 0.16)),
    rbind(c(0.50, 0.38, 0.29, 0.19, 0.13, 0.08),
          c(0.20, 0.24, 0.21, 0.16, 0.11, 0.07),
          c(0.25, 0.12, 0.12, 0.10, 0.08, 0.05),
          c(0.05, 0.26, 0.38, 0.55, 0.68, 0.80)),
    rbind(c(0.81, 0.72, 0.60, 0.45, 0.30, 0.15),
          c(0.15, 0.20, 0.25, 0.30, 0.35, 0.40),
          c(0.03, 0.05, 0.10, 0.15, 0.20, 0.25),
          c(0.01, 0.03, 0.05, 0.10, 0.15, 0.20)),
    rbind(c(0.84, 0.52, 0.36, 0.14, 0.05, 0.05),
          c(0.15, 0.47, 0.63, 0.85, 0.75, 0.65),
          c(0.01, 0.01, 0.01, 0.01, 0.10, 0.15),
          c(0.00, 0.00, 0.00, 0.00, 0.10, 0.15)),
    rbind(c(0.94, 0.87, 0.79, 0.68, 0.62, 0.50),
          c(0.05, 0.10, 0.15, 0.20, 0.20, 0.20),
          c(0.01, 0.03, 0.05, 0.10, 0.15, 0.25),
          c(0.00, 0.00, 0.01, 0.02, 0.03, 0.05))
  )
  eff <- list(
    rbind(c(0.64, 0.52, 0.45, 0.35, 0.20, 0.05),
          c(0.30, 0.40, 0.40, 0.40, 0.40, 0.15),
          c(0.05, 0.05, 0.10, 0.15, 0.20, 0.35),
          c(0.01, 0.03, 0.05, 0.10, 0.20, 0.45)),
    rbind(c(0.44, 0.37, 0.20, 0.10, 0.05, 0.05),
          c(0.40, 0.40, 0.50, 0.30, 0.20, 0.15),
          c(0.15, 0.20, 0.25, 0.35, 0.25, 0.30),
          c(0.01, 0.03, 0.05, 0.25, 0.50, 0.50)),
    rbind(c(0.60, 0.50, 0.30, 0.15, 0.10, 0.05),
          c(0.30, 0.20, 0.20, 0.15, 0.15, 0.15),
          c(0.05, 0.20, 0.35, 0.20, 0.20, 0.20),
          c(0.05, 0.10, 0.15, 0.50, 0.55, 0.60)),
    rbind(c(0.30, 0.20, 0.05, 0.05, 0.05, 0.05),
          c(0.35, 0.30, 0.25, 0.20, 0.15, 0.10),
          c(0.30, 0.40, 0.20, 0.25, 0.30, 0.30),
          c(0.05, 0.10, 0.50, 0.50, 0.50, 0.55)),
    rbind(c(0.40, 0.20, 0.10, 0.08, 0.06, 0.04),
          c(0.20, 0.20, 0.20, 0.20, 0.20, 0.20),
          c(0.30, 0.40, 0.20, 0.21, 0.22, 0.23),
          c(0.10, 0.20, 0.50, 0.51, 0.52, 0.53)),
    rbind(c(0.10, 0.05, 0.05, 0.05, 0.05, 0.05),
          c(0.20, 0.22, 0.21, 0.20, 0.19, 0.18),
          c(0.60, 0.24, 0.25, 0.26, 0.27, 0.28),
          c(0.10, 0.49, 0.49, 0.49, 0.49, 0.49)),
    rbind(c(0.15, 0.15, 0.10, 0.05, 0.05, 0.05),
          c(0.10, 0.10, 0.10, 0.10, 0.05, 0.02),
          c(0.25, 0.20, 0.20, 0.20, 0.20, 0.20),
          c(0.50, 0.55, 0.60, 0.65, 0.70, 0.73)),
    rbind(c(0.05, 0.10, 0.05, 0.10, 0.15, 0.25),
          c(0.30, 0.15, 0.25, 0.25, 0.30, 0.30),
          c(0.45, 0.20, 0.40, 0.40, 0.35, 0.30),
          c(0.20, 0.55, 0.30, 0.25, 0.20, 0.15)),
    rbind(c(0.70, 0.65, 0.60, 0.05, 0.05, 0.05),
          c(0.20, 0.20, 0.20, 0.45, 0.40, 0.35),
          c(0.00, 0.00, 0.00, 0.00, 0.05, 0.10),
          c(0.10, 0.15, 0.20, 0.50, 0.50, 0.50)),
    rbind(c(0.52, 0.45, 0.35, 0.20, 0.05, 0.05),
          c(0.40, 0.40, 0.40, 0.40, 0.15, 0.05),
          c(0.05, 0.10, 0.15, 0.20, 0.35, 0.10),
          c(0.03, 0.05, 0.10, 0.20, 0.45, 0.80))
  )
  od <- c(6, 5, 4, 3, 3, 2, 1, 2, 4, 6)
  lapply(seq_along(tox), function(s) {
    list(p_grade = tox[[s]], p_resp_grade = eff[[s]], od = od[s])
  })
}

# Boundaries printed in Section 3.1.
reproduce_takeda2023_boundaries <- function() {
  bd <- boin_et_boundary(0.47 / 1.5, 1.75 / 3)
  data.frame(
    item = "Section 3.1", design = "TITE-gBOIN-ET",
    quantity = c("lambda1", "lambda2", "eta1"),
    published = c(0.14, 0.37, 0.46),
    reproduced = c(bd$lambda1, bd$lambda2, bd$eta1), tolerance = 0.005,
    judged = TRUE, stringsAsFactors = FALSE
  )
}

# The quantities of the tables from one simulation, with the standard
# deviations over the trials that the allowances of the averages need.
takeda2023_quantities <- function(oc, od) {
  n <- oc$trials$n_pts
  n_doses <- ncol(n)
  total <- rowSums(n)
  above <- seq_len(n_doses) > od
  over <- rowSums(n[, above, drop = FALSE])
  pca <- sum(n[, od]) / sum(total)
  poa <- sum(over) / sum(total)
  dose_names <- paste0("DL", seq_len(n_doses))
  value <- c(
    unname(oc$sel_percent[od]), 100 * pca,
    sum(oc$sel_percent[above]), 100 * poa, oc$percent_no_obd, oc$duration,
    oc$percent_no_obd, unname(oc$sel_percent), unname(oc$n_pts_dose)
  )
  # NA marks a percentage of trials, whose allowance comes from the
  # binomial variance.
  spread <- c(
    NA, 100 * stats::sd(n[, od] - pca * total) / mean(total),
    NA, 100 * stats::sd(over - poa * total) / mean(total), NA,
    stats::sd(oc$trials$duration),
    NA, rep(NA, n_doses), apply(n, 2L, stats::sd)
  )
  names(value) <- c("PCS (%)", "PCA (%)", "POS (%)", "POA (%)", "PET (%)",
                    "Duration", "No OD (%)",
                    paste("Selection (%)", dose_names),
                    paste("Patients", dose_names))
  list(value = value, spread = spread)
}

# Half a unit of the last digit of a number as written in the file.
takeda2023_half_unit <- function(text) {
  digits <- ifelse(grepl(".", text, fixed = TRUE),
                   nchar(sub("^[^.]*\\.", "", text)), 0L)
  0.5 * 10^(-digits)
}

reproduce_takeda2023_tables <- function(n_trials = 10000, seed = 2023,
                                        n_pub = 10000,
                                        tables = file.path(
                                          "inst", "reproduce",
                                          "takeda2023-tables.csv"
                                        )) {
  pub <- utils::read.csv(tables, stringsAsFactors = FALSE,
                         colClasses = c("character", "numeric", "numeric",
                                        "integer", "character", "character",
                                        "character"))
  sc <- takeda2023_scenarios()
  factor_se <- sqrt(1 / n_pub + 1 / n_trials)
  settings <- unique(pub[, c("tau_e", "accrual")])
  runs <- list(
    list(label = "TITE-gBOIN-ET", compat = TRUE, cor = 0, cutoff = 0.95,
         judged = TRUE, all_settings = TRUE),
    list(label = "TITE-gBOIN-ET (rules of the package)", compat = FALSE,
         cor = 0.2, cutoff = 0.99, judged = FALSE, all_settings = FALSE)
  )
  rows <- list()
  for (run in runs) {
    for (k in seq_len(nrow(settings))) {
      if (!run$all_settings && k > 1L) next
      tau_e <- settings$tau_e[k]
      accrual <- settings$accrual[k]
      for (s in seq_along(sc)) {
        oc <- sim_tite_gboin_et(
          target = 0.313, target_resp = 0.583, p_grade = sc[[s]]$p_grade,
          p_resp_grade = sc[[s]]$p_resp_grade, score = c(0, 0.5, 1, 1.5),
          score_resp = c(0, 0.25, 1, 3), n_cohort = 12, cohort_size = 3,
          window = 30, accrual_rate = 1 / accrual, resp_window = tau_e,
          resp_cor = run$cor, n_trials = n_trials, lambda1 = 0.14,
          lambda2 = 0.37, eta1 = 0.46, cutoff_eli_resp = run$cutoff,
          boinet_compatible = run$compat, keep_trials = TRUE, seed = seed
        )
        q <- takeda2023_quantities(oc, sc[[s]]$od)
        cells <- pub[pub$tau_e == tau_e & pub$accrual == accrual &
                       pub$scenario == s, , drop = FALSE]
        at <- match(cells$quantity, names(q$value))
        reproduced <- unname(q$value[at])
        spread <- q$spread[at]
        allowance <- function(reference, half) {
          p <- (reference + reproduced) / 200
          ifelse(is.na(spread),
                 400 * sqrt(pmax(p * (1 - p), 0)) * factor_se,
                 4 * spread * factor_se) + half
        }
        published <- as.numeric(cells$published)
        rows[[length(rows) + 1L]] <- data.frame(
          item = paste0(cells$table, ", Scenario ", s), design = run$label,
          quantity = cells$quantity, published = published,
          reproduced = reproduced,
          tolerance = allowance(published,
                                takeda2023_half_unit(cells$published)),
          judged = run$judged, stringsAsFactors = FALSE
        )
        if (run$judged) {
          corrected <- as.numeric(cells$corrected)
          rows[[length(rows) + 1L]] <- data.frame(
            item = paste0(cells$table, ", Scenario ", s),
            design = "TITE-gBOIN-ET against the corrected program",
            quantity = cells$quantity, published = corrected,
            reproduced = reproduced,
            tolerance = allowance(corrected,
                                  takeda2023_half_unit(cells$corrected)),
            judged = TRUE, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
takeda2023_apply_explained <- function(res, file,
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

run_reproduce_takeda2023 <- function(n_trials = 10000, seed = 2023,
                                     out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_tite_gboin_et", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- rbind(reproduce_takeda2023_boundaries(),
               reproduce_takeda2023_tables(n_trials = n_trials, seed = seed))
  rownames(res) <- NULL
  within <- abs(res$published - res$reproduced) <= res$tolerance + 1e-12
  res$judgement <- ifelse(!res$judged, "INFO",
                          ifelse(within, "PASS", "FAIL"))
  res$judged <- NULL
  res <- takeda2023_apply_explained(res, "reproduce_Takeda2023.R")
  res$reproduced <- round(res$reproduced, 3)
  res$tolerance <- round(res$tolerance, 3)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Takeda2023.csv"),
                   row.names = FALSE)
  res
}
