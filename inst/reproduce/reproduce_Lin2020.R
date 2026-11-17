# Reproduction of Lin et al. (2020) and of the BOIN12 part of Zhou et al.
# (2022) with simFastBOIN
#
# Lin R, Zhou Y, Yan F, Li D, Yuan Y. BOIN12: Bayesian optimal interval phase
# I/II trial design for utility-based dose finding in immunotherapy and
# targeted therapies. JCO Precision Oncology 2020; 4: 1393-1402.
# doi:10.1200/PO.20.00257
#
# Zhou Y, Lin R, Lee JJ, Li D, Wang L, Li R, Yuan Y. TITE-BOIN12: A Bayesian
# phase I/II trial design to find the optimal biological dose with late-onset
# toxicity and efficacy. Statistics in Medicine 2022; 41: 1918-1931.
# doi:10.1002/sim.9337
#
# The deterministic quantities are compared exactly or to the printed digits:
#   - Table 3 of Lin et al.: the rank-based desirability score of the 166
#     outcomes of 0, 3, 6 and 9 patients (upper DLT limit 0.35, lower response
#     limit 0.25, utilities 100, 40, 60 and 0), read from the text of the
#     article into tests/testthat/fixtures/lin2020-table3.csv; PASS when the
#     score, or E, is the same;
#   - the escalation and de-escalation boundaries 0.276 and 0.419 of Figure 1
#     and the desirability 0.295 of an untried dose (Zhou et al., Section
#     2.2.4); PASS within half a unit of the last digit;
#   - the estimated utilities 40.0, 46.7, 53.3, 60.0 and 50.0 of the example of
#     Section 6 of the Supplement of Lin et al., which are observed means
#     (est_method = "observed"), and its OBD, dose 4;
#   - the desirabilities of the BOIN12 rows of Table S.3 of the Supplement of
#     Zhou et al., from the patient data of their Table 2. The value 0.199 of
#     dose 3 at the fifth analysis is not reproduced (0.191 from six patients
#     with two DLTs and three responses); the same row reports a response rate
#     of 0.50 for dose 4, whose three patients had no response, so that row
#     appears to hold a misprint; the cell is registered as EXPLAINED in
#     inst/validation/explained.csv.
#
# The simulation results of BOIN12 are compared within Monte Carlo error:
#   - Figure 3 (page 1399) for Scenarios 1-9 and Figures S3 and S4 of the
#     Supplement (pages 23 and 24) for Scenarios 10-16 of Table S3: (A)
#     percentage of correct selection of the OBD, (B) patients treated at
#     the OBD, (C) patients treated at overdoses (DLT probability above
#     0.35) and (D) risk of poor allocation, the percentage of trials with
#     fewer than 36 / 5 = 7.2 patients at the OBD. The values of Figure 3
#     are its printed labels, read by their positions and checked against
#     the heights of the bars and the rendered figure; those of Figures S3
#     and S4, which have no labels, are the plotted points of the BOIN12
#     line, read from the vector graphics of the PDF (their averages agree
#     with the printed averages, and Scenarios 1-9 agree with Figure 3 to
#     0.05). They are in inst/reproduce/lin2020-figure3.csv.
#     The article simulated 10,000 trials; the tolerance is three standard
#     errors of the difference with the simulation here, plus 0.05.
#     The trials follow the rules of the article where they differ from the
#     BOIN12 software: ties of desirability broken at random (note of Table
#     S2 of the Supplement; tie_break = "random") and the end of the trial
#     when no dose allowed for the next cohort is admissible (inferred from
#     the published results; no_candidate = "stop"). The same values with
#     the default rules of sim_boin12() (tie_break = "higher",
#     no_candidate = "nearest") are listed as INFO, with "(default rules)"
#     after the item.
#
# Settings of the article and of Section 7 of its Supplement: five doses,
# upper DLT limit 0.35, lower response limit 0.25, utilities 100, 40, 60 and
# 0, cohorts of 3, at most 36 patients, N* = 6, cutoffs 0.95 and 0.90,
# independent outcomes (correlation 0 of the latent normal variables), start
# at dose 1, no early stopping; the defaults of sim_boin12() except
# tie_break and no_candidate.
#
# The OBD of a scenario follows Table S4 of the Supplement: among the doses
# with p_T <= 0.35 and p_E >= 0.25, those whose utility
# 100 p_E (1 - p_T) + 40 (1 - p_T)(1 - p_E) + 60 p_T p_E is within 3 points
# of the highest (one dose in each of Scenarios 1-14). In Scenarios 15 and
# 16 no dose is the OBD, and the Supplement reports the percentage of early
# termination as (A) and the number of patients not enrolled as (B); (D) is
# then plotted as 0 and listed as INFO.
#
# Usage, from the package root (the simulations take about two minutes):
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Lin2020.R")
#   res <- run_reproduce_lin2020()
#
# The judgements are written to inst/reproduce/results/reproduce_Lin2020.csv.

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
apply_explained <- function(res, file, item, quantity,
                            path = file.path("inst", "validation",
                                             "explained.csv")) {
  if (is.null(res$reason)) res$reason <- ""
  if (!file.exists(path)) return(res)
  ex <- utils::read.csv(path, stringsAsFactors = FALSE)
  ex <- ex[ex$file == file, , drop = FALSE]
  hit <- match(paste(item, quantity, sep = "\r"),
               paste(ex$item, ex$quantity, sep = "\r"))
  use <- !is.na(hit) & res$judgement == "FAIL"
  res$judgement[use] <- "EXPLAINED"
  res$reason[use] <- ex$reason[hit[use]]
  res
}

reproduce_lin2020_table3 <- function(path = file.path("tests", "testthat",
                                                      "fixtures",
                                                      "lin2020-table3.csv")) {
  pub <- utils::read.csv(path, colClasses = c("integer", "integer",
                                              "integer", "character"))
  tab <- boin12_decision_table(target = 0.35, target_resp = 0.25,
                               n = c(0, 3, 6, 9))
  key_tab <- paste(tab$n, tab$n_tox, tab$n_resp)
  key_pub <- paste(pub$n, pub$n_tox, pub$n_resp)
  idx <- match(key_pub, key_tab)
  rep_score <- ifelse(tab$admissible[idx], as.character(tab$score[idx]), "E")
  data.frame(
    source = "Lin 2020 Table 3",
    item = paste0("n=", pub$n, " tox=", pub$n_tox, " eff=", pub$n_resp),
    published = pub$score,
    reproduced = rep_score,
    tolerance = NA_real_,
    judgement = ifelse(!is.na(idx) & rep_score == pub$score, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  )
}

reproduce_lin2020_numbers <- function() {
  tab <- boin12_decision_table(target = 0.35, target_resp = 0.25,
                               n = c(0, 3, 6, 9))
  des <- attr(tab, "design")
  pick <- function(n, y, x) {
    tab$desirability[tab$n == n & tab$n_tox == y & tab$n_resp == x]
  }
  ex <- boin12_select_obd(n_pts = c(3, 3, 6, 12, 6), n_tox = c(0, 1, 1, 3, 3),
                          n_resp = c(0, 1, 2, 6, 3), target = 0.35,
                          target_resp = 0.25, est_method = "observed")
  rows <- list(
    c("Lin 2020 Figure 1", "lambda_e", 0.276, des$lambda_e, 5e-4),
    c("Lin 2020 Figure 1", "lambda_d", 0.419, des$lambda_d, 5e-4),
    c("Zhou 2022 Section 2.2.4", "desirability of an untried dose", 0.295,
      pick(0, 0, 0), 5e-4),
    c("Lin 2020 SI Section 6", "utility dose 1", 40.0, ex$utility[1, 1], 0.05),
    c("Lin 2020 SI Section 6", "utility dose 2", 46.7, ex$utility[1, 2], 0.05),
    c("Lin 2020 SI Section 6", "utility dose 3", 53.3, ex$utility[1, 3], 0.05),
    c("Lin 2020 SI Section 6", "utility dose 4", 60.0, ex$utility[1, 4], 0.05),
    c("Lin 2020 SI Section 6", "utility dose 5", 50.0, ex$utility[1, 5], 0.05),
    c("Lin 2020 SI Section 6", "OBD", 4, ex$obd, 0),
    c("Zhou 2022 Table S.3", "analyses 1-6, dose 1 (3 patients, 0, 0)", 0.113,
      pick(3, 0, 0), 5e-4),
    c("Zhou 2022 Table S.3", "analyses 2-6, dose 2 (3 patients, 1, 1)", 0.156,
      pick(3, 1, 1), 5e-4),
    c("Zhou 2022 Table S.3", "analyses 3-4, dose 3 (3 patients, 1, 1)", 0.156,
      pick(3, 1, 1), 5e-4),
    c("Zhou 2022 Table S.3", "analyses 4-6, dose 4 (3 patients, 2, 0)", 0.022,
      pick(3, 2, 0), 5e-4),
    c("Zhou 2022 Table S.3", "analysis 5, dose 3 (6 patients, 2, 3)", 0.199,
      pick(6, 2, 3), 5e-4),
    c("Zhou 2022 Table S.3", "analysis 6, dose 3 (9 patients, 3, 5)", 0.208,
      pick(9, 3, 5), 5e-4)
  )
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("source", "item", "published", "reproduced", "tolerance")
  pub <- as.numeric(out$published)
  rep <- as.numeric(out$reproduced)
  tol <- as.numeric(out$tolerance)
  out$judgement <- ifelse(abs(pub - rep) <= tol + 1e-12, "PASS", "FAIL")
  out$reproduced <- format(round(rep, 4), trim = TRUE)
  out$tolerance <- tol
  out
}

# Scenarios 1-16 of Table S3 of the Supplement of Lin et al. (2020);
# Scenarios 1-9 are those of Figures 2 and 3 of the article.
lin2020_scenarios <- function() {
  p <- list(
    c(0.05, 0.20, 0.12, 0.35, 0.27, 0.36, 0.35, 0.37, 0.50, 0.38),
    c(0.05, 0.01, 0.12, 0.05, 0.27, 0.30, 0.35, 0.60, 0.50, 0.60),
    c(0.03, 0.05, 0.05, 0.10, 0.20, 0.50, 0.22, 0.68, 0.45, 0.70),
    c(0.05, 0.01, 0.10, 0.05, 0.15, 0.10, 0.20, 0.20, 0.27, 0.40),
    c(0.03, 0.10, 0.06, 0.20, 0.10, 0.40, 0.30, 0.45, 0.45, 0.50),
    c(0.02, 0.05, 0.05, 0.15, 0.10, 0.40, 0.20, 0.40, 0.30, 0.40),
    c(0.01, 0.45, 0.15, 0.45, 0.20, 0.45, 0.40, 0.65, 0.50, 0.80),
    c(0.15, 0.25, 0.20, 0.55, 0.25, 0.40, 0.35, 0.30, 0.45, 0.20),
    c(0.03, 0.45, 0.05, 0.30, 0.07, 0.25, 0.09, 0.20, 0.11, 0.10),
    c(0.01, 0.07, 0.08, 0.10, 0.10, 0.35, 0.20, 0.20, 0.25, 0.20),
    c(0.08, 0.10, 0.10, 0.20, 0.15, 0.70, 0.32, 0.70, 0.40, 0.75),
    c(0.05, 0.05, 0.07, 0.15, 0.10, 0.30, 0.20, 0.55, 0.35, 0.55),
    c(0.01, 0.10, 0.05, 0.15, 0.10, 0.25, 0.12, 0.50, 0.27, 0.50),
    c(0.02, 0.05, 0.05, 0.08, 0.07, 0.15, 0.10, 0.30, 0.15, 0.45),
    c(0.50, 0.30, 0.55, 0.35, 0.60, 0.45, 0.65, 0.66, 0.70, 0.70),
    c(0.25, 0.01, 0.55, 0.55, 0.67, 0.55, 0.87, 0.55, 0.95, 0.55)
  )
  # Each row of Table S3 lists the pairs (p_T, p_E) of doses 1-5.
  lapply(p, function(v) list(p_true = v[c(1, 3, 5, 7, 9)],
                             p_resp = v[c(2, 4, 6, 8, 10)]))
}

# The OBDs of a scenario as in Table S4 of the Supplement (none in
# Scenarios 15 and 16). The small allowances guard the comparisons with the
# limits and with the highest utility against rounding.
lin2020_true_obd <- function(p_true, p_resp) {
  u <- 100 * p_resp * (1 - p_true) + 40 * (1 - p_true) * (1 - p_resp) +
    60 * p_true * p_resp
  ok <- which(p_true <= 0.35 + 1e-12 & p_resp >= 0.25 - 1e-12)
  if (length(ok) == 0L) return(integer(0))
  ok[u[ok] >= max(u[ok]) - 3 - 1e-9]
}

# Metrics of Figure 3 and Figures S3 and S4 from the trials and OBDs of one
# simulation, with the standard deviation of the per-trial values for the
# tolerances.
lin2020_metrics <- function(oc, p_true, p_resp, n_max = 36) {
  obd <- lin2020_true_obd(p_true, p_resp)
  over <- which(p_true > 0.35)
  n_pts <- oc$trials$n_pts
  n_over <- if (length(over)) rowSums(n_pts[, over, drop = FALSE]) else
    rep(0, nrow(n_pts))
  if (length(obd)) {
    n_obd <- rowSums(n_pts[, obd, drop = FALSE])
    correct <- 100 * mean(oc$obd %in% obd)
    poor <- 100 * mean(n_obd < n_max / 5)
  } else {
    # No OBD: early termination and patients not enrolled.
    n_obd <- n_max - rowSums(n_pts)
    correct <- 100 * mean(oc$trials$early_stop)
    poor <- NA_real_
  }
  list(
    obd = obd,
    value = c(correct_obd_percent = correct, n_pts_obd = mean(n_obd),
              n_pts_overdose = mean(n_over), poor_allocation_percent = poor),
    sd = c(correct_obd_percent = NA_real_, n_pts_obd = stats::sd(n_obd),
           n_pts_overdose = stats::sd(n_over),
           poor_allocation_percent = NA_real_)
  )
}

# With label = NULL the values are judged; otherwise every row is INFO, the
# label is added to the item and the reason tells whether the value is
# within the tolerance.
reproduce_lin2020_figure3 <- function(
    n_trials = 10000, seed = 2020,
    path = file.path("inst", "reproduce", "lin2020-figure3.csv"),
    tie_break = "random", no_candidate = "stop", label = NULL) {
  pub <- utils::read.csv(path, stringsAsFactors = FALSE)
  # Figure 3 for Scenarios 1-9, Figures S3 and S4 for Scenarios 10-16.
  pub <- pub[(pub$scenario <= 9) == grepl("^Lin 2020 Figure 3", pub$source),
             , drop = FALSE]
  scen <- lin2020_scenarios()
  n_pub <- 10000
  metric_names <- c("correct_obd_percent", "n_pts_obd", "n_pts_overdose",
                    "poor_allocation_percent")
  rows <- list()
  for (s in seq_along(scen)) {
    oc <- sim_boin12(target = 0.35, target_resp = 0.25,
                     p_true = scen[[s]]$p_true, p_resp = scen[[s]]$p_resp,
                     n_cohort = 12, cohort_size = 3, n_trials = n_trials,
                     tie_break = tie_break, no_candidate = no_candidate,
                     keep_trials = TRUE, seed = seed)
    m <- lin2020_metrics(oc, scen[[s]]$p_true, scen[[s]]$p_resp)
    for (nm in metric_names) {
      hit <- pub$scenario == s & pub$metric == nm
      stopifnot(sum(hit) == 1L)
      p <- pub$value[hit]
      rep <- unname(m$value[nm])
      if (nm %in% c("correct_obd_percent", "poor_allocation_percent")) {
        q <- max(p, rep, 0.5) / 100
        se <- 100 * sqrt(q * (1 - q) * (1 / n_pub + 1 / n_trials))
      } else {
        se <- unname(m$sd[nm]) * sqrt(1 / n_pub + 1 / n_trials)
      }
      tol <- 3 * se + 0.05
      if (is.na(rep)) {
        judgement <- "INFO"
        reason <- "no dose is the OBD; the figure plots 0"
      } else {
        judgement <- if (abs(rep - p) <= tol) "PASS" else "FAIL"
        reason <- ""
      }
      item <- paste0("Scenario ", s, ", ", nm)
      if (!is.null(label)) {
        item <- paste0(item, " (", label, ")")
        if (judgement != "INFO") {
          reason <- if (judgement == "PASS") "within the tolerance" else
            "outside the tolerance"
          reason <- paste0("tie_break = \"", tie_break, "\", no_candidate = \"",
                           no_candidate, "\"; ", reason)
        }
        judgement <- "INFO"
      }
      rows[[length(rows) + 1L]] <- data.frame(
        source = sub(" \\(p\\. [0-9]+\\)$", "", pub$source[hit]),
        item = item,
        published = format(p),
        reproduced = if (is.na(rep)) NA_character_ else
          format(round(rep, 2), nsmall = 2),
        tolerance = if (is.na(rep)) NA_real_ else round(tol, 2),
        judgement = judgement,
        reason = reason,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

run_reproduce_lin2020 <- function(out_dir = "inst/reproduce/results",
                                  n_trials = 10000) {
  # The functions of the package under development must be loaded.
  if (!exists("boin12_decision_table", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- rbind(reproduce_lin2020_table3(), reproduce_lin2020_numbers())
  res$reason <- ""
  res <- rbind(res, reproduce_lin2020_figure3(n_trials = n_trials))
  res <- rbind(res, reproduce_lin2020_figure3(
    n_trials = n_trials, tie_break = "higher", no_candidate = "nearest",
    label = "default rules"))
  rownames(res) <- NULL
  res <- apply_explained(res, "reproduce_Lin2020.R", res$item, res$source)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Lin2020.csv"),
                   row.names = FALSE)
  res
}
