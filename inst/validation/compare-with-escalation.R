# Comparison of sim_boin12() with the BOIN12 design of the 'escalation'
# package
#
# escalation::get_boin12() (version 0.2.3) gives the dose for the next cohort
# and escalation::select_boin12_obd() the OBD at the end; both were written to
# agree with the article and the BOIN12 software. The two implementations do
# not draw the same random numbers, so they are compared in distribution: the
# percentage of trials selecting each dose as the OBD (or none) and the
# average number of patients at each dose, for the nine scenarios of Figure 2
# of Zhou et al. (2022) (Table S.4 of their Supplement; Scenario 1 is
# Scenario 1 of Lin et al., 2020), with an upper DLT limit of 0.35, a lower
# response limit of 0.25, utilities 100, 40, 60 and 0, 12 cohorts of 3 and
# independent outcomes. simFastBOIN uses est_method = "tail_probability",
# the estimate of 'escalation'.
#
# The rules differ in three places, documented in ?boin12_simulate: in
# simFastBOIN a dose made inadmissible by its DLT data also makes the higher
# doses inadmissible, a de-escalation never goes to an inadmissible dose, and
# the trial stops when no dose is admissible; 'escalation' checks only the
# doses next to the current one and keeps de-escalating to the dose below.
# The OBD of 'escalation' is selected without the admissibility check. These
# points matter only when a dose is inadmissible, so a cell is judged PASS
# when the two values lie within four standard errors of their difference,
# and FAIL otherwise; the FAILs are to be read against these differences.
#
# Usage, from the package root (the 'escalation' simulations take a few
# minutes):
#
# FAILs whose reason is registered in inst/validation/explained.csv are
# reported as EXPLAINED.
#
#   devtools::load_all()
#   source("inst/validation/compare-with-escalation.R")
#   res <- run_compare_with_escalation()
#
# The judgements are written to
# inst/validation/results/compare-with-escalation.csv.

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
apply_explained <- function(res, file, item, quantity,
                            path = file.path("inst", "validation",
                                             "explained.csv")) {
  res$reason <- ""
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

escalation_scenarios <- function() {
  tox <- list(
    c(0.05, 0.12, 0.27, 0.35, 0.50),
    c(0.05, 0.12, 0.27, 0.35, 0.50),
    c(0.03, 0.05, 0.20, 0.22, 0.45),
    c(0.01, 0.06, 0.13, 0.15, 0.20),
    c(0.03, 0.06, 0.10, 0.30, 0.45),
    c(0.02, 0.05, 0.10, 0.20, 0.30),
    c(0.01, 0.15, 0.20, 0.40, 0.50),
    c(0.15, 0.20, 0.25, 0.35, 0.45),
    c(0.15, 0.30, 0.35, 0.40, 0.45)
  )
  eff <- list(
    c(0.20, 0.35, 0.36, 0.37, 0.38),
    c(0.01, 0.05, 0.30, 0.60, 0.60),
    c(0.05, 0.10, 0.50, 0.68, 0.70),
    c(0.02, 0.15, 0.18, 0.25, 0.40),
    c(0.10, 0.20, 0.40, 0.45, 0.50),
    c(0.05, 0.15, 0.40, 0.40, 0.40),
    c(0.45, 0.45, 0.45, 0.65, 0.80),
    c(0.25, 0.55, 0.40, 0.30, 0.20),
    c(0.30, 0.30, 0.30, 0.30, 0.30)
  )
  lapply(seq_along(tox), function(i) {
    list(name = paste("Scenario", i), p_true = tox[[i]], p_resp = eff[[i]])
  })
}

# Selection percentages (doses and none) and patients per dose of one
# scenario from 'escalation'.
escalation_boin12_oc <- function(p_true, p_resp, n_sims) {
  model <- escalation::select_boin12_obd(
    escalation::stop_at_n(
      escalation::get_boin12(num_doses = length(p_true), phi_t = 0.35,
                             phi_e = 0.25, u2 = 40, u3 = 60, n_star = 6),
      n = 36
    )
  )
  sims <- escalation::simulate_trials(model, num_sims = n_sims,
                                      true_prob_tox = p_true,
                                      true_prob_eff = p_resp)
  rec <- escalation::recommended_dose(sims)
  npd <- as.matrix(escalation::n_at_dose(sims))
  list(rec = rec, npd = npd)
}

run_compare_with_escalation <- function(n_sims = 2000, n_trials = 20000,
                                        seed = 2026,
                                        out_dir = "inst/validation/results") {
  # The functions of the package under development and 'escalation' must be
  # available.
  if (!exists("sim_boin12", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  if (!requireNamespace("escalation", quietly = TRUE)) {
    stop("install the 'escalation' package first", call. = FALSE)
  }
  scenarios <- escalation_scenarios()
  set.seed(seed)
  rows <- list()
  for (s in seq_along(scenarios)) {
    sc <- scenarios[[s]]
    J <- length(sc$p_true)
    ref <- escalation_boin12_oc(sc$p_true, sc$p_resp, n_sims)
    oc <- sim_boin12(target = 0.35, target_resp = 0.25, p_true = sc$p_true,
                     p_resp = sc$p_resp, n_cohort = 12, cohort_size = 3,
                     n_trials = n_trials, est_method = "tail_probability",
                     keep_trials = TRUE, seed = seed + s)
    sel_ref <- vapply(seq_len(J), function(j) mean(ref$rec %in% j),
                      numeric(1))
    sel_ref <- c(sel_ref, mean(is.na(ref$rec)))
    sel_new <- c(oc$sel_percent, oc$percent_no_obd) / 100
    se_sel <- sqrt(sel_ref * (1 - sel_ref) / n_sims +
                     sel_new * (1 - sel_new) / n_trials)
    pts_ref <- colMeans(ref$npd)
    pts_new <- colMeans(oc$trials$n_pts)
    se_pts <- sqrt(apply(ref$npd, 2, stats::var) / n_sims +
                     apply(oc$trials$n_pts, 2, stats::var) / n_trials)
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = sc$name,
      quantity = c(paste0("OBD selected (%) DL", seq_len(J)),
                   "No OBD (%)", paste0("Patients DL", seq_len(J))),
      escalation = c(100 * sel_ref, pts_ref),
      simFastBOIN = c(100 * sel_new, pts_new),
      tolerance = c(400 * se_sel, 4 * se_pts),
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  res$judgement <- ifelse(abs(res$escalation - res$simFastBOIN) <=
                            res$tolerance + 1e-12, "PASS", "FAIL")
  res <- apply_explained(res, "compare-with-escalation.R", res$scenario,
                         res$quantity)
  res$escalation <- round(res$escalation, 2)
  res$simFastBOIN <- round(res$simFastBOIN, 2)
  res$tolerance <- round(res$tolerance, 2)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "compare-with-escalation.csv"),
                   row.names = FALSE)
  res
}
