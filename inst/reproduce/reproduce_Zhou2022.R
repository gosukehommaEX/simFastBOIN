# Reproduction of Zhou et al. (2022), TITE-BOIN12, with simFastBOIN
#
# Zhou Y, Lin R, Lee JJ, Li D, Wang L, Li R, Yuan Y. TITE-BOIN12: A Bayesian
# phase I/II trial design to find the optimal biological dose with late-onset
# toxicity and efficacy. Statistics in Medicine 2022; 41: 1918-1931.
# doi:10.1002/sim.9337
#
# Compared:
#   - Table S.3 of the Supplement, TITE-BOIN12 rows: the desirability and
#     the estimates (pi_T, pi_E) of every dose tried at the six analyses of
#     the example, from the patient data of Table 2 (TITE-BOIN12 columns).
#     The quantities of the table follow from replacing a pending outcome by
#     pi_hat itself (pending_estimate = "marginal"), which is judged: PASS
#     within half a unit of the printed digit. The values of Equation (11)
#     (pending_estimate = "conditional", the default) are listed as INFO.
#     Analysis 5, dose 3 cannot be obtained from Table 2: the table reports
#     pi_E = 0.55, but at day 315 two responses are observed and two patients
#     have completed the window without response at dose 3, so that pi_E is
#     at most 0.50. Analysis 6 reports pi_E = 0.56 for dose 4, whose three
#     patients had no response (pi_E = 0). Both are registered in
#     inst/validation/explained.csv.
#   - Figure 3 (page 1928): (A) percentage of correct selection of the OBD,
#     (B) patients treated at the OBD, (C) percentage of selecting a toxic
#     dose (DLT probability above 0.35) as the OBD and (D) patients treated
#     at toxic doses, for BOIN12 (sim_boin12()) and TITE-BOIN12AL
#     (sim_tite_boin12(), default settings), Scenarios 1-9 of Table S.4. The
#     values were read from the text of the figure by their positions and
#     checked against the rendered figure; they are in
#     inst/reproduce/zhou2022-figure3.csv. The article simulated 2000 trials;
#     the tolerance is three standard errors of the difference with the
#     simulation here, plus half a unit of the printed digit.
#     The article does not describe the accrual process, the distribution of
#     the times to event or how the duration is measured, so (E) the trial
#     duration is listed as INFO, as are the results of
#     pending_estimate = "marginal" and, for the duration of BOIN12,
#     TITE-BOIN12 waiting for every outcome (min_completed_share = 1).
#
# Settings of Section 3.1: five doses, upper DLT limit 0.35, lower response
# limit 0.25, utilities 100, 40, 60 and 0, cohorts of 3, at most 36 patients,
# three patients a month, windows of 45 and 60 days (1.5 and 2 months),
# cutoffs 0.95 and 0.90, independent outcomes (c = 0 of the Gumbel model),
# start at dose 1, no early stopping. Here the times between arrivals are
# uniform (accrual = "uniform", as in the code of Lin and Yuan, 2020) and the
# times to event are uniform over the windows (event_time = "uniform", the
# assumption of the approximated likelihood).
#
# The OBD of a scenario is the dose with the highest utility,
# 100 p_E (1 - p_T) + 40 (1 - p_T)(1 - p_E) + 60 p_T p_E, among the doses with
# p_T <= 0.35 (Scenario 1: dose 2; Scenario 2: dose 4, as in the text).
#
# Usage, from the package root (about two minutes):
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Zhou2022.R")
#   res <- run_reproduce_zhou2022()
#
# The judgements are written to inst/reproduce/results/reproduce_Zhou2022.csv.

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

# Table 2 of Zhou et al. (2022), TITE-BOIN12 columns: dose, DLT, response,
# day of enrollment and times to DLT and to response (days).
zhou2022_patients <- function() {
  data.frame(
    id = 1:18,
    dose = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 3, 3, 3, 3, 3, 3),
    y_tox = c(0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1),
    y_resp = c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0),
    day = c(1, 11, 21, 72, 82, 92, 143, 153, 163, 214, 224, 234, 285, 295,
            305, 316, 326, 336),
    x_tox = c(NA, NA, NA, NA, 30, NA, NA, 25, NA, NA, 40, 35, 25, NA, NA,
              NA, NA, 30),
    x_resp = c(NA, NA, NA, 40, NA, NA, 40, NA, NA, NA, NA, NA, 30, 55, NA,
               50, 50, NA)
  )
}

zhou2022_interim <- function(day, pending_estimate) {
  pts <- zhou2022_patients()
  des <- boin12_check_design(0.35, 0.25)
  tite_boin12_interim_cpp(
    dose = as.integer(pts$dose), enter = pts$day,
    end_tox = pts$day + ifelse(pts$y_tox == 1, pts$x_tox, 45),
    end_resp = pts$day + ifelse(pts$y_resp == 1, pts$x_resp, 60),
    dlt = as.integer(pts$y_tox), resp = as.integer(pts$y_resp), time = day,
    n_doses = 5L, window = 45, resp_window = 60, psi00 = 40, psi11 = 60,
    alpha = 1, beta = 1, u_benchmark = des$u_benchmark / 100, target = 0.35,
    target_resp = 0.25, cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
    pending = if (pending_estimate == "marginal") 1L else 0L
  )
}

reproduce_zhou2022_table_s3 <- function() {
  # Published TITE-BOIN12 rows of Table S.3: analysis, day, dose,
  # desirability, pi_T, pi_E (doses tried only).
  pub <- data.frame(
    analysis = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6),
    day = c(71, 142, 142, 213, 213, 213, 284, 284, 284, 284, 315, 315, 315,
            315, 396, 396, 396, 396),
    dose = c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
    desirability = c(0.113, 0.113, 0.211, 0.113, 0.156, 0.211, 0.113, 0.156,
                     0.156, 0.022, 0.113, 0.156, 0.193, 0.022, 0.113, 0.156,
                     0.208, 0.022),
    pi_tox = c(0, 0, 0.33, 0, 0.33, 0.33, 0, 0.33, 0.33, 0.67, 0, 0.33, 0.37,
               0.67, 0, 0.33, 0.33, 0.67),
    pi_resp = c(0, 0, 0.35, 0, 0.33, 0.35, 0, 0.33, 0.33, 0, 0, 0.33, 0.55,
                0, 0, 0.33, 0.56, 0.56)
  )
  rows <- list()
  for (method in c("marginal", "conditional")) {
    for (i in seq_len(nrow(pub))) {
      est <- zhou2022_interim(pub$day[i], method)
      j <- pub$dose[i]
      item <- paste0("analysis ", pub$analysis[i], " (day ", pub$day[i],
                     "), dose ", j)
      for (q in c("desirability", "pi_tox", "pi_resp")) {
        rep <- switch(q, desirability = est$desirability[j],
                      pi_tox = est$pi_tox[j], pi_resp = est$pi_resp[j])
        tol <- if (q == "desirability") 5e-4 else 5e-3
        rows[[length(rows) + 1L]] <- data.frame(
          source = if (method == "marginal") "Zhou 2022 Table S.3" else
            "Zhou 2022 Table S.3 (Equation 11)",
          item = item,
          quantity = q,
          published = pub[[q]][i],
          reproduced = round(rep, 4),
          tolerance = tol,
          judgement = if (method == "conditional") "INFO" else
            if (abs(rep - pub[[q]][i]) <= tol + 1e-12) "PASS" else "FAIL",
          reason = if (method == "conditional")
            "pending_estimate = \"conditional\" (Equation 11), the default" else
            "",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

zhou2022_scenarios <- function() {
  list(
    list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
         p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
    list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
         p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60)),
    list(p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
         p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70)),
    list(p_true = c(0.01, 0.06, 0.13, 0.15, 0.20),
         p_resp = c(0.02, 0.15, 0.18, 0.25, 0.40)),
    list(p_true = c(0.03, 0.06, 0.10, 0.30, 0.45),
         p_resp = c(0.10, 0.20, 0.40, 0.45, 0.50)),
    list(p_true = c(0.02, 0.05, 0.10, 0.20, 0.30),
         p_resp = c(0.05, 0.15, 0.40, 0.40, 0.40)),
    list(p_true = c(0.01, 0.15, 0.20, 0.40, 0.50),
         p_resp = c(0.45, 0.45, 0.45, 0.65, 0.80)),
    list(p_true = c(0.15, 0.20, 0.25, 0.35, 0.45),
         p_resp = c(0.25, 0.55, 0.40, 0.30, 0.20)),
    list(p_true = c(0.15, 0.30, 0.35, 0.40, 0.45),
         p_resp = c(0.30, 0.30, 0.30, 0.30, 0.30))
  )
}

zhou2022_true_obd <- function(p_true, p_resp) {
  u <- 100 * p_resp * (1 - p_true) + 40 * (1 - p_true) * (1 - p_resp) +
    60 * p_true * p_resp
  ok <- which(p_true <= 0.35)
  ok[which.max(u[ok])]
}

# Metrics of Figure 3 from the trials and OBDs of one simulation, with the
# standard deviation of the per-trial values for the tolerances.
zhou2022_metrics <- function(oc, p_true, p_resp) {
  obd <- zhou2022_true_obd(p_true, p_resp)
  toxic <- which(p_true > 0.35)
  sel <- oc$obd
  n_pts <- oc$trials$n_pts
  n_obd <- n_pts[, obd]
  n_toxic <- if (length(toxic)) rowSums(n_pts[, toxic, drop = FALSE]) else
    rep(0, nrow(n_pts))
  dur <- if (is.null(oc$trials$duration)) rep(NA_real_, nrow(n_pts)) else
    oc$trials$duration
  list(
    value = c(correct_obd_percent = 100 * mean(sel %in% obd),
              n_pts_obd = mean(n_obd),
              toxic_obd_percent = 100 * mean(sel %in% toxic),
              n_pts_toxic = mean(n_toxic),
              duration_months = mean(dur)),
    sd = c(correct_obd_percent = NA_real_, n_pts_obd = stats::sd(n_obd),
           toxic_obd_percent = NA_real_, n_pts_toxic = stats::sd(n_toxic),
           duration_months = stats::sd(dur))
  )
}

reproduce_zhou2022_figure3 <- function(
    n_trials = 10000, seed = 2022,
    path = file.path("inst", "reproduce", "zhou2022-figure3.csv")) {
  pub <- utils::read.csv(path, stringsAsFactors = FALSE)
  scen <- zhou2022_scenarios()
  n_pub <- 2000
  common <- list(target = 0.35, target_resp = 0.25, n_cohort = 12,
                 cohort_size = 3, n_trials = n_trials, keep_trials = TRUE,
                 seed = seed)
  timing <- list(window = 1.5, resp_window = 2, accrual_rate = 3,
                 accrual = "uniform", event_time = "uniform")
  runs <- list(
    list(label = "BOIN12", design = "BOIN12", fun = sim_boin12,
         extra = list(), judged = c(TRUE, TRUE, TRUE, TRUE, FALSE),
         note = "sim_boin12() has no time; see the timed BOIN12 rows"),
    list(label = "BOIN12 (min_completed_share = 1)", design = "BOIN12",
         fun = sim_tite_boin12,
         extra = c(timing, list(min_completed_share = 1)),
         judged = rep(FALSE, 5),
         note = paste("TITE-BOIN12 waiting for every outcome, which is",
                      "BOIN12 with time")),
    list(label = "TITE-BOIN12AL", design = "TITE-BOIN12AL",
         fun = sim_tite_boin12, extra = timing,
         judged = c(TRUE, TRUE, TRUE, TRUE, FALSE), note = ""),
    list(label = "TITE-BOIN12AL (marginal)", design = "TITE-BOIN12AL",
         fun = sim_tite_boin12,
         extra = c(timing, list(pending_estimate = "marginal")),
         judged = rep(FALSE, 5),
         note = "pending_estimate = \"marginal\"")
  )
  metric_names <- c("correct_obd_percent", "n_pts_obd", "toxic_obd_percent",
                    "n_pts_toxic", "duration_months")
  rows <- list()
  for (r in runs) {
    for (s in seq_along(scen)) {
      oc <- do.call(r$fun, c(common, list(p_true = scen[[s]]$p_true,
                                          p_resp = scen[[s]]$p_resp),
                             r$extra))
      m <- zhou2022_metrics(oc, scen[[s]]$p_true, scen[[s]]$p_resp)
      for (k in seq_along(metric_names)) {
        nm <- metric_names[k]
        p <- pub$value[pub$design == r$design & pub$scenario == s &
                         pub$metric == nm]
        rep <- unname(m$value[nm])
        if (nm %in% c("correct_obd_percent", "toxic_obd_percent")) {
          q <- max(p, rep, 0.5) / 100
          se <- 100 * sqrt(q * (1 - q) * (1 / n_pub + 1 / n_trials))
        } else {
          se <- unname(m$sd[nm]) * sqrt(1 / n_pub + 1 / n_trials)
        }
        tol <- 3 * se + 0.05
        if (is.na(rep)) {
          judgement <- "INFO"
        } else if (r$judged[k]) {
          judgement <- if (abs(rep - p) <= tol) "PASS" else "FAIL"
        } else {
          judgement <- "INFO"
        }
        reason <- r$note
        if (k == 5L && r$label != "BOIN12") {
          reason <- paste0(if (nzchar(reason)) paste0(reason, "; ") else "",
                           "the accrual process and the measure of the ",
                           "duration are not described in the article")
        }
        rows[[length(rows) + 1L]] <- data.frame(
          source = "Zhou 2022 Figure 3",
          item = paste0("Scenario ", s, ", ", r$label),
          quantity = nm,
          published = p,
          reproduced = if (is.na(rep)) NA_real_ else round(rep, 2),
          tolerance = if (is.na(rep)) NA_real_ else round(tol, 2),
          judgement = judgement,
          reason = reason,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

run_reproduce_zhou2022 <- function(out_dir = "inst/reproduce/results",
                                   n_trials = 10000) {
  # The functions of the package under development must be loaded.
  if (!exists("tite_boin12_interim_cpp", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- rbind(reproduce_zhou2022_table_s3(),
               reproduce_zhou2022_figure3(n_trials = n_trials))
  rownames(res) <- NULL
  res <- apply_explained(res, "reproduce_Zhou2022.R", res$item, res$quantity)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Zhou2022.csv"),
                   row.names = FALSE)
  res
}
