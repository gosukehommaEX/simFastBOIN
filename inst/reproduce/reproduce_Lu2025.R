# Reproduction of Lu et al. (2025) with simFastBOIN
#
# Lu M, Zhang J, Yuan Y, Lin R. Comb-BOIN12: A utility-based Bayesian optimal
# interval design for dose optimization in cancer drug-combination trials.
# Statistics in Biopharmaceutical Research 2025; 17: 266-276.
# doi:10.1080/19466315.2024.2370403
#
# The published values are in inst/reproduce: lu2025-tables.csv (Tables 1 and
# 2 and Supplementary Tables S2 and S3, from the text of the PDFs; the bold
# OBDCs read from the page images), lu2025-table3-4.csv (the Comb-BOIN12 rows
# of Tables 3 and 4), and lu2025-figure2.csv, lu2025-figureS1.csv and
# lu2025-figureS2.csv (the bars and curves read from the vector graphics of
# the PDFs in two independent ways, value_a and value_b, which agree within
# 0.001 for the bars and 0.0011 for the curves; value_b is used).
#
# The deterministic quantities are compared to the printed digits:
#   - the expected utilities x 100 of Tables 1 and 2 (wT = 0.11) and of
#     Supplementary Tables S2 and S3 (wT = 0.05 and 0.17), 792 values, from
#     the true utility of oc_scenarios() with psi00 = 100 wT / (1 + wT) and
#     psi11 = 100 / (1 + wT), on the scale of the article:
#     u = U (1 + wT) / 100 - wT; PASS within 0.05;
#   - the OBDCs of Tables 1 and 2 (bold), the combinations with the highest
#     utility among those with a DLT probability at most phi_T (true_obd of
#     oc_scenarios() with tox_limit = target);
#   - the tail probabilities of Supplementary Figure S1, 164 points, from the
#     desirability of boin12_decision_table() with the same psi00 and psi11
#     and benchmark_weight = w_u; PASS within 0.001 (half a point of the
#     plot is 0.0009).
#
# The simulation results of Comb-BOIN12 are compared within Monte Carlo error
# (5,000 trials in the article; the tolerance is three standard errors of the
# difference plus 0.05, the standard error of a selection percentage being
# taken at the larger of the published and the simulated percentage, and at
# no less than 0.5 percent):
#   - Tables 3 and 4: the percentages of patients treated at the OBDCs, at
#     the promising combinations (PDCs: DLT probability at most phi_T and
#     response probability at least 0.4) and at the overly toxic
#     combinations (DLT probability above phi_T). They do not depend on the
#     final selection. A percentage is the average over trials of the number
#     of patients as a percentage of the maximum sample size; in Scenario 9,
#     where no combination is acceptable, the places left untreated count as
#     the "patients allocated to the OBDC", so that the two percentages of the
#     article add up to 100;
#   - Supplementary Figure S2 (page 5), whose bars for wT = 0.11 are also
#     those of Comb-BOIN12 in Figure 2: the selection of the OBDCs (in
#     Scenario 9, of no combination) and of the overly toxic combinations and
#     the percentages of patients there, for wT = 0.11, 0.05 and 0.17, with
#     obd_bound = "contour".
#   - Figure 2 (page 274): the bars A-D of Comb-BOIN12, the same for both
#     pairs of windows and the same as those of Supplementary Figure S2 for
#     wT = 0.11, against the same simulated values, so that a value of
#     Comb-BOIN12 that fails is listed up to three times; the bars A-D of
#     TITE-Comb-BOIN12 from tite_comb_boin12_simulate() with the windows
#     (AT, AE) = (1, 2) and (2, 3) months, two patients a month, uniform
#     gaps between arrivals (accrual = "uniform"), uniform event times
#     (event_time = "uniform", as Section 2.2 of the article assumes), no
#     correlation, the accrual suspension of TITE-BOIN12
#     (min_completed_share = 0.5), the pending outcomes counted only in the
#     utility, as Sections 2.2 and 3.1 of the article state the rules
#     (pending_use = "utility"; with the default pending_use = "all" more of
#     the bars of panels C and D differ), and obd_bound = "contour".
#   The selection percentages of Tables 3 and 4 are listed as INFO, with
#   obd_bound = "rectangle" and "contour": the trials of Tables 3 and 4 and
#   of Figure 2 have the same allocation percentages (they are the same
#   simulated trials), but different selection percentages, and no selection
#   rule that was tried reproduces those of Tables 3 and 4
#   (Validation/comb_boin12/README.md of the development folder). The trial
#   durations of Figure 2 (panel E) are listed as INFO: those of
#   TITE-Comb-BOIN12 come from the trials above, those of Comb-BOIN12 from
#   tite_comb_boin12_simulate() with min_completed_share = 1 (every decision
#   waits for all outcomes); the article does not describe the accrual and
#   the waiting rules, and under none of the six combinations of accrual
#   (uniform, exponential, fixed) and event times (uniform, Weibull) that
#   were tried do the durations agree (Comb-BOIN12 about 0.5 to 1.5 months
#   longer, TITE-Comb-BOIN12 about 1.3 to 2.3 months shorter on average).
#
# Settings of Section 4.1 of the article: drugs A and B with the dose values
# (80, 160, 240, 320) and (75, 150, 225, 300) mg in Scenarios 1-15 and
# (100, 200) and (50, 100, 200) mg in Scenarios 16-19; phi_T 0.35, 0.30 and
# 0.33, phi_E 0.2 and at most 51, 45 and 42 patients in Scenarios 1-9, 10-15
# and 16-19; cohorts of 3 from (1, 1); wT = 0.11; w_u = 0.6; cutoffs 0.95
# (0.85 at (1, 1), extrasafe = TRUE, offset = 0.10) and 0.9; N* = 6; Cauchy
# priors of scale 2.5. The article does not state the correlation of the DLT
# and the response (0 here), nor how ties are broken and what happens when no
# candidate is admissible: the defaults of sim_comb_boin12(), tie_break =
# "higher" and no_candidate = "nearest", agree with the published allocation
# percentages, whereas tie_break = "random" does not.
#
# Usage, from the package root (the simulations take about three hours with
# devtools::load_all(), which compiles the C++ code without optimization;
# most of the time goes to the 570,000 fits of the logistic model, and a
# message is shown for every simulation):
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Lu2025.R")
#   res <- run_reproduce_lu2025()
#
# The judgements are written to inst/reproduce/results/reproduce_Lu2025.csv.

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
lu2025_apply_explained <- function(res, file,
                                   path = file.path("inst", "validation",
                                                    "explained.csv")) {
  if (!file.exists(path)) return(res)
  ex <- utils::read.csv(path, stringsAsFactors = FALSE)
  ex <- ex[ex$file == file, , drop = FALSE]
  hit <- match(paste(res$item, res$source, sep = "\r"),
               paste(ex$item, ex$quantity, sep = "\r"))
  use <- !is.na(hit) & res$judgement == "FAIL"
  res$judgement[use] <- "EXPLAINED"
  res$reason[use] <- ex$reason[hit[use]]
  res
}

lu2025_read <- function(name) {
  utils::read.csv(file.path("inst", "reproduce", name),
                  stringsAsFactors = FALSE)
}

lu2025_psi <- function(w_tox) {
  c(psi00 = 100 * w_tox / (1 + w_tox), psi11 = 100 / (1 + w_tox))
}

# Scenarios 1-19 with the settings of Section 4.1.
lu2025_scenarios <- function(tab = lu2025_read("lu2025-tables.csv")) {
  lapply(1:19, function(s) {
    d <- tab[tab$scenario == s, , drop = FALSE]
    n_a <- max(d$dose_a)
    n_b <- max(d$dose_b)
    cell <- function(v) {
      m <- matrix(NA_real_, n_a, n_b)
      m[cbind(d$dose_a, d$dose_b)] <- v
      m
    }
    list(
      scenario = s,
      p_true = cell(d$p_true),
      p_resp = cell(d$p_resp),
      target = if (s <= 9) 0.35 else if (s <= 15) 0.30 else 0.33,
      n_cohort = if (s <= 9) 17L else if (s <= 15) 15L else 14L,
      dose_a = if (s <= 15) c(80, 160, 240, 320) else c(100, 200),
      dose_b = if (s <= 15) c(75, 150, 225, 300) else c(50, 100, 200),
      cells = d
    )
  })
}

lu2025_truth <- function(sc, w_tox) {
  psi <- lu2025_psi(w_tox)
  oc_scenarios(list(list(p_true = sc$p_true, p_resp = sc$p_resp)),
               target = sc$target, psi00 = psi[["psi00"]],
               psi11 = psi[["psi11"]])
}

lu2025_label <- function(m) {
  if (!any(m)) return("none")
  w <- which(m, arr.ind = TRUE)
  w <- w[order(w[, 2L], w[, 1L]), , drop = FALSE]
  paste0("A", w[, 1L], "B", w[, 2L], collapse = ", ")
}

reproduce_lu2025_utilities <- function() {
  rows <- list()
  for (sc in lu2025_scenarios()) {
    for (w_tox in c(0.11, 0.05, 0.17)) {
      truth <- lu2025_truth(sc, w_tox)
      col <- c("0.11" = "utility_011", "0.05" = "utility_005",
               "0.17" = "utility_017")[[format(w_tox)]]
      src <- if (w_tox == 0.11) sc$cells$source[1] else sc$cells$source_si[1]
      key <- paste(truth$dose_a, truth$dose_b)
      idx <- match(paste(sc$cells$dose_a, sc$cells$dose_b), key)
      u <- 100 * (truth$utility[idx] / 100 * (1 + w_tox) - w_tox)
      pub <- sc$cells[[col]]
      rows[[length(rows) + 1L]] <- data.frame(
        source = sub(" \\(p\\. [0-9]+\\)$", "", src),
        item = paste0("Scenario ", sc$scenario, ", A", sc$cells$dose_a, "B",
                      sc$cells$dose_b, ", utility x 100, wT = ", w_tox),
        published = sprintf("%.1f", pub),
        reproduced = sprintf("%.2f", u),
        tolerance = 0.05,
        judgement = ifelse(abs(u - pub) <= 0.05 + 1e-9, "PASS", "FAIL"),
        reason = "",
        stringsAsFactors = FALSE
      )
      if (w_tox == 0.11) {
        obd <- matrix(FALSE, nrow(sc$p_true), ncol(sc$p_true))
        obd[cbind(truth$dose_a, truth$dose_b)] <- truth$true_obd
        pub_obd <- matrix(FALSE, nrow(sc$p_true), ncol(sc$p_true))
        pub_obd[cbind(sc$cells$dose_a, sc$cells$dose_b)] <- sc$cells$obdc
        rows[[length(rows) + 1L]] <- data.frame(
          source = sub(" \\(p\\. [0-9]+\\)$", "", src),
          item = paste0("Scenario ", sc$scenario, ", OBDC (bold)"),
          published = lu2025_label(pub_obd),
          reproduced = lu2025_label(obd),
          tolerance = NA_real_,
          judgement = if (identical(obd, pub_obd)) "PASS" else "FAIL",
          reason = "",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

reproduce_lu2025_figure_s1 <- function() {
  pub <- lu2025_read("lu2025-figureS1.csv")
  data <- list("dose 1" = c(6, 1, 0), "dose 2" = c(9, 3, 1))
  rep <- numeric(nrow(pub))
  for (i in seq_len(nrow(pub))) {
    w_u <- as.numeric(sub("^wu=([0-9.]+),.*$", "\\1", pub$setting[i]))
    d <- data[[sub("^.*, ", "", pub$setting[i])]]
    psi <- lu2025_psi(pub$w_tox[i])
    tab <- boin12_decision_table(target = 0.3, target_resp = 0.25,
                                 n = d[1L], psi00 = psi[["psi00"]],
                                 psi11 = psi[["psi11"]],
                                 benchmark_weight = w_u)
    rep[i] <- tab$desirability[tab$n_tox == d[2L] & tab$n_resp == d[3L]]
  }
  data.frame(
    source = "Lu 2025 SI Figure S1",
    item = paste0(pub$setting, ", wT = ", sprintf("%.2f", pub$w_tox)),
    published = sprintf("%.4f", pub$prob_b),
    reproduced = sprintf("%.4f", rep),
    tolerance = 0.001,
    judgement = ifelse(abs(rep - pub$prob_b) <= 0.001, "PASS", "FAIL"),
    reason = "",
    stringsAsFactors = FALSE
  )
}

# Metrics of Tables 3 and 4 and of Figure S2 from the trials and the
# selection of one simulation, with the standard deviations of the per-trial
# percentages of patients for the tolerances.
lu2025_metrics <- function(trials, sel, sc, obdc) {
  n_max <- 3 * sc$n_cohort
  n_trials <- dim(trials$n_pts)[1L]
  pdc <- sc$p_true <= sc$target + 1e-12 & sc$p_resp >= 0.4 - 1e-12
  over <- sc$p_true > sc$target + 1e-12
  n_cells <- matrix(trials$n_pts, n_trials)
  n_total <- rowSums(n_cells)
  chosen <- rep(NA_integer_, n_trials)
  ok <- !is.na(sel$obd_a)
  chosen[ok] <- sel$obd_a[ok] + nrow(sc$p_true) * (sel$obd_b[ok] - 1L)
  pick <- function(m) ok & as.vector(m)[ifelse(ok, chosen, 1L)]
  at <- function(m) {
    if (!any(m)) return(rep(0, n_trials))
    rowSums(n_cells[, which(as.vector(m)), drop = FALSE])
  }
  if (any(obdc)) {
    sel_obdc <- 100 * mean(pick(obdc))
    pts_obdc <- 100 * at(obdc) / n_max
  } else {
    sel_obdc <- 100 * mean(!ok)
    pts_obdc <- 100 * (n_max - n_total) / n_max
  }
  pts_pdc <- 100 * at(pdc) / n_max
  pts_over <- 100 * at(over) / n_max
  list(
    value = c(sel_obdc = sel_obdc, pts_obdc = mean(pts_obdc),
              sel_pdc = 100 * mean(pick(pdc)), pts_pdc = mean(pts_pdc),
              sel_overtox = 100 * mean(pick(over)),
              pts_overtox = mean(pts_over)),
    sd = c(sel_obdc = NA_real_, pts_obdc = stats::sd(pts_obdc),
           sel_pdc = NA_real_, pts_pdc = stats::sd(pts_pdc),
           sel_overtox = NA_real_, pts_overtox = stats::sd(pts_over))
  )
}

lu2025_row <- function(source, item, p, m, nm, n_trials, judge = TRUE,
                       reason = "") {
  rep <- unname(m$value[nm])
  n_pub <- 5000
  if (grepl("^sel_", nm)) {
    q <- max(p, rep, 0.5) / 100
    se <- 100 * sqrt(q * (1 - q) * (1 / n_pub + 1 / n_trials))
  } else {
    se <- unname(m$sd[nm]) * sqrt(1 / n_pub + 1 / n_trials)
  }
  tol <- 3 * se + 0.05
  within <- abs(rep - p) <= tol
  if (judge) {
    judgement <- if (within) "PASS" else "FAIL"
  } else {
    judgement <- "INFO"
    reason <- paste0(reason, if (within) "within" else "outside",
                     " the tolerance")
  }
  data.frame(source = source, item = item, published = format(p, nsmall = 1),
             reproduced = format(round(rep, 2), nsmall = 2),
             tolerance = round(tol, 2), judgement = judgement,
             reason = reason, stringsAsFactors = FALSE)
}

reproduce_lu2025_simulation <- function(n_trials = 5000, seed = 2025) {
  tab34 <- lu2025_read("lu2025-table3-4.csv")
  fig <- lu2025_read("lu2025-figureS2.csv")
  fig2 <- lu2025_read("lu2025-figure2.csv")
  panels <- c(A_sel_obdc = "sel_obdc", B_pts_obdc = "pts_obdc",
              C_sel_overtox = "sel_overtox", D_pts_overtox = "pts_overtox")
  rows <- list()
  for (sc in lu2025_scenarios()) {
    s <- sc$scenario
    for (w_tox in c(0.11, 0.05, 0.17)) {
      message("Comb-BOIN12, Scenario ", s, ", wT = ", w_tox, ": ",
              format(Sys.time(), "%H:%M:%S"))
      psi <- lu2025_psi(w_tox)
      trials <- comb_boin12_simulate(
        target = sc$target, target_resp = 0.2, p_true = sc$p_true,
        p_resp = sc$p_resp, n_cohort = sc$n_cohort, cohort_size = 3,
        n_trials = n_trials, extrasafe = TRUE, offset = 0.10,
        psi00 = psi[["psi00"]], psi11 = psi[["psi11"]],
        benchmark_weight = 0.6, seed = seed)
      truth <- lu2025_truth(sc, w_tox)
      obdc <- matrix(FALSE, nrow(sc$p_true), ncol(sc$p_true))
      obdc[cbind(truth$dose_a, truth$dose_b)] <- truth$true_obd
      bounds <- if (w_tox == 0.11) c("contour", "rectangle") else "contour"
      for (b in bounds) {
        sel <- comb_boin12_select_obd(
          trials$n_pts, trials$n_tox, trials$n_resp, target = sc$target,
          target_resp = 0.2, extrasafe = TRUE, offset = 0.10,
          psi00 = psi[["psi00"]], psi11 = psi[["psi11"]], obd_bound = b,
          dose_a = sc$dose_a, dose_b = sc$dose_b,
          early_stop = trials$early_stop)
        m <- lu2025_metrics(trials, sel, sc, obdc)
        if (w_tox == 0.11) {
          for (nm in names(m$value)) {
            hit <- tab34$scenario == s & tab34$metric == nm
            stopifnot(sum(hit) == 1L)
            src <- sub(" \\(p\\. [0-9]+\\)$", "", tab34$source[hit])
            if (grepl("^pts_", nm)) {
              if (b != "contour") next
              rows[[length(rows) + 1L]] <- lu2025_row(
                src, paste0("Scenario ", s, ", ", nm), tab34$value[hit], m,
                nm, n_trials)
            } else {
              rows[[length(rows) + 1L]] <- lu2025_row(
                src, paste0("Scenario ", s, ", ", nm, " (obd_bound = \"", b,
                            "\")"), tab34$value[hit], m, nm, n_trials,
                judge = FALSE,
                reason = paste0("the selection of Tables 3 and 4 differs ",
                                "from that of Figure 2 for the same trials ",
                                "and is not reproduced by the rules tried; "))
            }
          }
        }
        if (b == "contour") {
          for (pn in names(panels)) {
            hit <- fig$panel == pn & fig$scenario == s &
              fig$method == paste0("wT=", w_tox)
            stopifnot(sum(hit) == 1L)
            rows[[length(rows) + 1L]] <- lu2025_row(
              "Lu 2025 SI Figure S2",
              paste0("Scenario ", s, ", ", panels[[pn]], ", wT = ", w_tox),
              round(fig$value_b[hit], 2), m, panels[[pn]], n_trials)
          }
          if (w_tox == 0.11) {
            # The Comb-BOIN12 bars of Figure 2 (both pairs of windows).
            for (wn in names(lu2025_windows())) {
              method <- paste0("Comb-BOIN12 (", wn, ")")
              for (pn in names(panels)) {
                hit <- fig2$panel == pn & fig2$scenario == s &
                  fig2$method == method
                stopifnot(sum(hit) == 1L)
                rows[[length(rows) + 1L]] <- lu2025_row(
                  "Lu 2025 Figure 2",
                  paste0("Scenario ", s, ", ", panels[[pn]], ", ", method),
                  round(fig2$value_b[hit], 2), m, panels[[pn]], n_trials)
              }
            }
          }
        }
      }
    }
  }
  do.call(rbind, rows)
}

# The windows (AT, AE) of Figure 2, in months.
lu2025_windows <- function() {
  list("AT=1, AE=2" = c(1, 2), "AT=2, AE=3" = c(2, 3))
}

# Mean trial duration, with the standard deviation for the tolerance.
lu2025_duration <- function(trials) {
  list(value = c(duration = mean(trials$duration)),
       sd = c(duration = stats::sd(trials$duration)))
}

# TITE-Comb-BOIN12 (panels A-E) and the durations of Comb-BOIN12 (panel E)
# of Figure 2.
reproduce_lu2025_figure2 <- function(n_trials = 5000, seed = 2025) {
  fig2 <- lu2025_read("lu2025-figure2.csv")
  panels <- c(A_sel_obdc = "sel_obdc", B_pts_obdc = "pts_obdc",
              C_sel_overtox = "sel_overtox", D_pts_overtox = "pts_overtox")
  psi <- lu2025_psi(0.11)
  reason_e <- paste0("the accrual and waiting rules are not described in ",
                     "the article; no timing assumption that was tried ",
                     "reproduces the durations; ")
  rows <- list()
  for (sc in lu2025_scenarios()) {
    s <- sc$scenario
    truth <- lu2025_truth(sc, 0.11)
    obdc <- matrix(FALSE, nrow(sc$p_true), ncol(sc$p_true))
    obdc[cbind(truth$dose_a, truth$dose_b)] <- truth$true_obd
    for (wn in names(lu2025_windows())) {
      win <- lu2025_windows()[[wn]]
      for (method in c("TITE-Comb-BOIN12", "Comb-BOIN12")) {
        message(method, " (", wn, "), Scenario ", s, ": ",
                format(Sys.time(), "%H:%M:%S"))
        trials <- tite_comb_boin12_simulate(
          target = sc$target, target_resp = 0.2, p_true = sc$p_true,
          p_resp = sc$p_resp, n_cohort = sc$n_cohort, cohort_size = 3,
          window = win[1L], resp_window = win[2L], accrual_rate = 2,
          accrual = "uniform", event_time = "uniform",
          min_completed_share = if (method == "Comb-BOIN12") 1 else 0.5,
          pending_use = "utility", n_trials = n_trials, extrasafe = TRUE,
          offset = 0.10, psi00 = psi[["psi00"]], psi11 = psi[["psi11"]],
          benchmark_weight = 0.6, seed = seed)
        label <- paste0(method, " (", wn, ")")
        hit <- fig2$panel == "E_duration" & fig2$scenario == s &
          fig2$method == label
        stopifnot(sum(hit) == 1L)
        rows[[length(rows) + 1L]] <- lu2025_row(
          "Lu 2025 Figure 2", paste0("Scenario ", s, ", duration, ", label),
          round(fig2$value_b[hit], 2), lu2025_duration(trials), "duration",
          n_trials, judge = FALSE, reason = reason_e)
        if (method == "Comb-BOIN12") next
        sel <- comb_boin12_select_obd(
          trials$n_pts, trials$n_tox, trials$n_resp, target = sc$target,
          target_resp = 0.2, extrasafe = TRUE, offset = 0.10,
          psi00 = psi[["psi00"]], psi11 = psi[["psi11"]],
          obd_bound = "contour", dose_a = sc$dose_a, dose_b = sc$dose_b,
          early_stop = trials$early_stop)
        m <- lu2025_metrics(trials, sel, sc, obdc)
        for (pn in names(panels)) {
          hit <- fig2$panel == pn & fig2$scenario == s & fig2$method == label
          stopifnot(sum(hit) == 1L)
          rows[[length(rows) + 1L]] <- lu2025_row(
            "Lu 2025 Figure 2",
            paste0("Scenario ", s, ", ", panels[[pn]], ", ", label),
            round(fig2$value_b[hit], 2), m, panels[[pn]], n_trials)
        }
      }
    }
  }
  do.call(rbind, rows)
}

run_reproduce_lu2025 <- function(out_dir = "inst/reproduce/results",
                                 n_trials = 5000) {
  # The functions of the package under development must be loaded.
  if (!exists("tite_comb_boin12_simulate", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- rbind(reproduce_lu2025_utilities(), reproduce_lu2025_figure_s1(),
               reproduce_lu2025_simulation(n_trials = n_trials),
               reproduce_lu2025_figure2(n_trials = n_trials))
  rownames(res) <- NULL
  res <- lu2025_apply_explained(res, "reproduce_Lu2025.R")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Lu2025.csv"),
                   row.names = FALSE)
  res
}
