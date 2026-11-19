# Reproduction of Kakizume et al. (2024) with simFastBOIN
#
# Kakizume T, Takeda K, Taguri M, Morita S. BOIN-ETC: A Bayesian optimal
# interval design considering efficacy and toxicity to identify the optimal
# dose combinations. Statistical Methods in Medical Research 2024; 33:
# 716-727. doi:10.1177/09622802241236936
#
# The published values are in inst/reproduce: kakizume2024-scenarios.csv
# (Table 2: the DLT and response probabilities of the 10 scenarios and the
# true optimal dose combinations, ODCs, in boldface), kakizume2024-summary.csv
# (Table 3 and Supplementary Tables S1 and S6 to S9: %ODC, %nODC and #N of
# every scenario and design) and kakizume2024-by-combination.csv
# (Supplementary Tables S2 to S5: the percentages of trials selecting every
# combination as an ODC and the mean numbers of patients there). They were
# read from the text of the PDFs in two independent ways, which agree on all
# 3,240 values, and the boldface from the page images.
#
# The BOIN-ETC1, 2 and 3 columns are reproduced; Cai's design is not part of
# the package. The Mean rows, averages of the rows above them, are not
# compared. The trials of the article (1,000 for each setting) are compared
# with n_trials trials of the package, with a tolerance of four standard
# errors of the difference between two independent simulations plus the
# rounding of the published value (half of its last digit, times the number
# of combinations it adds up for %ODC and %nODC, which the code of the
# authors rounds combination by combination):
#   - the allocation, PASS or FAIL: the mean numbers of patients at every
#     combination (Tables S4 and S5), %nODC (the average over trials of the
#     percentage of the patients of a trial treated at an ODC) and #N (the
#     mean number of patients) of Tables 3, S1 and S6 to S9;
#   - the selection, INFO: %ODC of Tables 3, S1 and S6 to S9 and the
#     selection percentages of Tables S2 and S3. The code of the authors
#     estimates the DLT probabilities with Iso::biviso(), whose convergence
#     errors decide which of the combinations tied in the exact isotonic
#     regression is the MTD combination of a row; the package computes the
#     exact regression, so that its selection differs from the published one
#     by more than the Monte Carlo error in some settings.
# Values of the article that disagree with the other values of the same
# tables (misprints) are listed with the reason in
# inst/validation/explained.csv, and their judgements of FAIL become
# EXPLAINED.
#
# Settings of Section 3.1 of the article: a target DLT probability of 0.30
# and a target response probability of 0.60, the boundaries
# (lambda1, lambda2, eta1) = (0.14, 0.35, 0.48) of boin_et_boundary(), cohorts
# of 3 from A1B1, n_stop (n_earlystop) of 6 and 9 with at most 24 and 36
# patients in Scenarios 1 and 2 and 48 and 72 in Scenarios 3 to 10, the
# exclusion of a combination when Pr(p > 0.30) > 0.95 with a beta(1, 1)
# prior, and the utility scores (u1, u2, u3, u4) = (100, 40, 60, 0) in Table
# 3 and Tables S1 to S5, (100, 25, 75, 0) in Tables S6 and S8 and
# (100, 50, 50, 0) in Tables S7 and S9 (psi00 = u2, psi11 = u3). The DLT and
# the response are independent (resp_cor = 0), as in the code of the
# authors. BOIN-ETC1 is simulated with etc1_start = "below" and rule4 =
# "authors_code", the rules of the code of the authors that produced the
# published values (with the defaults "diagonal" and "paper", the rules of
# the text of the article, the allocation differs; see ?boin_etc_simulate).
# Tied ODCs are all counted (odc_tie = "all"), as in the article, so that
# the selection percentages of a row of Tables S2 and S3 add up to more than
# 100.
#
# Usage, from the package root (a few minutes):
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Kakizume2024.R")
#   res <- run_reproduce_kakizume2024()
#
# The judgements are written to
# inst/reproduce/results/reproduce_Kakizume2024.csv.

kakizume2024_read <- function(name) {
  utils::read.csv(file.path("inst", "reproduce", name),
                  stringsAsFactors = FALSE)
}

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
kakizume2024_apply_explained <- function(res, file,
                                         path = file.path("inst",
                                                          "validation",
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

# Scenarios 1-10 of Table 2, with the true ODCs (boldface).
kakizume2024_scenarios <- function(tab = kakizume2024_read(
                                     "kakizume2024-scenarios.csv")) {
  lapply(sort(unique(tab$scenario)), function(s) {
    d <- tab[tab$scenario == s, , drop = FALSE]
    n_a <- d$n_a[1L]
    n_b <- d$n_b[1L]
    cell <- function(v) {
      m <- matrix(NA_real_, n_a, n_b)
      m[cbind(d$dose_a, d$dose_b)] <- v
      m
    }
    list(scenario = s, p_true = cell(d$p_tox), p_resp = cell(d$p_eff),
         odc = cell(d$odc) == 1,
         n_cohort = if (n_a == 2L && n_b == 3L) c("6" = 8L, "9" = 12L)
                    else c("6" = 16L, "9" = 24L))
  })
}

kakizume2024_designs <- function() {
  list(
    ETC1 = list(subtrials = "etc1", etc1_start = "below",
                rule4 = "authors_code"),
    ETC2 = list(subtrials = "etc2", etc1_start = "diagonal", rule4 = "paper"),
    ETC3 = list(subtrials = "etc3", etc1_start = "diagonal", rule4 = "paper")
  )
}

kakizume2024_utilities <- function() {
  list("100-40-60-0" = c(psi00 = 40, psi11 = 60),
       "100-25-75-0" = c(psi00 = 25, psi11 = 75),
       "100-50-50-0" = c(psi00 = 50, psi11 = 50))
}

kakizume2024_label <- function(n_a, n_b) {
  paste0("A", rep(seq_len(n_a), n_b), "B", rep(seq_len(n_b), each = n_a))
}

# Allocation and selection of one setting (scenario, n_stop, design), with
# the standard deviations of the per-trial values for the tolerances.
kakizume2024_simulate <- function(sc, n_stop, design, n_trials, seed) {
  des <- kakizume2024_designs()[[design]]
  trials <- boin_etc_simulate(
    target = 0.30, target_resp = 0.60, p_true = sc$p_true,
    p_resp = sc$p_resp, n_cohort = sc$n_cohort[[as.character(n_stop)]],
    cohort_size = 3, n_trials = n_trials, subtrials = des$subtrials,
    n_earlystop = n_stop, etc1_start = des$etc1_start, rule4 = des$rule4,
    seed = seed)
  n_cells <- matrix(trials$n_pts, n_trials)
  n_total <- rowSums(n_cells)
  odc <- which(as.vector(sc$odc))
  pct_odc_pts <- 100 * rowSums(n_cells[, odc, drop = FALSE]) / n_total
  out <- list(
    mean_n_pts = colMeans(n_cells),
    sd_n_pts = apply(n_cells, 2L, stats::sd),
    pct_n_odc = mean(pct_odc_pts), sd_pct_n_odc = stats::sd(pct_odc_pts),
    mean_n = mean(n_total), sd_mean_n = stats::sd(n_total),
    selection = list()
  )
  for (un in names(kakizume2024_utilities())) {
    psi <- kakizume2024_utilities()[[un]]
    sel <- boin_etc_select_obd(
      trials$n_pts, trials$n_tox, trials$n_resp, trials$n_tox_resp,
      target = 0.30, psi00 = psi[["psi00"]], psi11 = psi[["psi11"]],
      odc_tie = "all", excluded = trials$excluded,
      early_stop = trials$early_stop)
    tied <- matrix(sel$tied, n_trials)
    n_odc <- 100 * rowSums(tied[, odc, drop = FALSE])
    out$selection[[un]] <- list(pct_selected = 100 * colMeans(tied),
                                pct_odc = mean(n_odc),
                                sd_pct_odc = stats::sd(n_odc))
  }
  out
}

kakizume2024_row <- function(source, item, published, reproduced, se,
                             rounding, judge, n_digits) {
  tol <- 4 * se + rounding
  within <- abs(reproduced - published) <= tol + 1e-9
  if (judge) {
    judgement <- if (within) "PASS" else "FAIL"
    reason <- ""
  } else {
    judgement <- "INFO"
    reason <- paste0("the selection of the code of the authors depends on ",
                     "the convergence errors of Iso::biviso(); ",
                     if (within) "within" else "outside", " the tolerance")
  }
  data.frame(source = source, item = item,
             published = formatC(published, format = "f", digits = n_digits),
             reproduced = formatC(reproduced, format = "f",
                                  digits = n_digits + 1L),
             tolerance = round(tol, 3), judgement = judgement,
             reason = reason, stringsAsFactors = FALSE)
}

reproduce_kakizume2024_simulation <- function(n_trials = 10000,
                                              seed = 2024) {
  summ <- kakizume2024_read("kakizume2024-summary.csv")
  bycomb <- kakizume2024_read("kakizume2024-by-combination.csv")
  factor_se <- sqrt(1 / 1000 + 1 / n_trials)
  rows <- list()
  add <- function(r) rows[[length(rows) + 1L]] <<- r
  for (sc in kakizume2024_scenarios()) {
    s <- sc$scenario
    lab <- kakizume2024_label(nrow(sc$p_true), ncol(sc$p_true))
    n_odc <- sum(sc$odc)
    for (n_stop in c(6L, 9L)) {
      for (design in names(kakizume2024_designs())) {
        message("BOIN-", design, ", Scenario ", s, ", n_stop = ", n_stop,
                ": ", format(Sys.time(), "%H:%M:%S"))
        m <- kakizume2024_simulate(sc, n_stop, design, n_trials, seed)
        prefix <- paste0("Scenario ", s, ", ", design, ", n_stop = ",
                         n_stop, ", ")
        # Tables S2 to S5.
        b <- bycomb[bycomb$scenario == s & bycomb$design == design &
                      bycomb$n_stop == n_stop, , drop = FALSE]
        stopifnot(nrow(b) == 2L * length(lab))
        for (i in seq_len(nrow(b))) {
          k <- b$dose_a[i] + nrow(sc$p_true) * (b$dose_b[i] - 1L)
          src <- paste0("Kakizume 2024 SI ", b$table[i])
          if (b$metric[i] == "mean_n_pts") {
            add(kakizume2024_row(
              src, paste0(prefix, lab[k], ", mean patients"), b$value[i],
              m$mean_n_pts[k], m$sd_n_pts[k] * factor_se, 0.005, TRUE, 2L))
          } else {
            rep <- m$selection[["100-40-60-0"]]$pct_selected[k]
            q <- max(b$value[i], rep, 0.5) / 100
            add(kakizume2024_row(
              src, paste0(prefix, lab[k], ", selected (%)"), b$value[i],
              rep, 100 * sqrt(q * (1 - q)) * factor_se, 0.05, FALSE, 1L))
          }
        }
        # Tables 3, S1 and S6 to S9.
        t <- summ[summ$scenario == as.character(s) &
                    summ$design == design & summ$n_stop == n_stop, ,
                  drop = FALSE]
        stopifnot(nrow(t) == 3L)
        for (i in seq_len(nrow(t))) {
          src <- paste0("Kakizume 2024 ",
                        if (t$file[i] == "SI") "SI " else "", t$table[i])
          util <- m$selection[[t$utility[i]]]
          suffix <- paste0(", utility ", t$utility[i])
          add(kakizume2024_row(
            src, paste0(prefix, "%ODC", suffix), t$pct_odc[i], util$pct_odc,
            util$sd_pct_odc * factor_se, 0.05 * n_odc, FALSE, 1L))
          add(kakizume2024_row(
            src, paste0(prefix, "%nODC", suffix), t$pct_n_odc[i],
            m$pct_n_odc, m$sd_pct_n_odc * factor_se, 0.05 * n_odc, TRUE,
            1L))
          add(kakizume2024_row(
            src, paste0(prefix, "#N", suffix), t$mean_n[i], m$mean_n,
            m$sd_mean_n * factor_se, 0.05, TRUE, 1L))
        }
      }
    }
  }
  do.call(rbind, rows)
}

run_reproduce_kakizume2024 <- function(out_dir = "inst/reproduce/results",
                                       n_trials = 10000, seed = 2024) {
  # The functions of the package under development must be loaded.
  if (!exists("boin_etc_simulate", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- reproduce_kakizume2024_simulation(n_trials = n_trials, seed = seed)
  rownames(res) <- NULL
  res <- kakizume2024_apply_explained(res, "reproduce_Kakizume2024.R")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Kakizume2024.csv"),
                   row.names = FALSE)
  res
}
