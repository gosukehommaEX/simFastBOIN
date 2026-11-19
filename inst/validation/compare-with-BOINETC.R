# Comparison of simFastBOIN with the code of Kakizume et al. (2024)
#
# Kakizume T, Takeda K, Taguri M, Morita S. BOIN-ETC: A Bayesian optimal
# interval design considering efficacy and toxicity to identify the optimal
# dose combinations. Statistical Methods in Medical Research 2024; 33:
# 716-727. doi:10.1177/09622802241236936
#
# The code of the authors (https://github.com/kakizume/BOINETC, folder prog)
# is not a package. Its main program, SIM_Main.R, does not run as published
# (sc.list, u.min, OUTFILE and OUT are not defined, and a parenthesis of a
# call of source() is unbalanced). This script sources the files of prog
# with the functions and the scenarios into an environment, defines there the
# objects that SIM_Main.R defines or misses (the settings, OUTFILE and OUT as
# a temporary folder, u.min, which is unused) and runs the loops of
# SIM_Main.R and FUNC_BOINETC.study.R itself, in the order of item 1. It
# replaces readr::read_csv() by utils::read.csv(), which reads the files of
# the code in the same way, and takes biviso() from the package Iso, which
# must be installed. The functions of the code are not changed.
#
# 1. The published values (PASS or FAIL to the printed digits). The trials
#    of the article are those of the code with R's random numbers after
#    set.seed(1234), once for each size of the dose matrix, with the runs in
#    the order BOIN-ETC1, 2, 3 (the functions of get.oc.BOINETC_m1_v02.R,
#    m2 and m3), n_stop 6 then 9, and the scenarios of that size: 1 and 2
#    (2 x 3), 3, 4 and 5 (4 x 4), and 7, 8 and 9 (3 x 5). This order was
#    found with a translation of the code into Python; Scenarios 6 and 10
#    are not in these streams and are not compared here (they are compared
#    in distribution in inst/reproduce/reproduce_Kakizume2024.R). The 48 runs
#    are summarized by the function CALC.SUM() of the code for the three
#    utilities of the article and compared with Table 3 and Supplementary
#    Tables S1 to S9: the mean numbers of patients at every combination
#    (S4, S5), the selection percentages (S2, S3, utility 100-40-60-0) and
#    %ODC, %nODC and #N (Table 3, S1, S6 to S9).
# 2. The selection of simFastBOIN on the same trials (INFO). For every run
#    and utility, the optimal dose combinations (ODCs) of
#    boin_etc_select_obd() with odc_tie = "all" are compared trial by trial
#    with those of the code, which are computed by a copy of the lines of
#    CALC.SUM() (checked against the percentages of CALC.SUM(), PASS or
#    FAIL). The two differ when combinations tie for the MTD combination of
#    a row in the exact isotonic regression of simFastBOIN: the code keeps
#    only those that the convergence errors of Iso::biviso() leave closest
#    to the target. They also differ when the code selects an untreated
#    combination (whose utility it estimates at 0.5) or a combination
#    excluded by its DLT data, which simFastBOIN does not.
# Values of the article that disagree with the other values of the same
# tables (misprints) are listed with the reason in
# inst/validation/explained.csv, and their judgements of FAIL become
# EXPLAINED.
#
# Usage, from the package root, with the code of the authors in a folder
# (about ten minutes):
#
#   devtools::load_all()
#   source("inst/validation/compare-with-BOINETC.R")
#   res <- run_compare_with_boinetc(prog_dir = "path/to/BOINETC-main/prog")
#
# The judgements are written to
# inst/validation/results/compare-with-BOINETC.csv.

boinetc_read <- function(name) {
  utils::read.csv(file.path("inst", "reproduce", name),
                  stringsAsFactors = FALSE)
}

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
boinetc_apply_explained <- function(res, file,
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

boinetc_groups <- function() {
  list("2 x 3" = c(1L, 2L), "4 x 4" = c(3L, 4L, 5L), "3 x 5" = c(7L, 8L, 9L))
}

boinetc_methods <- function() c(m1 = "ETC1", m2 = "ETC2", m3 = "ETC3")

boinetc_utilities <- function() {
  list("100-40-60-0" = c(100, 40, 60, 0), "100-25-75-0" = c(100, 25, 75, 0),
       "100-50-50-0" = c(100, 50, 50, 0))
}

boinetc_n_cohort <- function(sc, n_stop) {
  if (sc <= 2L) c("6" = 8L, "9" = 12L)[[as.character(n_stop)]]
  else c("6" = 16L, "9" = 24L)[[as.character(n_stop)]]
}

# The functions of the code of the authors in an environment, with the
# objects that SIM_Main.R defines or loads.
boinetc_env <- function(prog_dir, work_dir) {
  files <- c("BOINETC.scenario.R", "FUNC_BOINETC.subtrial.R",
             "get.oc.BOINETC_m1_v02.R", "get.oc.BOINETC_m2_v02.R",
             "get.oc.BOINETC_m3_v02.R", "FUNC_CALC.SUM.R")
  miss <- files[!file.exists(file.path(prog_dir, files))]
  if (length(miss) > 0L) {
    stop("not found in 'prog_dir': ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  e <- new.env(parent = globalenv())
  for (f in files) sys.source(file.path(prog_dir, f), envir = e)
  e$biviso <- Iso::biviso
  e$read_csv <- function(file) utils::read.csv(file)
  e$OUTFILE <- e$OUT <- paste0(work_dir, "/")
  e$list.p.truetox <- lapply(1:10, function(s) {
    get(paste0("pt.true.mat", s), envir = e)
  })
  e$list.p.trueeff <- lapply(1:10, function(s) {
    get(paste0("pe.true.mat", s), envir = e)
  })
  e$list.tdose <- list(c(4, 5), c(2, 3), c(8, 11, 14), c(7, 10), c(4, 7),
                       c(7, 10), c(9, 11, 13), c(6, 8, 10), c(9, 11),
                       c(8, 10))
  e$phi <- 0.30
  e$delta <- 0.60
  e$cohortsize <- 3
  e$ntrial <- 1000
  e$lambda1 <- 0.14
  e$lambda2 <- 0.35
  e$eta1 <- 0.48
  e$alpha <- 1
  e$beta <- 1
  e
}

# The scenarios of the code must be those of Table 2.
boinetc_check_scenarios <- function(e) {
  tab <- boinetc_read("kakizume2024-scenarios.csv")
  for (s in 1:10) {
    d <- tab[tab$scenario == s, , drop = FALSE]
    pt <- e$list.p.truetox[[s]]
    pe <- e$list.p.trueeff[[s]]
    idx <- cbind(d$dose_a, d$dose_b)
    odc <- sort(d$dose_a[d$odc == 1] + nrow(pt) * (d$dose_b[d$odc == 1] - 1))
    if (!identical(dim(pt), c(d$n_a[1L], d$n_b[1L])) ||
        max(abs(pt[idx] - d$p_tox)) > 1e-12 ||
        max(abs(pe[idx] - d$p_eff)) > 1e-12 ||
        !identical(as.numeric(odc), as.numeric(e$list.tdose[[s]]))) {
      stop("Scenario ", s, " of the code differs from Table 2", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# One run of the code (1,000 trials), written to the folder OUTFILE.
boinetc_run <- function(e, sc, n_stop, method) {
  fun <- switch(method, m1 = e$get.oc.BOINETC, m2 = e$get.oc.BOINETC_m2,
                m3 = e$get.oc.BOINETC_m3)
  filename <- paste0("BOINETC_", method, ".new.sc", sc, "_cs3_ns", n_stop)
  fun(ncohort = boinetc_n_cohort(sc, n_stop), cohortsize = 3,
      n.earlystop = n_stop, ntrial = 1000, phi = 0.30, delta = 0.60,
      lambda1 = 0.14, lambda2 = 0.35, eta1 = 0.48,
      tdose = e$list.tdose[[sc]], u.min = NA, n.scr = sc,
      pt.true.mat = e$list.p.truetox[[sc]],
      pe.true.mat = e$list.p.trueeff[[sc]], filename = filename)
  filename
}

# CALC.SUM() of the code for one run and one utility.
boinetc_calc <- function(e, filename, sc, n_stop, method, u) {
  e$u <- u
  e$ncohort <- boinetc_n_cohort(sc, n_stop)
  e$n.es <- n_stop
  e$method <- paste0("BOINETC_", method)
  n_sink <- sink.number()
  on.exit(while (sink.number() > n_sink) sink(), add = TRUE)
  e$CALC.SUM(FILENAME = filename, SC = sc)
}

# The trials of one run, as arrays (trial by level of drug A by level of
# drug B).
boinetc_trials <- function(e, filename, n_a, n_b) {
  d <- utils::read.csv(paste0(e$OUT, filename, "_data.csv"))
  d <- as.matrix(d[-(1:2), -1, drop = FALSE])
  nd <- n_a * n_b
  block <- function(k) {
    x <- array(round(d[, (k - 1L) * nd + seq_len(nd)]), c(nrow(d), n_a, n_b))
    storage.mode(x) <- "integer"
    x
  }
  list(n_pts = block(1L), n_tox = block(2L), n_resp = block(3L),
       n_cat = lapply(7:10, block))
}

# The ODCs of one trial as CALC.SUM() of the code selects them (a copy of
# its lines).
boinetc_code_odc <- function(n_pts, n_tox, n_cat, u, phi = 0.30, alpha = 1,
                             beta = 1) {
  nd1 <- nrow(n_pts)
  phat_tox <- Iso::biviso(matrix((as.vector(n_tox) + 0.05) /
                                   (as.vector(n_pts) + 0.1), nrow = nd1))
  dist <- abs(phat_tox - phi)
  dl_mtd <- which(dist == apply(dist, 1, min))
  ind_tol <- matrix(0, nd1, ncol(n_pts))
  for (k in dl_mtd) {
    ind_tol[(k - 1L) %% nd1 + 1L, seq_len((k - 1L) %/% nd1 + 1L)] <- 1
  }
  x <- (u[1] * n_cat[[1]] + u[2] * n_cat[[2]] + u[3] * n_cat[[3]] +
          u[4] * n_cat[[4]]) / 100
  util <- (x + alpha) / (n_pts + alpha + beta) * ind_tol
  which(c(util) == max(util, na.rm = TRUE) & util > 0)
}

boinetc_row <- function(source, item, published, reproduced, digits,
                        judge = TRUE, reason = "") {
  same <- abs(published - reproduced) < 1e-6
  data.frame(source = source, item = item,
             published = formatC(published, format = "f", digits = digits),
             reproduced = formatC(reproduced, format = "f", digits = digits),
             tolerance = 0,
             judgement = if (!judge) "INFO" else if (same) "PASS" else "FAIL",
             reason = reason, stringsAsFactors = FALSE)
}

# Parts 1 and 2 for one run.
boinetc_compare_run <- function(e, filename, sc, n_stop, method, summ,
                                bycomb) {
  design <- boinetc_methods()[[method]]
  pt <- e$list.p.truetox[[sc]]
  n_a <- nrow(pt)
  n_b <- ncol(pt)
  lab <- paste0("A", rep(seq_len(n_a), n_b), "B", rep(seq_len(n_b),
                                                       each = n_a))
  tdose <- e$list.tdose[[sc]]
  prefix <- paste0("Scenario ", sc, ", ", design, ", n_stop = ", n_stop, ", ")
  tr <- boinetc_trials(e, filename, n_a, n_b)
  n_trials <- dim(tr$n_pts)[1L]
  rows <- list()
  add <- function(r) rows[[length(rows) + 1L]] <<- r
  t <- summ[summ$scenario == as.character(sc) & summ$design == design &
              summ$n_stop == n_stop, , drop = FALSE]
  b <- bycomb[bycomb$scenario == sc & bycomb$design == design &
                bycomb$n_stop == n_stop, , drop = FALSE]
  stopifnot(nrow(t) == 3L, nrow(b) == 2L * n_a * n_b)
  for (un in names(boinetc_utilities())) {
    u <- boinetc_utilities()[[un]]
    out <- boinetc_calc(e, filename, sc, n_stop, method, u)
    s1 <- out$out.sum1
    s2 <- out$out.sum2
    # Values by combination, in the order of lab.
    by_comb <- function(stem) {
      m <- as.matrix(s1[, paste0(stem, seq_len(n_b)), drop = FALSE])
      m <- m[match(seq_len(n_a), s1$doseA), , drop = FALSE]
      as.vector(m)
    }
    p_odc <- by_comb("pODC")
    if (un == "100-40-60-0") {
      m_npts <- by_comb("m.NPTS")
      for (i in seq_len(nrow(b))) {
        k <- b$dose_a[i] + n_a * (b$dose_b[i] - 1L)
        src <- paste0("Kakizume 2024 SI ", b$table[i])
        if (b$metric[i] == "mean_n_pts") {
          add(boinetc_row(src, paste0(prefix, lab[k], ", mean patients"),
                          b$value[i], m_npts[k], 2L))
        } else {
          add(boinetc_row(src, paste0(prefix, lab[k], ", selected (%)"),
                          b$value[i], p_odc[k], 1L))
        }
      }
    }
    r <- t[t$utility == un, , drop = FALSE]
    stopifnot(nrow(r) == 1L)
    src <- paste0("Kakizume 2024 ", if (r$file == "SI") "SI " else "",
                  r$table)
    suffix <- paste0(", utility ", un)
    add(boinetc_row(src, paste0(prefix, "%ODC", suffix), r$pct_odc,
                    s2$correct.pODC, 1L))
    add(boinetc_row(src, paste0(prefix, "%nODC", suffix), r$pct_n_odc,
                    s2$target.PPTS, 1L))
    add(boinetc_row(src, paste0(prefix, "#N", suffix), r$mean_n,
                    s2$mean.npat, 1L))

    # Part 2: the copy of the selection of the code, then simFastBOIN.
    code_odc <- lapply(seq_len(n_trials), function(i) {
      boinetc_code_odc(tr$n_pts[i, , ], tr$n_tox[i, , ],
                       lapply(tr$n_cat, function(x) x[i, , ]), u)
    })
    copy_pct <- round(as.vector(table(factor(unlist(code_odc),
                                             levels = seq_len(n_a * n_b)))) /
                        n_trials * 100, 1)
    add(data.frame(
      source = "Copy of the selection of CALC.SUM()",
      item = paste0(prefix, "selected (%) of every combination", suffix),
      published = paste(formatC(p_odc, format = "f", digits = 1),
                        collapse = " "),
      reproduced = paste(formatC(copy_pct, format = "f", digits = 1),
                         collapse = " "),
      tolerance = 0,
      judgement = if (max(abs(copy_pct - p_odc)) < 1e-6) "PASS" else "FAIL",
      reason = "", stringsAsFactors = FALSE))
    sel <- boin_etc_select_obd(
      tr$n_pts, tr$n_tox, tr$n_resp, tr$n_cat[[3L]], target = 0.30,
      psi00 = u[2], psi11 = u[3], odc_tie = "all")
    tied <- matrix(sel$tied, n_trials)
    pkg_odc <- lapply(seq_len(n_trials), function(i) which(tied[i, ]))
    differ <- !mapply(identical, lapply(code_odc, as.integer),
                      lapply(pkg_odc, as.integer))
    tie_mtd <- vapply(seq_len(n_trials), function(i) {
      d <- matrix(abs(sel$tox_est[i, , ] - 0.30), n_a)
      any(rowSums(d - apply(d, 1L, min) <= 1e-8) > 1L)
    }, logical(1))
    # Untreated, or excluded by the DLT data (Pr(p > 0.30) > 0.95).
    other <- mapply(function(i, k) {
      n <- tr$n_pts[i, , ][k]
      y <- tr$n_tox[i, , ][k]
      any(n == 0L | stats::pbeta(0.30, y + 1, n - y + 1,
                                  lower.tail = FALSE) > 0.95)
    }, seq_len(n_trials), code_odc)
    pkg_pct <- 100 * colMeans(tied)
    add(data.frame(
      source = "simFastBOIN on the trials of the code",
      item = paste0(prefix, "%ODC", suffix),
      published = formatC(s2$correct.pODC, format = "f", digits = 1),
      reproduced = formatC(sum(pkg_pct[tdose]), format = "f", digits = 1),
      tolerance = NA_real_, judgement = "INFO",
      reason = sprintf(paste0(
        "ODCs differ in %.1f%% of the trials: %.1f%% with tied MTD ",
        "combinations in the exact isotonic regression, %.1f%% with an ",
        "untreated or excluded ODC of the code, %.1f%% otherwise"),
        100 * mean(differ), 100 * mean(differ & tie_mtd),
        100 * mean(differ & !tie_mtd & other),
        100 * mean(differ & !tie_mtd & !other)),
      stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

compare_with_boinetc <- function(prog_dir) {
  if (!requireNamespace("Iso", quietly = TRUE)) {
    stop("the package Iso is needed", call. = FALSE)
  }
  work_dir <- tempfile("boinetc")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  e <- boinetc_env(normalizePath(prog_dir, winslash = "/", mustWork = TRUE),
                   normalizePath(work_dir, winslash = "/"))
  boinetc_check_scenarios(e)
  summ <- boinetc_read("kakizume2024-summary.csv")
  bycomb <- boinetc_read("kakizume2024-by-combination.csv")
  old_seed <- get_random_seed()
  on.exit(restore_random_seed(old_seed), add = TRUE)
  rows <- list()
  for (g in names(boinetc_groups())) {
    set.seed(1234)
    for (method in names(boinetc_methods())) {
      for (n_stop in c(6L, 9L)) {
        for (sc in boinetc_groups()[[g]]) {
          message("Code of the authors, ", boinetc_methods()[[method]],
                  ", Scenario ", sc, ", n_stop = ", n_stop, ": ",
                  format(Sys.time(), "%H:%M:%S"))
          filename <- boinetc_run(e, sc, n_stop, method)
          rows[[length(rows) + 1L]] <- boinetc_compare_run(
            e, filename, sc, n_stop, method, summ, bycomb)
        }
      }
    }
  }
  do.call(rbind, rows)
}

run_compare_with_boinetc <- function(prog_dir,
                                     out_dir = "inst/validation/results") {
  # The functions of the package under development must be loaded.
  if (!exists("boin_etc_select_obd", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- compare_with_boinetc(prog_dir)
  rownames(res) <- NULL
  res <- boinetc_apply_explained(res, "compare-with-BOINETC.R")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "compare-with-BOINETC.csv"),
                   row.names = FALSE)
  res
}
