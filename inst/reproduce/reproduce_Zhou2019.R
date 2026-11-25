# Reproduction of Zhou, Lee and Yuan (2019) with simFastBOIN
#
# Zhou Y, Lee JJ, Yuan Y. A utility-based Bayesian optimal interval (U-BOIN)
# phase I/II design to identify the optimal biological dose for targeted and
# immune therapies. Statistics in Medicine 2019; 38(28): S5299-S5316.
# doi:10.1002/sim.8361
#
# There is no reference implementation (the article offers a web
# application). The script compares the package with the article in two
# ways.
#
# Values that do not depend on random numbers:
#   - Table 2: the escalation and de-escalation boundaries lambda_e and
#     lambda_d of BOIN for six targets, to the three printed decimals. The
#     article prints 0.358 for lambda_d at the target 0.30, which is 0.35852;
#     this cell is registered in inst/validation/explained.csv.
#   - Table S2 of the supporting information: the boundaries of stage I for
#     1 to 24 patients (target 0.25, elimination at the DLT limit 0.30).
#   - Tables S4 to S6: the utility of a dose with 3, 6 and 9 patients and
#     the utilities 0, 15, 25, 100, to the printed decimal, and whether the
#     dose is admissible. The tables mark a dose without any response
#     inadmissible, which the efficacy prior Beta(0.05, 0.05) reproduces
#     (rows "U-BOIN, efficacy prior Beta(0.05, 0.05)"). With the default
#     Beta(1, 1) prior, which the simulations of the article support, 12
#     counts differ; they are registered in inst/validation/explained.csv.
#     Table S6 does not print the count of 6 responses, 2 DLTs and 6
#     responses without a DLT, which is not compared.
#   - The true utilities of Table 3, Table S9 and Table S8, from the Gumbel
#     model with the association parameter 0.2, rounded to an integer.
#
# Simulations: Table 3 (Simulation A, eight scenarios), Table S9 (scenarios
# A1 to A8) and the U-BOIN-CD rows of Table S8 (N = 39): the selection
# percentages, the percentage of early stopping (no OBD) and the average
# numbers of patients, 264 cells, with 54 patients in cohorts of 3 (39 for
# Table S8), the DLT limit 0.30, the response limit 0.20, the utilities of
# example 1, s1 = 12, s2 = N and the Gumbel model with c = 0.2. A cell is
# judged PASS when it lies within four standard errors of the difference
# between the 2000 trials of the article and the n_trials simulated here,
# plus half a unit of the printed digit, and FAIL otherwise.
#
# Two points of the article were found before the package code was written,
# with an implementation in Python (02_Rpackage/Validation/uboin of the
# development files, not part of the package):
#   - The results of Scenario 7 of Table 3 agree with the rates of Scenario
#     7 of Table S8 (DLT 0.03, 0.08, 0.25, 0.40, 0.55; efficacy 0.05, 0.25,
#     0.35, 0.40, 0.45), not with the rates printed in Table 3. The judged
#     rows use the rates of Table S8; the rows "U-BOIN, rates of Table 3 as
#     printed" use the printed rates, and their cells outside the allowance
#     are registered in inst/validation/explained.csv.
#   - The text of the supporting information gives s1 = 9 for Table S8, but
#     its U-BOIN-CD rows agree with s1 = 12. The judged rows use s1 = 12; the
#     rows "U-BOIN, s1 = 9 as in the text" use s1 = 9, and their cells outside
#     the allowance are registered in inst/validation/explained.csv.
#
# The published values are in inst/reproduce/zhou2019-tables.csv and
# inst/reproduce/zhou2019-utility-tables.csv, written from the text of the
# article and of its supporting information.
#
# Usage, from the package root:
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Zhou2019.R")
#   res <- run_reproduce_zhou2019()
#   table(res$design, res$judgement)
#
# The judgements are written to
# inst/reproduce/results/reproduce_Zhou2019.csv.

zhou2019_tables <- function(path = file.path("inst", "reproduce",
                                             "zhou2019-tables.csv")) {
  # Every column is read as text; "NA" (no elimination boundary) and empty
  # doses stay text.
  utils::read.csv(path, stringsAsFactors = FALSE, na.strings = character(0),
                  colClasses = c("character", "character", "character",
                                 "character", "character"))
}

zhou2019_utility_tables <- function(path = file.path(
                                      "inst", "reproduce",
                                      "zhou2019-utility-tables.csv"
                                    )) {
  utils::read.csv(path, stringsAsFactors = FALSE)
}

# The DLT and efficacy rates of a scenario.
zhou2019_rates <- function(pub, table, scenario) {
  cells <- pub[pub$table == table & pub$scenario == scenario, , drop = FALSE]
  get <- function(q) {
    x <- cells[cells$quantity == q, , drop = FALSE]
    as.numeric(x$published[order(as.integer(x$dose))])
  }
  list(p_true = get("DLT rate"), p_resp = get("Efficacy rate"))
}

# Half a unit of the last digit of a number as written in the file.
zhou2019_half_unit <- function(text) {
  digits <- ifelse(grepl(".", text, fixed = TRUE),
                   nchar(sub("^[^.]*\\.", "", text)), 0L)
  0.5 * 10^(-digits)
}

zhou2019_row <- function(item, design, quantity, published, reproduced,
                         tolerance) {
  data.frame(item = item, design = design, quantity = quantity,
             published = published, reproduced = reproduced,
             tolerance = tolerance, judged = TRUE, stringsAsFactors = FALSE)
}

reproduce_zhou2019_boundaries <- function(pub) {
  rows <- list()
  t2 <- pub[pub$table == "Table 2", , drop = FALSE]
  for (i in seq_len(nrow(t2))) {
    target <- as.numeric(t2$scenario[i])
    lam <- boin_lambda(target)
    rows[[length(rows) + 1L]] <- zhou2019_row(
      "Table 2", "BOIN", paste0(t2$quantity[i], ", target ", t2$scenario[i]),
      as.numeric(t2$published[i]), lam[[t2$quantity[i]]],
      zhou2019_half_unit(t2$published[i])
    )
  }
  # The boundaries that uboin_simulate() uses, for 24 patients at a dose.
  trials <- uboin_simulate(tox_upper = 0.30, resp_lower = 0.20,
                           p_true = rep(0.1, 5), p_resp = rep(0.5, 5),
                           n_cohort = 8, cohort_size = 3, n_trials = 1,
                           seed = 1)
  b <- trials$boundary
  s2 <- pub[pub$table == "Table S2", , drop = FALSE]
  ours <- list(escalate = b$b_esc, deescalate = b$b_deesc,
               eliminate = b$b_elim)
  n <- as.integer(s2$dose)
  value <- mapply(function(q, k) ours[[q]][k], s2$quantity, n)
  # A missing elimination boundary is NA in both; it is compared as -1.
  value[is.na(value)] <- -1
  published <- rep(-1, nrow(s2))
  shown <- s2$published != "NA"
  published[shown] <- as.numeric(s2$published[shown])
  rows[[length(rows) + 1L]] <- zhou2019_row(
    "Table S2", "U-BOIN", paste0(s2$quantity, ", ", n, " patients"),
    published, as.numeric(value), 0
  )
  do.call(rbind, rows)
}

reproduce_zhou2019_utility_tables <- function(upub) {
  key <- function(d) paste(d$n_resp, d$n_tox, d$n_resp_no_tox)
  rows <- list()
  for (n in sort(unique(upub$n))) {
    cells <- upub[upub$n == n, , drop = FALSE]
    item <- paste0(cells$table[1L], ", ", n, " patients")
    where <- paste0("responses ", cells$n_resp, ", DLTs ", cells$n_tox,
                    ", responses without DLT ", cells$n_resp_no_tox)
    for (prior in list(c(0.05, 0.05), c(1, 1))) {
      tab <- uboin_utility_table(n, tox_upper = 0.30, resp_lower = 0.20,
                                 utility = c(0, 15, 25, 100),
                                 prior_resp = prior)
      at <- match(key(cells), key(tab))
      design <- if (prior[1L] == 1) {
        "U-BOIN"
      } else {
        "U-BOIN, efficacy prior Beta(0.05, 0.05)"
      }
      rows[[length(rows) + 1L]] <- zhou2019_row(
        item, design, paste0("Admissible, ", where),
        as.numeric(cells$admissible), as.numeric(tab$admissible[at]), 0
      )
      if (prior[1L] == 1) {
        shown <- !is.na(cells$utility)
        rows[[length(rows) + 1L]] <- zhou2019_row(
          item, "U-BOIN", paste0("Utility, ", where[shown]),
          cells$utility[shown], tab$utility[at][shown], 0.05
        )
      }
    }
  }
  do.call(rbind, rows)
}

# The quantities of the tables from one simulation, with the standard
# deviations over the trials that the allowances of the averages need.
zhou2019_quantities <- function(oc) {
  n <- oc$trials$n_pts
  dose_names <- paste0("DL", seq_len(ncol(n)))
  value <- c(unname(oc$sel_percent), oc$percent_no_obd,
             unname(oc$n_pts_dose))
  # NA marks a percentage of trials, whose allowance comes from the
  # binomial variance.
  spread <- c(rep(NA_real_, ncol(n) + 1L), apply(n, 2L, stats::sd))
  names(value) <- c(paste("Selection (%)", dose_names), "Early stopping (%)",
                    paste("Patients", dose_names))
  list(value = value, spread = spread)
}

reproduce_zhou2019_simulations <- function(pub, n_trials = 10000,
                                           seed = 2019, n_pub = 2000) {
  factor_se <- sqrt(1 / n_pub + 1 / n_trials)
  rows <- list()
  for (table in c("Table 3", "Table S9", "Table S8")) {
    scen <- unique(pub$scenario[pub$table == table])
    n_cohort <- if (table == "Table S8") 13 else 18
    for (s in scen) {
      runs <- list(list(design = "U-BOIN", source = c(table, s),
                        n_stage1 = 12))
      if (table == "Table 3" && s == "7") {
        runs[[1L]]$source <- c("Table S8", "7")
        runs[[2L]] <- list(design = "U-BOIN, rates of Table 3 as printed",
                           source = c(table, s), n_stage1 = 12)
      }
      if (table == "Table S8") {
        runs[[2L]] <- list(design = "U-BOIN, s1 = 9 as in the text",
                           source = c(table, s), n_stage1 = 9)
      }
      cells <- pub[pub$table == table & pub$scenario == s, , drop = FALSE]
      item <- paste0(table, ", Scenario ", s)

      # The true utilities, from the rates as printed.
      printed <- zhou2019_rates(pub, table, s)
      utility <- cells[cells$quantity == "Utility", , drop = FALSE]
      utility <- utility[order(as.integer(utility$dose)), , drop = FALSE]
      oc <- sim_uboin(tox_upper = 0.30, resp_lower = 0.20,
                      p_true = printed$p_true, p_resp = printed$p_resp,
                      n_cohort = n_cohort, cohort_size = 3,
                      resp_assoc = 0.2, n_trials = 1, seed = seed)
      rows[[length(rows) + 1L]] <- zhou2019_row(
        item, "U-BOIN", paste0("Utility DL", utility$dose),
        as.numeric(utility$published), round(unname(oc$true_utility)), 0
      )

      sim <- cells[cells$quantity %in% c("Selection (%)",
                                         "Early stopping (%)", "Patients"), ,
                   drop = FALSE]
      label <- ifelse(sim$dose == "", sim$quantity,
                      paste0(sim$quantity, " DL", sim$dose))
      published <- as.numeric(sim$published)
      for (run in runs) {
        rates <- zhou2019_rates(pub, run$source[1L], run$source[2L])
        oc <- sim_uboin(tox_upper = 0.30, resp_lower = 0.20,
                        p_true = rates$p_true, p_resp = rates$p_resp,
                        n_cohort = n_cohort, cohort_size = 3,
                        resp_assoc = 0.2, n_stage1 = run$n_stage1,
                        n_trials = n_trials, keep_trials = TRUE, seed = seed)
        q <- zhou2019_quantities(oc)
        at <- match(label, names(q$value))
        reproduced <- unname(q$value[at])
        spread <- q$spread[at]
        p <- (published + reproduced) / 200
        tolerance <- ifelse(is.na(spread),
                            400 * sqrt(pmax(p * (1 - p), 0)) * factor_se,
                            4 * spread * factor_se) +
          zhou2019_half_unit(sim$published)
        rows[[length(rows) + 1L]] <- zhou2019_row(
          item, run$design, label, published, reproduced, tolerance
        )
      }
    }
  }
  do.call(rbind, rows)
}

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
zhou2019_apply_explained <- function(res, file,
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

run_reproduce_zhou2019 <- function(n_trials = 10000, seed = 2019,
                                   out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_uboin", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  pub <- zhou2019_tables()
  res <- rbind(reproduce_zhou2019_boundaries(pub),
               reproduce_zhou2019_utility_tables(zhou2019_utility_tables()),
               reproduce_zhou2019_simulations(pub, n_trials = n_trials,
                                              seed = seed))
  rownames(res) <- NULL
  within <- abs(res$published - res$reproduced) <= res$tolerance + 1e-9
  res$judgement <- ifelse(!res$judged, "INFO",
                          ifelse(within, "PASS", "FAIL"))
  res$judged <- NULL
  res <- zhou2019_apply_explained(res, "reproduce_Zhou2019.R")
  res$reproduced <- round(res$reproduced, 3)
  res$tolerance <- round(res$tolerance, 3)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Zhou2019.csv"),
                   row.names = FALSE)
  res
}
