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
# Usage, from the package root:
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

run_reproduce_lin2020 <- function(out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("boin12_decision_table", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- rbind(reproduce_lin2020_table3(), reproduce_lin2020_numbers())
  rownames(res) <- NULL
  res <- apply_explained(res, "reproduce_Lin2020.R", res$item, res$source)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Lin2020.csv"),
                   row.names = FALSE)
  res
}
