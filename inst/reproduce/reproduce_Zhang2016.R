# Reproduction of Zhang and Yuan (2016) with simFastBOIN
#
# Zhang L, Yuan Y. A practical Bayesian design to identify the maximum
# tolerated dose contour for drug combination trials. Statistics in Medicine
# 2016; 35(27): 4924-4936. doi:10.1002/sim.7095
#
# Section 4.1 of the article simulates the waterfall design in the 14
# scenarios of Table III (dose matrices of 2 x 3, 4 x 4 and 3 x 5
# combinations) with a target DLT probability of 0.3, cohorts of 3, the end
# of a subtrial at 12 patients at a dose and at most four patients per dose
# of a subtrial, and 1000 trials. Table IV gives the percentage of trials
# selecting each combination as an MTD, Table V the percentage of trials
# selecting exactly the MTD contour and the percentages of patients treated
# above and on the contour, and Table S1 of the supporting information the
# average number of patients at each combination. These are 398 cells for
# the waterfall design; the cells of the PIPE design are not reproduced.
#
# BOIN::get.oc.comb(mtd.contour = TRUE) of 'BOIN' 2.7.2 was read and written
# again in Python with R's random numbers
# (02_Rpackage/Validation/waterfall of the development files, not part of
# the package). With the default seed of that function, 6, its default
# numbers of cohorts and 1000 trials, the transcription gives all 398 cells
# as printed: the article was produced with these rules and this seed. The
# cells are therefore judged PASS when they are equal to the printed digits,
# and FAIL otherwise; there is no allowance for simulation error.
#
#   - The percentages are rounded to one decimal.
#   - The average numbers of patients are rounded as the reference
#     implementation prints them, to two decimals, and then half up to the
#     one decimal of Table S1 (2.649 becomes 2.65 and then 2.7).
#
# The judged rows use boin_compatible = TRUE: when a subtrial other than the
# one added after the first ends without a candidate MTD, the reference
# implementation drops its data, and so do the published numbers of
# patients. The rows with the default rules of the package, which keep
# these data, are INFO.
#
# The sample size of Table V, which the article matched between the two
# designs, is compared with the multiple of 3 closest to the average number
# of patients treated, the patients of the dropped subtrials included; these
# rows are INFO, because the article does not say how the sizes were
# matched.
#
# Scenario 7 of Table III is not monotone (0.66 at A2B4 and 0.64 at A3B4);
# it is simulated as printed.
#
# Usage, from the package root:
#
#   devtools::load_all()
#   source("inst/reproduce/reproduce_Zhang2016.R")
#   res <- run_reproduce_zhang2016()
#   table(res$design, res$judgement)
#
# The judgements are written to
# inst/reproduce/results/reproduce_Zhang2016.csv.

# The cells of the article: inst/reproduce/zhang2016-tables.csv, written
# from the text of the article and of its supporting information.
zhang2016_tables <- function(path = file.path("inst", "reproduce",
                                              "zhang2016-tables.csv")) {
  utils::read.csv(path, stringsAsFactors = FALSE,
                  colClasses = c("character", "integer", "character",
                                 "integer", "integer", "character"))
}

# The matrix of true DLT probabilities of a scenario of Table III.
zhang2016_scenario <- function(pub, scenario) {
  cells <- pub[pub$table == "Table III" & pub$scenario == scenario, ,
               drop = FALSE]
  p_true <- matrix(NA_real_, max(cells$dose_a), max(cells$dose_b))
  p_true[cbind(cells$dose_a, cells$dose_b)] <- as.numeric(cells$published)
  p_true
}

# An average number of patients as Table S1 prints it: two decimals as the
# reference implementation formats it, then half up to one decimal.
zhang2016_patients <- function(x) {
  two <- round(as.numeric(formatC(x, digits = 2, format = "f")) * 100)
  floor((two + 5) / 10) / 10
}

# The quantities of the tables from one simulation, as printed.
zhang2016_quantities <- function(oc, cells) {
  value <- rep(NA_real_, nrow(cells))
  at <- cbind(cells$dose_a, cells$dose_b)
  sel <- cells$quantity == "Selection (%)"
  pts <- cells$quantity == "Patients"
  value[sel] <- round(oc$sel_percent[at[sel, , drop = FALSE]], 1)
  value[pts] <- zhang2016_patients(oc$n_pts_dose[at[pts, , drop = FALSE]])
  value[cells$quantity == "Contour selected (%)"] <-
    round(oc$pct_contour_selected, 1)
  value[cells$quantity == "Patients above contour (%)"] <-
    round(oc$pct_patients_above, 1)
  value[cells$quantity == "Patients on contour (%)"] <-
    round(oc$pct_patients_contour, 1)
  value[cells$quantity == "Sample size"] <-
    3 * round(oc$total_n_enrolled / 3)
  value
}

# sim_waterfall() warns when a matrix is not monotone, which Scenario 7 is
# as printed; only that warning is silenced.
zhang2016_simulate <- function(p_true, n_trials, seed, compatible) {
  withCallingHandlers(
    sim_waterfall(target = 0.3, p_true = p_true, cohort_size = 3,
                  n_trials = n_trials, n_earlystop = 12,
                  boin_compatible = compatible, seed = seed),
    warning = function(w) {
      if (grepl("not non-decreasing", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

reproduce_zhang2016_tables <- function(n_trials = 1000, seed = 6,
                                       tables = file.path(
                                         "inst", "reproduce",
                                         "zhang2016-tables.csv"
                                       )) {
  pub <- zhang2016_tables(tables)
  runs <- list(
    list(label = "Waterfall", compat = TRUE, judged = TRUE),
    list(label = "Waterfall (rules of the package)", compat = FALSE,
         judged = FALSE)
  )
  rows <- list()
  for (run in runs) {
    for (s in sort(unique(pub$scenario))) {
      p_true <- zhang2016_scenario(pub, s)
      oc <- zhang2016_simulate(p_true, n_trials, seed, run$compat)
      cells <- pub[pub$table != "Table III" & pub$scenario == s, ,
                   drop = FALSE]
      where <- ifelse(is.na(cells$dose_a), "",
                      paste0(" A", cells$dose_a, "B", cells$dose_b))
      rows[[length(rows) + 1L]] <- data.frame(
        item = paste0(cells$table, ", Scenario ", s), design = run$label,
        quantity = paste0(cells$quantity, where),
        published = as.numeric(cells$published),
        reproduced = zhang2016_quantities(oc, cells), tolerance = 0,
        judged = run$judged & cells$quantity != "Sample size",
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# Judgements of FAIL listed in inst/validation/explained.csv, with the
# reason, become EXPLAINED.
zhang2016_apply_explained <- function(res, file,
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

run_reproduce_zhang2016 <- function(n_trials = 1000, seed = 6,
                                    out_dir = "inst/reproduce/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_waterfall", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  res <- reproduce_zhang2016_tables(n_trials = n_trials, seed = seed)
  rownames(res) <- NULL
  equal <- abs(res$published - res$reproduced) <= res$tolerance + 1e-9
  res$judgement <- ifelse(!res$judged, "INFO",
                          ifelse(equal, "PASS", "FAIL"))
  res$judged <- NULL
  res <- zhang2016_apply_explained(res, "reproduce_Zhang2016.R")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "reproduce_Zhang2016.csv"),
                   row.names = FALSE)
  res
}
