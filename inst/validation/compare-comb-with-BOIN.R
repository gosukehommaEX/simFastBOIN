# Comparison of the drug-combination design of simFastBOIN with the BOIN package
# -----------------------------------------------------------------------------
#
# The two implementations are expected to agree EXACTLY, not within a Monte
# Carlo tolerance. sim_comb_boin() draws its random variates in the order of
# BOIN::get.oc.comb() and calls the same bivariate isotonic regression,
# Iso::biviso(), so for a given seed the two produce the same trials and the
# same selections. BOIN::get.oc.comb() rounds its matrices to two decimals and
# its totals to one; the comparison rounds the results of simFastBOIN in the
# same way, so any non-zero difference in this report is a defect to
# investigate, not sampling noise.
#
# Usage
#
#   library(simFastBOIN)
#   source(system.file("validation", "compare-comb-with-BOIN.R",
#                      package = "simFastBOIN"))
#
#   results <- compare_comb_with_boin(quick = TRUE)   # 16 configurations
#   results <- compare_comb_with_boin()               # 200 configurations
#   summarize_comb_comparison(results)
#
#   write.csv(results, "simFastBOIN-comb-vs-BOIN.csv", row.names = FALSE)
#
# BOIN::get.oc.comb() requires no more rows (levels of drug A) than columns
# (levels of drug B), so every matrix below has that shape.

# -----------------------------------------------------------------------------
# Grid
# -----------------------------------------------------------------------------

comb_validation_scenarios <- function() {
  list(
    "Yan et al. 2020" = matrix(c(0.02, 0.04, 0.08, 0.14,
                                 0.08, 0.25, 0.42, 0.48,
                                 0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE),
    "Lin and Yin 1" = matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
                               0.10, 0.15, 0.30, 0.45, 0.55,
                               0.15, 0.30, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE),
    "Lin and Yin 4" = matrix(c(0.30, 0.45, 0.60, 0.70, 0.80,
                               0.45, 0.55, 0.65, 0.75, 0.85,
                               0.50, 0.60, 0.70, 0.80, 0.90), nrow = 3, byrow = TRUE),
    "All safe" = matrix(c(0.01, 0.02, 0.03, 0.04,
                          0.02, 0.03, 0.04, 0.05,
                          0.03, 0.04, 0.05, 0.06), nrow = 3, byrow = TRUE),
    "Two by three" = matrix(c(0.05, 0.20, 0.35,
                              0.15, 0.30, 0.50), nrow = 2, byrow = TRUE)
  )
}

comb_validation_designs <- function(quick = FALSE) {
  designs <- list(
    list(target = 0.25, n_cohort = 16, cohort_size = 3, n_earlystop = 100,
         start_dose = c(1, 1)),
    list(target = 0.30, n_cohort = 20, cohort_size = 3, n_earlystop = 12,
         start_dose = c(1, 1)),
    list(target = 0.30, n_cohort = 15, cohort_size = 2, n_earlystop = 9,
         start_dose = c(1, 1)),
    list(target = 0.25, n_cohort = 12, cohort_size = 3, n_earlystop = 18,
         start_dose = c(2, 2)),
    list(target = 0.30, n_cohort = 30, cohort_size = 1, n_earlystop = 12,
         start_dose = c(1, 1))
  )
  if (quick) designs[2] else designs
}

comb_validation_options <- function() {
  expand.grid(
    titration = c(FALSE, TRUE),
    extrasafe = c(FALSE, TRUE),
    bound_mtd = c(FALSE, TRUE),
    KEEP.OUT.ATTRS = FALSE
  )
}

# -----------------------------------------------------------------------------
# One configuration
# -----------------------------------------------------------------------------

compare_comb_one <- function(scenario_name, p_true, design, options, n_trials,
                             seed) {

  row <- data.frame(
    scenario = scenario_name,
    target = design$target,
    n_cohort = design$n_cohort,
    cohort_size = design$cohort_size,
    n_earlystop = design$n_earlystop,
    start_dose = paste(design$start_dose, collapse = ","),
    titration = options$titration,
    extrasafe = options$extrasafe,
    bound_mtd = options$bound_mtd,
    stringsAsFactors = FALSE
  )

  reference <- tryCatch(
    BOIN::get.oc.comb(
      target = design$target, p.true = p_true,
      ncohort = design$n_cohort, cohortsize = design$cohort_size,
      n.earlystop = design$n_earlystop, startdose = design$start_dose,
      titration = options$titration, cutoff.eli = 0.95,
      extrasafe = options$extrasafe, offset = 0.05, ntrial = n_trials,
      mtd.contour = FALSE, boundMTD = options$bound_mtd, seed = seed
    ),
    error = function(e) e
  )

  if (inherits(reference, "error")) {
    row$status <- paste("BOIN error:", conditionMessage(reference))
    row$max_diff_sel <- NA_real_
    row$max_diff_npts <- NA_real_
    row$max_diff_ntox <- NA_real_
    row$diff_no_mtd <- NA_real_
    row$diff_total_n <- NA_real_
    row$diff_total_tox <- NA_real_
    row$same_pcs <- NA
    row$identical <- NA
    return(row)
  }

  ours <- sim_comb_boin(
    target = design$target, p_true = p_true,
    n_cohort = design$n_cohort, cohort_size = design$cohort_size,
    n_trials = n_trials, start_dose = design$start_dose,
    n_earlystop = design$n_earlystop, cutoff_eli = 0.95,
    extrasafe = options$extrasafe, offset = 0.05,
    titration = options$titration, bound_mtd = options$bound_mtd,
    seed = seed
  )

  row$status <- "ok"
  row$max_diff_sel <- max(abs(round(unname(ours$sel_percent), 2) -
                                unname(reference$selpercent)))
  row$max_diff_npts <- max(abs(round(unname(ours$n_pts_dose), 2) -
                                 unname(reference$npatients)))
  row$max_diff_ntox <- max(abs(round(unname(ours$n_tox_dose), 2) -
                                 unname(reference$ntox)))
  row$diff_no_mtd <- abs(ours$percent_no_mtd - reference$percentstop)
  row$diff_total_n <- abs(round(ours$total_n_pts, 1) - reference$totaln)
  row$diff_total_tox <- abs(round(ours$total_n_tox, 1) - reference$totaltox)
  row$same_pcs <- identical(paste0(round(ours$pct_true_mtd_selected, 1), "%"),
                            reference$pcs) &&
    identical(paste0(round(ours$pct_patients_true_mtd, 1), "%"),
              reference$npercent)

  # percentstop is 100 minus a sum of rounded percentages, so it is compared
  # with a floating point tolerance; everything else must be exactly equal.
  row$identical <- row$max_diff_sel == 0 && row$max_diff_npts == 0 &&
    row$max_diff_ntox == 0 && row$diff_no_mtd < 1e-8 &&
    row$diff_total_n == 0 && row$diff_total_tox == 0 && row$same_pcs
  row
}

# -----------------------------------------------------------------------------
# The whole grid
# -----------------------------------------------------------------------------

compare_comb_with_boin <- function(n_trials = 1000, seed = 6, quick = FALSE,
                                   verbose = TRUE) {

  if (!requireNamespace("BOIN", quietly = TRUE)) {
    stop("The BOIN package is needed for this comparison", call. = FALSE)
  }

  scenarios <- comb_validation_scenarios()
  if (quick) scenarios <- scenarios[1:2]
  designs <- comb_validation_designs(quick)
  options_grid <- comb_validation_options()

  total <- length(scenarios) * length(designs) * nrow(options_grid)
  if (verbose) {
    message("Comparing ", total, " configurations at ", n_trials,
            " trials each. This calls BOIN::get.oc.comb() once per configuration.")
  }

  rows <- vector("list", total)
  k <- 0
  started <- Sys.time()

  for (s in seq_along(scenarios)) {
    for (d in seq_along(designs)) {
      for (o in seq_len(nrow(options_grid))) {
        k <- k + 1
        rows[[k]] <- compare_comb_one(
          scenario_name = names(scenarios)[s],
          p_true = scenarios[[s]],
          design = designs[[d]],
          options = as.list(options_grid[o, ]),
          n_trials = n_trials,
          seed = seed
        )
        if (verbose && k %% 10 == 0) {
          message("  ", k, " of ", total, " done (",
                  format(round(difftime(Sys.time(), started, units = "secs"))), ")")
        }
      }
    }
  }

  out <- do.call(rbind, rows)
  attr(out, "n_trials") <- n_trials
  attr(out, "seed") <- seed
  attr(out, "elapsed") <- difftime(Sys.time(), started, units = "secs")
  out
}

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------

summarize_comb_comparison <- function(results) {

  cat("simFastBOIN drug-combination design against BOIN\n")
  cat("  configurations   : ", nrow(results), "\n", sep = "")
  cat("  trials each      : ", attr(results, "n_trials"), "\n", sep = "")
  cat("  seed             : ", attr(results, "seed"), "\n", sep = "")
  if (!is.null(attr(results, "elapsed"))) {
    cat("  elapsed          : ",
        format(round(attr(results, "elapsed"))), "\n", sep = "")
  }
  cat("\n")

  failed <- results$status != "ok"
  if (any(failed)) {
    cat("BOIN::get.oc.comb() could not be run for ", sum(failed),
        " configurations:\n", sep = "")
    print(unique(results$status[failed]))
    cat("\n")
  }

  usable <- results[!failed, , drop = FALSE]
  if (nrow(usable) == 0L) {
    cat("No configuration could be compared.\n")
    return(invisible(results))
  }

  quantities <- c("max_diff_sel", "max_diff_npts", "max_diff_ntox",
                  "diff_no_mtd", "diff_total_n", "diff_total_tox")
  worst <- vapply(usable[quantities], max, numeric(1))
  names(worst) <- c("MTD selection (%)", "Patients per combination",
                    "DLTs per combination", "No MTD (%)", "Total patients",
                    "Total DLTs")

  cat("Largest absolute difference over all configurations:\n")
  print(format(worst, scientific = TRUE, digits = 3), quote = FALSE)
  cat("\n")

  n_identical <- sum(usable$identical)
  status <- if (n_identical == nrow(usable)) "PASS" else "FAIL"
  cat(status, ": ", n_identical, " of ", nrow(usable),
      " configurations agree exactly.\n", sep = "")

  if (n_identical < nrow(usable)) {
    cat("\nConfigurations that do not agree:\n")
    print(usable[!usable$identical, , drop = FALSE])
    cat("\nThese are not Monte Carlo differences. The two implementations",
        "consume\nthe same random variates, so any disagreement is a defect.\n")
  }

  invisible(results)
}
