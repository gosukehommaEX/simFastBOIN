# Comprehensive comparison of simFastBOIN against the BOIN package
# -----------------------------------------------------------------------------
#
# The two implementations are expected to agree EXACTLY, not within a Monte
# Carlo tolerance. simFastBOIN draws one uniform variate per patient in
# enrollment order and applies the decision rules in the same order as
# BOIN::get.oc(), so for a given seed the two produce the same trials, patient
# by patient. Any non-zero difference in this report is therefore a defect to
# investigate, not sampling noise.
#
# Usage
#
#   library(simFastBOIN)
#   source(system.file("validation", "compare-with-BOIN.R", package = "simFastBOIN"))
#
#   results <- compare_with_boin(quick = TRUE)    # about half a minute
#   results <- compare_with_boin()                # the full grid, a few minutes
#   summarize_comparison(results)
#   check_inert_options()
#
#   write.csv(results, "simFastBOIN-vs-BOIN.csv", row.names = FALSE)
#
# Note on the scenarios: none of them place two doses at exactly the target
# rate. BOIN::get.oc() builds its overdosing summary with
# `if (which(p.true == target) == ndose - 1)`, which is not written for a
# condition of length greater than one, so such a scenario stops with an error
# inside BOIN before any comparison can be made. simFastBOIN handles those
# scenarios, but they cannot be cross-checked here.

# -----------------------------------------------------------------------------
# Grid
# -----------------------------------------------------------------------------

validation_scenarios <- function() {
  list(
    "MTD at dose 3"   = c(0.05, 0.15, 0.25, 0.45, 0.60),
    "MTD at dose 1"   = c(0.28, 0.42, 0.55, 0.68, 0.80),
    "MTD at dose 5"   = c(0.02, 0.05, 0.09, 0.16, 0.28),
    "All doses toxic" = c(0.45, 0.55, 0.65, 0.75, 0.85),
    "All doses safe"  = c(0.01, 0.03, 0.06, 0.10, 0.16)
  )
}

validation_designs <- function(quick = FALSE) {
  designs <- list(
    list(target = 0.30, n_cohort = 20, cohort_size = 3, n_earlystop = 100, start_dose = 1),
    list(target = 0.30, n_cohort = 20, cohort_size = 3, n_earlystop = 18,  start_dose = 1),
    list(target = 0.25, n_cohort = 12, cohort_size = 3, n_earlystop = 12,  start_dose = 2),
    list(target = 0.20, n_cohort = 15, cohort_size = 2, n_earlystop = 10,  start_dose = 1),
    list(target = 0.40, n_cohort = 24, cohort_size = 1, n_earlystop = 12,  start_dose = 1)
  )
  if (quick) designs[2] else designs
}

validation_options <- function() {
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

compare_one <- function(scenario_name, p_true, design, options, n_trials, seed) {

  row <- data.frame(
    scenario = scenario_name,
    target = design$target,
    n_cohort = design$n_cohort,
    cohort_size = design$cohort_size,
    n_earlystop = design$n_earlystop,
    start_dose = design$start_dose,
    titration = options$titration,
    extrasafe = options$extrasafe,
    bound_mtd = options$bound_mtd,
    stringsAsFactors = FALSE
  )

  reference <- tryCatch(
    BOIN::get.oc(
      target = design$target, p.true = p_true,
      ncohort = design$n_cohort, cohortsize = design$cohort_size,
      n.earlystop = design$n_earlystop, startdose = design$start_dose,
      titration = options$titration, cutoff.eli = 0.95,
      extrasafe = options$extrasafe, offset = 0.05,
      boundMTD = options$bound_mtd, ntrial = n_trials, seed = seed
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
    row$identical <- NA
    return(row)
  }

  ours <- sim_boin(
    target = design$target, p_true = p_true,
    n_cohort = design$n_cohort, cohort_size = design$cohort_size,
    n_trials = n_trials, start_dose = design$start_dose,
    n_earlystop = design$n_earlystop, cutoff_eli = 0.95,
    extrasafe = options$extrasafe, offset = 0.05,
    titration = options$titration, bound_mtd = options$bound_mtd,
    seed = seed
  )

  row$status <- "ok"
  row$max_diff_sel <- max(abs(unname(ours$sel_percent) - reference$selpercent))
  row$max_diff_npts <- max(abs(unname(ours$n_pts_dose) - reference$npatients))
  row$max_diff_ntox <- max(abs(unname(ours$n_tox_dose) - reference$ntox))
  row$diff_no_mtd <- abs(ours$percent_no_mtd - reference$percentstop)
  row$diff_total_n <- abs(ours$total_n_pts - reference$totaln)
  row$diff_total_tox <- abs(ours$total_n_tox - reference$totaltox)

  # A floating point tolerance, not a Monte Carlo one: the quantities are means
  # of identical integers, accumulated in a different order by each package.
  row$identical <- max(row$max_diff_sel, row$max_diff_npts, row$max_diff_ntox,
                       row$diff_no_mtd, row$diff_total_n, row$diff_total_tox) < 1e-9
  row
}

# -----------------------------------------------------------------------------
# The whole grid
# -----------------------------------------------------------------------------

compare_with_boin <- function(n_trials = 2000, seed = 6, quick = FALSE,
                              verbose = TRUE) {

  if (!requireNamespace("BOIN", quietly = TRUE)) {
    stop("The BOIN package is needed for this comparison", call. = FALSE)
  }

  scenarios <- validation_scenarios()
  if (quick) scenarios <- scenarios[1:2]
  designs <- validation_designs(quick)
  options_grid <- validation_options()

  total <- length(scenarios) * length(designs) * nrow(options_grid)
  if (verbose) {
    message("Comparing ", total, " configurations at ", n_trials,
            " trials each. This calls BOIN::get.oc() once per configuration.")
  }

  rows <- vector("list", total)
  k <- 0
  started <- Sys.time()

  for (s in seq_along(scenarios)) {
    for (d in seq_along(designs)) {
      for (o in seq_len(nrow(options_grid))) {
        k <- k + 1
        rows[[k]] <- compare_one(
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

summarize_comparison <- function(results) {

  cat("simFastBOIN against BOIN\n")
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
    cat("BOIN::get.oc() could not be run for ", sum(failed),
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
  names(worst) <- c("MTD selection (%)", "Patients per dose", "DLTs per dose",
                    "No MTD (%)", "Total patients", "Total DLTs")

  cat("Largest absolute difference over all configurations:\n")
  print(format(worst, scientific = TRUE, digits = 3), quote = FALSE)
  cat("\n")

  n_identical <- sum(usable$identical)
  cat(n_identical, " of ", nrow(usable),
      " configurations agree to within floating point.\n", sep = "")

  if (n_identical < nrow(usable)) {
    cat("\nConfigurations that do not agree:\n")
    print(usable[!usable$identical, , drop = FALSE])
    cat("\nThese are not Monte Carlo differences. The two implementations",
        "consume\nthe same random variates, so any disagreement is a defect.\n")
  } else {
    cat("The two implementations produce the same trials, patient by patient.\n")
  }

  invisible(results)
}

# -----------------------------------------------------------------------------
# Options with no counterpart in BOIN must be inert when switched off, and must
# leave the trials untouched when they only affect the MTD selection.
# -----------------------------------------------------------------------------

check_inert_options <- function(n_trials = 2000, seed = 6) {

  p_true <- c(0.05, 0.15, 0.25, 0.45, 0.60)
  base <- list(target = 0.30, p_true = p_true, n_cohort = 20, cohort_size = 3,
               n_trials = n_trials, keep_trials = TRUE, seed = seed)

  plain <- do.call(sim_boin, base)

  # The one-of-three modification does not apply at a target of 0.30, so the
  # simulation must be bit for bit the same.
  stayed <- do.call(sim_boin, c(base, list(stay_on_1_of_3 = TRUE)))
  inert_ok <- identical(plain$sel_percent, stayed$sel_percent) &&
    identical(plain$trials$n_pts, stayed$trials$n_pts)

  # Capping the MTD estimate changes the selection only, never the trials.
  capped <- do.call(sim_boin, c(base, list(mtd_max_estimate = 0.30)))
  trials_ok <- identical(plain$trials$n_pts, capped$trials$n_pts) &&
    identical(plain$trials$n_tox, capped$trials$n_tox)

  # The overdose cutoff is a summary, so it cannot touch the trials either.
  recut <- do.call(sim_boin, c(base, list(overdose_cutoff = 0.33)))
  summary_ok <- identical(plain$sel_percent, recut$sel_percent) &&
    identical(plain$trials$n_pts, recut$trials$n_pts)

  cat("Options with no counterpart in BOIN\n")
  cat("  stay_on_1_of_3 inert at target 0.30 : ", inert_ok, "\n", sep = "")
  cat("  mtd_max_estimate leaves trials alone: ", trials_ok, "\n", sep = "")
  cat("  overdose_cutoff is a summary only   : ", summary_ok, "\n", sep = "")

  invisible(c(stay_on_1_of_3 = inert_ok, mtd_max_estimate = trials_ok,
              overdose_cutoff = summary_ok))
}
