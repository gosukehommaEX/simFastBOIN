# Comparison of the waterfall design of simFastBOIN with the BOIN package
#
# Yan F, Zhang L, Zhou Y, Pan H, Liu S, Yuan Y. BOIN: An R Package for
# Designing Single-Agent and Drug-Combination Dose-Finding Trials Using
# Bayesian Optimal Interval Designs. Journal of Statistical Software 2020;
# 94(13): 1-32. doi:10.18637/jss.v094.i13
#
# The two implementations are expected to agree EXACTLY, not within a Monte
# Carlo tolerance. With boin_compatible = TRUE, sim_waterfall() draws its
# random variates in the order of BOIN::get.oc.comb(mtd.contour = TRUE),
# applies the rules of 'BOIN' 2.7.2 and calls the same bivariate isotonic
# regression, Iso::biviso(), so for a given seed the two produce the same
# trials and the same selections. BOIN::get.oc.comb() formats its matrices
# with two decimals and its totals with one; the comparison formats the
# results of simFastBOIN in the same way. A cell is judged PASS when the two
# values are equal and FAIL otherwise.
#
# The rows with the default rules of the package, which keep the data of a
# subtrial that ends without a candidate MTD, are INFO.
#
# BOIN::get.oc.comb(mtd.contour = TRUE) takes the escalation and
# de-escalation boundaries from the default values of p.saf and p.tox
# whatever values are given, whereas simFastBOIN uses the values given, so
# the comparison keeps the defaults.
#
# The scenarios are six of Table III of Zhang and Yuan (2016) (see
# inst/reproduce/reproduce_Zhang2016.R) and a matrix with steep DLT
# probabilities, in which subtrials often end without a candidate MTD.
#
# Usage, from the package root after devtools::load_all(), with BOIN
# installed:
#
#   source("inst/validation/compare-waterfall-with-BOIN.R")
#   res <- run_compare_waterfall_with_boin()
#   table(res$design, res$judgement)
#
# The judgements are written to
# inst/validation/results/compare-waterfall-with-BOIN.csv.

waterfall_validation_scenarios <- function() {
  source_file <- "inst/reproduce/reproduce_Zhang2016.R"
  if (!exists("zhang2016_scenario", mode = "function")) source(source_file)
  pub <- zhang2016_tables()
  keep <- c(1, 4, 5, 7, 9, 12)
  out <- lapply(keep, function(s) zhang2016_scenario(pub, s))
  names(out) <- paste("Scenario", keep)
  out[["Steep"]] <- matrix(c(0.05, 0.25, 0.45, 0.60, 0.70,
                             0.20, 0.40, 0.55, 0.70, 0.80,
                             0.30, 0.50, 0.65, 0.80, 0.90),
                           nrow = 3, byrow = TRUE)
  out
}

waterfall_validation_configs <- function() {
  base <- list(target = 0.30, cohort_size = 3, n_earlystop = 12,
               cutoff_eli = 0.95, titration = FALSE, extrasafe = FALSE,
               bound_mtd = FALSE, n_cohort = NULL)
  config_of <- function(name, scenario, ...) {
    utils::modifyList(c(list(name = name, scenario = scenario), base),
                      list(...), keep.null = TRUE)
  }
  list(
    config_of("plain", "Scenario 1"),
    config_of("plain", "Scenario 4"),
    config_of("plain", "Scenario 5"),
    config_of("plain", "Scenario 7"),
    config_of("plain", "Scenario 9"),
    config_of("plain", "Scenario 12"),
    config_of("plain", "Steep"),
    config_of("cohorts of 1", "Scenario 5", cohort_size = 1),
    config_of("cohorts of 2, n_earlystop 9", "Scenario 9", cohort_size = 2,
         n_earlystop = 9),
    config_of("target 0.25", "Scenario 12", target = 0.25),
    config_of("cutoff_eli 0.90", "Scenario 7", cutoff_eli = 0.90),
    config_of("n_cohort 6, 3, 3", "Scenario 9", n_cohort = c(6, 3, 3)),
    config_of("extrasafe", "Scenario 4", extrasafe = TRUE),
    config_of("extrasafe", "Steep", extrasafe = TRUE),
    config_of("bound_mtd", "Scenario 4", bound_mtd = TRUE),
    config_of("bound_mtd", "Steep", bound_mtd = TRUE),
    config_of("titration", "Scenario 1", titration = TRUE),
    config_of("titration", "Scenario 12", titration = TRUE),
    config_of("titration, extrasafe, bound_mtd, cohorts of 2", "Steep",
         cohort_size = 2, titration = TRUE, extrasafe = TRUE,
         bound_mtd = TRUE),
    config_of("titration, extrasafe, bound_mtd", "Scenario 4", titration = TRUE,
         extrasafe = TRUE, bound_mtd = TRUE)
  )
}

# sim_waterfall() warns when a matrix is not monotone, which Scenario 7 of
# Zhang and Yuan (2016) is as printed; only that warning is silenced.
waterfall_validation_quiet <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("not non-decreasing", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# The results of sim_waterfall() in the layout and with the formatting of
# BOIN::get.oc.comb(mtd.contour = TRUE).
waterfall_as_reference <- function(oc) {
  two <- function(x) {
    matrix(as.numeric(formatC(x, digits = 2, format = "f")), nrow = nrow(x))
  }
  one <- function(x) as.numeric(formatC(x, digits = 1, format = "f"))
  list(
    selpercent = round(unname(oc$sel_percent), 2),
    npatients = two(oc$n_pts_dose),
    ntox = two(oc$n_tox_dose),
    totaltox = one(sum(oc$n_tox_dose)),
    totaln = one(sum(oc$n_pts_dose)),
    npercent.contour = round(oc$pct_patients_contour, 1),
    npercent.above.contour = round(oc$pct_patients_above, 1),
    npercent.below.contour = round(oc$pct_patients_below, 1),
    pcs.contour = round(oc$pct_contour_selected, 1)
  )
}

# One value per cell of a result in the layout of the reference.
waterfall_cells <- function(x, n_rows, n_cols) {
  percent <- function(v) {
    if (is.null(v)) NA_real_ else as.numeric(sub("%", "", v, fixed = TRUE))
  }
  cell <- paste0("A", rep(seq_len(n_rows), times = n_cols), "B",
                 rep(seq_len(n_cols), each = n_rows))
  data.frame(
    quantity = c(rep("MTD selected (%)", n_rows * n_cols),
                 rep("patients treated", n_rows * n_cols),
                 rep("DLTs", n_rows * n_cols), "total DLTs",
                 "total patients", "patients on contour (%)",
                 "patients above contour (%)", "patients below contour (%)",
                 "contour selected (%)"),
    cell = c(cell, cell, cell, rep(NA_character_, 6L)),
    value = c(as.numeric(unname(x$selpercent)),
              as.numeric(unname(x$npatients)), as.numeric(unname(x$ntox)),
              x$totaltox, x$totaln, percent(x$npercent.contour),
              percent(x$npercent.above.contour),
              percent(x$npercent.below.contour), percent(x$pcs.contour)),
    stringsAsFactors = FALSE
  )
}

compare_one_waterfall_with_boin <- function(config, p_true, n_trials = 500,
                                            seed = 6) {
  n_rows <- nrow(p_true)
  n_cols <- ncol(p_true)
  run <- function(compatible) {
    waterfall_validation_quiet(sim_waterfall(
      target = config$target, p_true = p_true,
      cohort_size = config$cohort_size, n_cohort = config$n_cohort,
      n_trials = n_trials, n_earlystop = config$n_earlystop,
      cutoff_eli = config$cutoff_eli, extrasafe = config$extrasafe,
      offset = 0.05, titration = config$titration,
      bound_mtd = config$bound_mtd, boin_compatible = compatible,
      seed = seed
    ))
  }
  ours <- run(TRUE)
  default <- run(FALSE)

  reference <- BOIN::get.oc.comb(
    target = config$target, p.true = p_true,
    ncohort = ours$settings$n_cohort, cohortsize = config$cohort_size,
    n.earlystop = config$n_earlystop, titration = config$titration,
    cutoff.eli = config$cutoff_eli, extrasafe = config$extrasafe,
    offset = 0.05, ntrial = n_trials, mtd.contour = TRUE,
    boundMTD = config$bound_mtd, seed = seed
  )

  ref <- waterfall_cells(reference, n_rows, n_cols)
  blocks <- lapply(c(TRUE, FALSE), function(compatible) {
    oc <- if (compatible) ours else default
    mine <- waterfall_cells(waterfall_as_reference(oc), n_rows, n_cols)
    out <- data.frame(
      design = if (compatible) {
        "Waterfall"
      } else {
        "Waterfall (rules of the package)"
      },
      scenario = config$scenario, configuration = config$name,
      quantity = ref$quantity, cell = ref$cell, BOIN = ref$value,
      simFastBOIN = mine$value, stringsAsFactors = FALSE
    )
    out$difference <- out$simFastBOIN - out$BOIN
    out$judgement <- if (compatible) {
      ifelse(!is.na(out$difference) & abs(out$difference) < 1e-9, "PASS",
             "FAIL")
    } else {
      "INFO"
    }
    out
  })
  do.call(rbind, blocks)
}

run_compare_waterfall_with_boin <- function(n_trials = 500, seed = 6,
                                            out_dir =
                                              "inst/validation/results") {
  # The functions of the package under development must be loaded.
  if (!exists("sim_waterfall", mode = "function")) {
    stop("run devtools::load_all() from the package root first", call. = FALSE)
  }
  if (!requireNamespace("BOIN", quietly = TRUE)) {
    stop("the 'BOIN' package is needed", call. = FALSE)
  }
  scenarios <- waterfall_validation_scenarios()
  configs <- waterfall_validation_configs()
  runs <- lapply(configs, function(config) {
    compare_one_waterfall_with_boin(config, scenarios[[config$scenario]],
                                    n_trials = n_trials, seed = seed)
  })
  res <- do.call(rbind, runs)
  rownames(res) <- NULL
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, file.path(out_dir, "compare-waterfall-with-BOIN.csv"),
                   row.names = FALSE)
  res
}
