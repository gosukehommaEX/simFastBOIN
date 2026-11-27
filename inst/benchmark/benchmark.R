# Computation time of simFastBOIN and of the reference implementations
# -----------------------------------------------------------------------------
#
# Every design that was compared with a reference implementation, in another
# R package or in the code of the authors of the design, is run by both with
# the same settings and the same number of trials: the scenario and the
# design of the comparison in inst/validation/ or inst/reproduce/ (or of the
# README for the BOIN design), with the number of trials used there for the
# reference. The other designs are run by simFastBOIN alone, with the
# settings of the example of their vignette and 10,000 trials.
#
# A call that is fast is repeated until the repetitions take at least
# min_time seconds, and the elapsed (wall-clock) time per call is recorded.
# Every measurement is taken reps times, alternating between the two
# implementations. Printed output and warnings of the calls are discarded.
#
# The package must be installed, not loaded with devtools::load_all():
# load_all() compiles the C++ code without optimization, which makes it
# several times slower. Usage, from the package root, in a new R session:
#
#   devtools::install()
#   # restart R
#   library(simFastBOIN)
#   source("inst/benchmark/benchmark.R")
#   res <- run_benchmark(boinetc_zip = "path/to/BOINETC-main.zip",
#                        takeda2023_zip = "path/to/Takeda_2023_SI.zip")
#
# The reference packages are 'BOIN', 'boinet', 'bfboinet',
# 'UnifiedDoseFinding' and 'escalation'; a design whose package is not
# installed is recorded with the status "not installed". The code of the
# authors of BOIN-ETC (the zip file of https://github.com/kakizume/BOINETC)
# and the R program of the supporting information of Takeda et al. (2023)
# (its zip file) are unpacked into a temporary folder. The code of BOIN-ETC
# is read by boinetc_env() of inst/validation/compare-with-BOINETC.R, which
# needs 'Iso'; the program of Takeda et al. (2023) is run as published, for
# one setting and the given number of trials. Without the zip files the two
# designs are recorded with the status "not available". The times are
# written, with the date, the operating system, the processor and the
# versions of R and of the packages, to inst/benchmark/results/benchmark.csv,
# after every design, so that an interrupted run keeps the designs already
# measured. The reference packages take most of the time (several minutes
# for each of the slowest designs and each repetition).
# -----------------------------------------------------------------------------

# The functions of a reproduction script of the installed package, which
# define the scenarios of the articles.
benchmark_source <- function(file) {
  env <- new.env(parent = globalenv())
  sys.source(system.file("reproduce", file, package = "simFastBOIN",
                         mustWork = TRUE), envir = env)
  env
}

# The designs. Each has a label, the setting, the number of trials, the call
# of simFastBOIN and, when there is one, the call of the reference package.
benchmark_cases <- function(boinetc_prog = NULL, takeda2023_prog = NULL) {
  t18 <- benchmark_source("reproduce_Takeda2018.R")$takeda2018_scenarios()[[1]]
  t22 <- benchmark_source("reproduce_Takeda2022.R")$takeda2022_scenarios(
    "table")[[1]]
  t23 <- benchmark_source("reproduce_Takeda2023.R")$takeda2023_scenarios()[[1]]
  mu <- benchmark_source("reproduce_Mu2019.R")$mu2019_scenarios()[[1]]
  zhang <- benchmark_source("reproduce_Zhang2016.R")
  wf <- zhang$zhang2016_scenario(zhang$zhang2016_tables(
    system.file("reproduce", "zhang2016-tables.csv",
                package = "simFastBOIN")), 9)
  wf_cohorts <- sim_waterfall(target = 0.30, p_true = wf, cohort_size = 3,
                              n_trials = 1, boin_compatible = TRUE,
                              seed = 6)$settings$n_cohort
  comb <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
                   0.10, 0.15, 0.30, 0.45, 0.55,
                   0.15, 0.30, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE)
  boin12_tox <- c(0.05, 0.12, 0.27, 0.35, 0.50)
  boin12_eff <- c(0.20, 0.35, 0.36, 0.37, 0.38)
  bf_tox <- c(0.03, 0.05, 0.20, 0.22, 0.45)
  bf_eff <- c(0.05, 0.10, 0.50, 0.68, 0.70)
  score <- c(0, 0.5, 1, 1.5)
  comb12_tox <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
  comb12_eff <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
  etc_tox <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
  etc_eff <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
  safe_curve <- c(0.01, 0.02, 0.05, 0.12, 0.30)
  tite_curve <- c(0.05, 0.15, 0.30, 0.45, 0.60)
  bf_curve <- c(0.04, 0.12, 0.25, 0.43, 0.63)
  bf_resp <- c(0.10, 0.20, 0.30, 0.45, 0.58)
  n10k <- 10000
  etc_code <- benchmark_boinetc(boinetc_prog)
  si_run <- benchmark_takeda2023(takeda2023_prog, n_trials = 1000)

  # A reference given as code (available = TRUE or FALSE) rather than as a
  # package (available = NULL) has no version.
  case <- function(design, setting, n_trials, ours, fun, reference = NULL,
                   package = NA_character_, reference_fun = NA_character_,
                   available = NULL) {
    list(design = design, setting = setting, n_trials = n_trials,
         ours = ours, fun = fun, reference = reference, package = package,
         reference_fun = reference_fun, available = available)
  }

  list(
    case("BOIN",
         "Five doses, 20 cohorts of three, n_earlystop = 18 (README)",
         10000,
         function(n) {
           sim_boin(target = 0.30, p_true = c(0.05, 0.15, 0.25, 0.45, 0.60),
                    n_cohort = 20, cohort_size = 3, n_earlystop = 18,
                    n_trials = n, seed = 6)
         }, "sim_boin()",
         function(n) {
           BOIN::get.oc(target = 0.30,
                        p.true = c(0.05, 0.15, 0.25, 0.45, 0.60),
                        ncohort = 20, cohortsize = 3, n.earlystop = 18,
                        ntrial = n, seed = 6)
         }, "BOIN", "BOIN::get.oc()"),
    case("BOIN combination",
         paste("Lin and Yin 1 (3 x 5), 20 cohorts of three,",
               "n_earlystop = 12 (compare-comb-with-BOIN.R)"),
         1000,
         function(n) {
           sim_comb_boin(target = 0.30, p_true = comb, n_cohort = 20,
                         cohort_size = 3, n_trials = n, start_dose = c(1, 1),
                         n_earlystop = 12, cutoff_eli = 0.95,
                         extrasafe = FALSE, offset = 0.05, titration = FALSE,
                         bound_mtd = FALSE, seed = 6)
         }, "sim_comb_boin()",
         function(n) {
           BOIN::get.oc.comb(target = 0.30, p.true = comb, ncohort = 20,
                             cohortsize = 3, n.earlystop = 12,
                             startdose = c(1, 1), titration = FALSE,
                             cutoff.eli = 0.95, extrasafe = FALSE,
                             offset = 0.05, ntrial = n, mtd.contour = FALSE,
                             boundMTD = FALSE, seed = 6)
         }, "BOIN", "BOIN::get.oc.comb()"),
    case("Waterfall",
         paste("Scenario 9 of Zhang and Yuan (2016) (3 x 5), cohorts of",
               "three, boin_compatible = TRUE",
               "(compare-waterfall-with-BOIN.R)"),
         500,
         function(n) {
           sim_waterfall(target = 0.30, p_true = wf, cohort_size = 3,
                         n_trials = n, n_earlystop = 12, cutoff_eli = 0.95,
                         extrasafe = FALSE, offset = 0.05, titration = FALSE,
                         bound_mtd = FALSE, boin_compatible = TRUE, seed = 6)
         }, "sim_waterfall()",
         function(n) {
           BOIN::get.oc.comb(target = 0.30, p.true = wf, ncohort = wf_cohorts,
                             cohortsize = 3, n.earlystop = 12,
                             titration = FALSE, cutoff.eli = 0.95,
                             extrasafe = FALSE, offset = 0.05, ntrial = n,
                             mtd.contour = TRUE, boundMTD = FALSE, seed = 6)
         }, "BOIN", "BOIN::get.oc.comb(mtd.contour = TRUE)"),
    case("BOIN-ET",
         paste("Scenario 1 of Takeda et al. (2018), 12 cohorts of three",
               "(compare-with-boinet.R)"),
         5000,
         function(n) {
           sim_boin_et(target = 0.30, target_resp = 0.60, p_true = t18$p_true,
                       p_resp = t18$p_resp, n_cohort = 12, cohort_size = 3,
                       n_trials = n, resp_cor = 0.2, est_method = "observed",
                       obd_method = "max_resp", boinet_compatible = TRUE,
                       seed = 1)
         }, "sim_boin_et()",
         function(n) {
           boinet::boinet(n.dose = length(t18$p_true), start.dose = 1,
                          size.cohort = 3, n.cohort = 12,
                          toxprob = t18$p_true, effprob = t18$p_resp,
                          phi = 0.30, delta = 0.60, tau.T = 30, tau.E = 30,
                          te.corr = 0.2, accrual = 10,
                          estpt.method = "obs.prob",
                          obd.method = "max.effprob", n.sim = n,
                          seed.sim = 1)
         }, "boinet", "boinet::boinet()"),
    case("TITE-BOIN-ET",
         paste("Scenario 1 of Takeda et al. (2018), windows of 30 and 45",
               "days, a patient every 10 days (compare-with-boinet.R)"),
         2000,
         function(n) {
           sim_tite_boin_et(target = 0.30, target_resp = 0.60,
                            p_true = t18$p_true, p_resp = t18$p_resp,
                            n_cohort = 12, cohort_size = 3, window = 30,
                            resp_window = 45, accrual_rate = 1 / 10,
                            accrual = "uniform", event_time = "weibull",
                            n_trials = n, resp_cor = 0.2,
                            est_method = "observed", obd_method = "max_resp",
                            boinet_compatible = TRUE, seed = 1)
         }, "sim_tite_boin_et()",
         function(n) {
           boinet::tite.boinet(n.dose = length(t18$p_true), start.dose = 1,
                               size.cohort = 3, n.cohort = 12,
                               toxprob = t18$p_true, effprob = t18$p_resp,
                               phi = 0.30, delta = 0.60, tau.T = 30,
                               tau.E = 45, te.corr = 0.2,
                               gen.event.time = "weibull", accrual = 10,
                               gen.enroll.time = "uniform",
                               estpt.method = "obs.prob",
                               obd.method = "max.effprob", n.sim = n,
                               seed.sim = 1)
         }, "boinet", "boinet::tite.boinet()"),
    case("gBOIN-ET",
         paste("Scenario 1 of Takeda et al. (2022), probabilities of Table",
               "S1, 10 cohorts of three (compare-with-gboinet.R)"),
         5000,
         function(n) {
           sim_gboin_et(target = 0.47 / 1.5, target_resp = 0.95 / 1.5,
                        p_grade = t22$p_grade, p_resp_grade = t22$p_resp_grade,
                        score = score, score_resp = c(0, 1, 1.5),
                        n_cohort = 10, cohort_size = 3, n_trials = n,
                        resp_cor = 0.2, obd_method = "max_resp",
                        boinet_compatible = TRUE, seed = 1)
         }, "sim_gboin_et()",
         function(n) {
           boinet::gboinet(n.dose = ncol(t22$p_grade), start.dose = 1,
                           size.cohort = 3, n.cohort = 10,
                           toxprob = t22$p_grade, effprob = t22$p_resp_grade,
                           sev.weight = score, res.weight = c(0, 1, 1.5),
                           phi = 0.47 / 1.5, delta = 0.95 / 1.5, tau.T = 30,
                           tau.E = 30, te.corr = 0.2, accrual = 10,
                           estpt.method = "obs.prob",
                           obd.method = "max.effprob", n.sim = n,
                           seed.sim = 1)
         }, "boinet", "boinet::gboinet()"),
    case("TITE-gBOIN-ET",
         paste("Scenario 1 of Takeda et al. (2023), windows of 30 and 60",
               "days, 12 cohorts of three (compare-with-tite-gboinet.R)"),
         5000,
         function(n) {
           sim_tite_gboin_et(target = 0.47 / 1.5, target_resp = 1.75 / 3,
                             p_grade = t23$p_grade,
                             p_resp_grade = t23$p_resp_grade, score = score,
                             score_resp = c(0, 0.25, 1, 3), n_cohort = 12,
                             cohort_size = 3, window = 30,
                             accrual_rate = 1 / 10, resp_window = 60,
                             resp_cor = 0.2, n_trials = n,
                             obd_method = "max_resp",
                             boinet_compatible = TRUE, seed = 1)
         }, "sim_tite_gboin_et()",
         function(n) {
           boinet::tite.gboinet(n.dose = ncol(t23$p_grade), start.dose = 1,
                                size.cohort = 3, n.cohort = 12,
                                toxprob = t23$p_grade,
                                effprob = t23$p_resp_grade,
                                sev.weight = score,
                                res.weight = c(0, 0.25, 1, 3),
                                phi = 0.47 / 1.5, delta = 1.75 / 3,
                                tau.T = 30, tau.E = 60, te.corr = 0.2,
                                accrual = 10, estpt.method = "obs.prob",
                                obd.method = "max.effprob", n.sim = n,
                                seed.sim = 1)
         }, "boinet", "boinet::tite.gboinet()"),
    case("TITE-gBOIN-ET, program of the article",
         paste("Scenario 1 of Takeda et al. (2023), windows of 30 and 60",
               "days, a patient every 10 days, 12 cohorts of three; the",
               "program of the supporting information as published",
               "(reproduce_Takeda2023.R)"),
         1000,
         function(n) {
           sim_tite_gboin_et(target = 0.47 / 1.5, target_resp = 1.75 / 3,
                             p_grade = t23$p_grade,
                             p_resp_grade = t23$p_resp_grade, score = score,
                             score_resp = c(0, 0.25, 1, 3), n_cohort = 12,
                             cohort_size = 3, window = 30,
                             accrual_rate = 1 / 10, resp_window = 60,
                             resp_cor = 0.2, n_trials = n,
                             obd_method = "max_resp",
                             boinet_compatible = TRUE, seed = 1)
         }, "sim_tite_gboin_et()",
         si_run, "program of Takeda et al. (2023)",
         "simulation1_TITE-gBOIN-ET.R of Takeda et al. (2023)",
         available = !is.null(si_run)),
    case("BF-BOIN-ET",
         paste("Example of 'bfboinet', 10 cohorts of three, windows of one",
               "month, three patients a month (compare-with-bfboinet.R)"),
         1000,
         function(n) {
           sim_bf_boin_et(target = 0.30, target_resp = 0.25, p_true = bf_tox,
                          p_resp = bf_eff, n_cohort = 10, cohort_size = 3,
                          window = 1, resp_window = 1, accrual_rate = 3,
                          accrual = "uniform", event_time = "weibull",
                          suspend = 0, n_trials = n, seed = 1)
         }, "sim_bf_boin_et()",
         function(n) {
           bfboinet::get.oc.backboinet(
             target_T = 0.30, toxprob = bf_tox, target_E = 0.25,
             effprob = bf_eff, n.dose = length(bf_tox), startdose = 1,
             ncohort = 10, cohortsize = 3, pT.saf = 0.6 * 0.30,
             pT.tox = 1.4 * 0.30, pE.saf = 0.6 * 0.25, alpha.T1 = 0.5,
             alpha.E1 = 0.5, tau.T = 1, tau.E = 1, te.corr = 0.2,
             gen.event.time = "weibull", accrual = 3,
             gen.enroll.time = "uniform", n.elimination = 6,
             stopping.npts = 12, suspend = 0, stopping.prob.T = 0.95,
             stopping.prob.E = 0.90, ppsi01 = 0, ppsi00 = 40, ppsi11 = 60,
             ppsi10 = 100, n.sim = n, seed.sim = 1
           )
         }, "bfboinet", "bfboinet::get.oc.backboinet()"),
    case("gBOIN",
         paste("Scenario 1 of Mu et al. (2019), 10 cohorts of three,",
               "elim_rule = \"count\" (compare-with-UnifiedDoseFinding.R)"),
         4000,
         function(n) {
           sim_gboin(target = 0.47 / 1.5, p_grade = mu, score = score / 1.5,
                     n_cohort = 10, cohort_size = 3, n_trials = n,
                     n_earlystop = 100, elim_rule = "count", seed = 2019)
         }, "sim_gboin()",
         function(n) {
           UnifiedDoseFinding::get_oc_QuasiBOIN(
             target = 0.47 / 1.5, p.true = t(mu), score = score / 1.5,
             ncohort = 10, cohortsize = 3, n.earlystop = 100, startdose = 1,
             ntrial = n, seed = 2020
           )
         }, "UnifiedDoseFinding", "UnifiedDoseFinding::get_oc_QuasiBOIN()"),
    case("BOIN12",
         paste("Scenario 1 of Lin et al. (2020), 12 cohorts of three",
               "(compare-with-escalation.R)"),
         2000,
         function(n) {
           sim_boin12(target = 0.35, target_resp = 0.25, p_true = boin12_tox,
                      p_resp = boin12_eff, n_cohort = 12, cohort_size = 3,
                      n_trials = n, est_method = "tail_probability",
                      seed = 2027)
         }, "sim_boin12()",
         function(n) {
           set.seed(2026)
           model <- escalation::select_boin12_obd(
             escalation::stop_at_n(
               escalation::get_boin12(num_doses = length(boin12_tox),
                                      phi_t = 0.35, phi_e = 0.25, u2 = 40,
                                      u3 = 60, n_star = 6),
               n = 36
             )
           )
           escalation::simulate_trials(model, num_sims = n,
                                       true_prob_tox = boin12_tox,
                                       true_prob_eff = boin12_eff)
         }, "escalation", "escalation::simulate_trials()"),
    case("iBOIN",
         paste("Scenario 3 of Zhou et al. (2021), prior_ess = 3, 10 cohorts",
               "of three (vignette)"),
         n10k,
         function(n) {
           sim_iboin(target = 0.30, p_true = c(0.08, 0.15, 0.31, 0.45, 0.55),
                     skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
                     prior_ess = 3, n_cohort = 10, cohort_size = 3,
                     n_earlystop = 100, n_trials = n, seed = 123)
         }, "sim_iboin()"),
    case("Accelerated BOIN",
         paste("titration_doses = 2, n_earlystop = 6, eight cohorts of",
               "three (vignette)"),
         n10k,
         function(n) {
           sim_boin(target = 0.30, p_true = safe_curve, n_cohort = 8,
                    cohort_size = 3, n_trials = n, titration = TRUE,
                    titration_doses = 2, n_earlystop = 6, seed = 1)
         }, "sim_boin()"),
    case("TITE-BOIN, imputation",
         paste("Window of three months, two patients a month, 10 cohorts",
               "of three (vignette)"),
         n10k,
         function(n) {
           sim_tite_boin(target = 0.30, p_true = tite_curve, n_cohort = 10,
                         cohort_size = 3, window = 3, accrual_rate = 2,
                         method = "imputation", n_trials = n, seed = 123)
         }, "sim_tite_boin()"),
    case("TITE-BOIN, ESS",
         paste("Window of three months, two patients a month, 10 cohorts",
               "of three (vignette)"),
         n10k,
         function(n) {
           sim_tite_boin(target = 0.30, p_true = tite_curve, n_cohort = 10,
                         cohort_size = 3, window = 3, accrual_rate = 2,
                         method = "ess", n_trials = n, seed = 123)
         }, "sim_tite_boin()"),
    case("BF-BOIN",
         paste("Setting of Zhao et al. (2024), window of one month, three",
               "patients a month (vignette)"),
         n10k,
         function(n) {
           sim_bf_boin(target = 0.25, p_true = bf_curve, p_resp = bf_resp,
                       n_cohort = 10, cohort_size = 3, window = 1,
                       accrual_rate = 3, n_earlystop = 9,
                       stay_on_1_of_3 = TRUE, accrual = "uniform",
                       no_slot = "leave", n_trials = n, seed = 123)
         }, "sim_bf_boin()"),
    case("BE-BOIN",
         paste("Scenario of Zhao et al. (2024), window of three months, two",
               "patients a month (vignette)"),
         n10k,
         function(n) {
           sim_be_boin(target = 0.25, p_true = bf_curve, p_resp = bf_resp,
                       n_cohort = 10, cohort_size = 3, window = 3,
                       accrual_rate = 2, n_earlystop = 9, n_trials = n,
                       seed = 123)
         }, "sim_be_boin()"),
    case("TITE-BOIN12",
         paste("Scenario 1 of Lin et al. (2020), windows of 45 and 60",
               "days, three patients every 30 days, 12 cohorts of three",
               "(vignette)"),
         n10k,
         function(n) {
           sim_tite_boin12(target = 0.35, target_resp = 0.25,
                           p_true = boin12_tox, p_resp = boin12_eff,
                           n_cohort = 12, cohort_size = 3, window = 45,
                           resp_window = 60, accrual_rate = 3 / 30,
                           n_trials = n, seed = 123)
         }, "sim_tite_boin12()"),
    case("U-BOIN",
         paste("Scenario 1 of Simulation A of Zhou, Lee and Yuan (2019),",
               "18 cohorts of three (vignette)"),
         n10k,
         function(n) {
           sim_uboin(tox_upper = 0.30, resp_lower = 0.20,
                     p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
                     p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65), n_cohort = 18,
                     cohort_size = 3, resp_assoc = 0.2, n_trials = n,
                     seed = 123)
         }, "sim_uboin()"),
    case("Comb-BOIN12",
         paste("Scenario 16 of Lu et al. (2025) (2 x 3), 14 cohorts of",
               "three (vignette)"),
         n10k,
         function(n) {
           sim_comb_boin12(target = 0.33, target_resp = 0.2,
                           p_true = comb12_tox, p_resp = comb12_eff,
                           n_cohort = 14, cohort_size = 3,
                           psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11,
                           benchmark_weight = 0.6, extrasafe = TRUE,
                           offset = 0.10, n_trials = n, seed = 123)
         }, "sim_comb_boin12()"),
    case("TITE-Comb-BOIN12",
         paste("Scenario 16 of Lu et al. (2025), windows of one and two",
               "months, two patients a month (vignette)"),
         n10k,
         function(n) {
           sim_tite_comb_boin12(target = 0.33, target_resp = 0.2,
                                p_true = comb12_tox, p_resp = comb12_eff,
                                n_cohort = 14, cohort_size = 3,
                                psi00 = 100 * 0.11 / 1.11,
                                psi11 = 100 / 1.11, benchmark_weight = 0.6,
                                extrasafe = TRUE, offset = 0.10, window = 1,
                                resp_window = 2, accrual_rate = 2,
                                event_time = "uniform", n_trials = n,
                                seed = 123)
         }, "sim_tite_comb_boin12()"),
    case("BOIN-ETC",
         paste("Scenario 1 of Kakizume et al. (2024) (2 x 3), BOIN-ETC1,",
               "eight cohorts of three, n_earlystop = 6, the rules of the",
               "code of the authors (compare-with-BOINETC.R)"),
         1000,
         function(n) {
           sim_boin_etc(target = 0.30, target_resp = 0.60, p_true = etc_tox,
                        p_resp = etc_eff, n_cohort = 8, cohort_size = 3,
                        n_trials = n, subtrials = "etc1", n_earlystop = 6,
                        etc1_start = "below", rule4 = "authors_code",
                        seed = 123)
         }, "sim_boin_etc()",
         if (is.null(etc_code)) NULL else function(n) {
           etc_code$get.oc.BOINETC(
             ncohort = 8, cohortsize = 3, n.earlystop = 6, ntrial = n,
             phi = 0.30, delta = 0.60, lambda1 = 0.14, lambda2 = 0.35,
             eta1 = 0.48, tdose = etc_code$list.tdose[[1]], u.min = NA,
             n.scr = 1, pt.true.mat = etc_code$list.p.truetox[[1]],
             pe.true.mat = etc_code$list.p.trueeff[[1]],
             filename = "benchmark"
           )
         }, "code of Kakizume et al. (2024)",
         "get.oc.BOINETC() of Kakizume et al. (2024)",
         available = !is.null(etc_code))
  )
}

# The code of the authors of BOIN-ETC in an environment, read by
# boinetc_env() of inst/validation/compare-with-BOINETC.R; NULL without it.
benchmark_boinetc <- function(prog_dir) {
  if (is.null(prog_dir)) return(NULL)
  v <- new.env(parent = globalenv())
  sys.source(system.file("validation", "compare-with-BOINETC.R",
                         package = "simFastBOIN", mustWork = TRUE),
             envir = v)
  work <- tempfile("boinetc")
  dir.create(work)
  v$boinetc_env(prog_dir, work)
}

# A function that runs the R program of the supporting information of
# Takeda et al. (2023) for Scenario 1, an efficacy window of 60 days and a
# patient every 10 days, with n_trials trials: a copy of the program with
# these four lines changed, in the folder of the program; NULL without it.
benchmark_takeda2023 <- function(prog_dir, n_trials) {
  if (is.null(prog_dir)) return(NULL)
  script <- file.path(prog_dir, "simulation1", "simulation1_TITE-gBOIN-ET.R")
  src <- readLines(script, warn = FALSE)
  change <- list(
    c("ls1   <- length(unique(toxprob$Scenario))", "ls1   <- 1"),
    c("tau.E <- c(45,60,90)", "tau.E <- 60"),
    c("accrual <- c(5,10)", "accrual <- 10"),
    c("lss <- 10000", paste0("lss <- ", n_trials))
  )
  for (ch in change) {
    hit <- which(src == ch[1])
    if (length(hit) != 1L) {
      stop("not found once in ", script, ": ", ch[1], call. = FALSE)
    }
    src[hit] <- ch[2]
  }
  copy <- file.path(dirname(script), "benchmark_one_setting.R")
  writeLines(src, copy)
  function(n) {
    stopifnot(n == n_trials)
    source(copy, local = new.env(parent = globalenv()), chdir = TRUE)
  }
}

# Unpack a zip file into a temporary folder and return the folder that holds
# the given subfolder; NULL without a zip file.
benchmark_unzip <- function(zip, subfolder) {
  if (is.null(zip)) return(NULL)
  if (!file.exists(zip)) stop("not found: ", zip, call. = FALSE)
  exdir <- tempfile("benchmark")
  utils::unzip(zip, exdir = exdir)
  hit <- list.dirs(exdir, recursive = TRUE)
  hit <- hit[basename(hit) == basename(subfolder) &
               basename(dirname(hit)) == basename(dirname(subfolder))]
  if (length(hit) != 1L) {
    stop("the folder ", subfolder, " is not in ", zip, call. = FALSE)
  }
  hit
}

# Evaluate a call, discarding its printed output and its warnings.
benchmark_quiet <- function(expr) {
  utils::capture.output(out <- suppressWarnings(expr))
  invisible(out)
}

# The elapsed time of one call, in seconds: the call is repeated until the
# repetitions take at least min_time seconds.
benchmark_time <- function(run, n_trials, min_time) {
  invisible(gc(verbose = FALSE))
  runs <- 0L
  start <- proc.time()[["elapsed"]]
  repeat {
    benchmark_quiet(run(n_trials))
    runs <- runs + 1L
    total <- proc.time()[["elapsed"]] - start
    if (total >= min_time) break
  }
  c(seconds = total / runs, runs = runs)
}

# The processor, read from the operating system; NA when it cannot be read.
benchmark_cpu <- function() {
  sys <- Sys.info()[["sysname"]]
  out <- tryCatch({
    if (sys == "Windows") {
      x <- system2("powershell",
                   c("-NoProfile", "-Command",
                     shQuote("(Get-CimInstance Win32_Processor).Name")),
                   stdout = TRUE, stderr = FALSE)
    } else if (sys == "Darwin") {
      x <- system2("sysctl", c("-n", "machdep.cpu.brand_string"),
                   stdout = TRUE, stderr = FALSE)
    } else {
      x <- sub("^model name\\s*:\\s*", "",
               grep("^model name", readLines("/proc/cpuinfo"), value = TRUE))
    }
    x <- trimws(x)
    x[nzchar(x)][1]
  }, error = function(e) NA_character_, warning = function(w) NA_character_)
  if (length(out) == 0L || is.na(out)) {
    out <- Sys.getenv("PROCESSOR_IDENTIFIER", NA_character_)
  }
  out
}

# The date, the operating system, the processor and the version of R.
benchmark_environment <- function() {
  data.frame(
    date = format(Sys.time(), "%Y-%m-%d"),
    os = utils::sessionInfo()$running,
    machine = Sys.info()[["machine"]],
    cpu = benchmark_cpu(),
    cores = parallel::detectCores(logical = FALSE),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    stringsAsFactors = FALSE
  )
}

run_benchmark <- function(reps = 3, min_time = 1,
                          out_file = file.path("inst", "benchmark", "results",
                                               "benchmark.csv"),
                          boinetc_zip = NULL, takeda2023_zip = NULL) {
  if (!requireNamespace("simFastBOIN", quietly = TRUE)) {
    stop("install the package first", call. = FALSE)
  }
  if (requireNamespace("pkgload", quietly = TRUE) &&
      pkgload::is_dev_package("simFastBOIN")) {
    stop("simFastBOIN was loaded with devtools::load_all(); install it, ",
         "restart R and call library(simFastBOIN)", call. = FALSE)
  }
  if (!exists("sim_boin", mode = "function")) {
    stop("call library(simFastBOIN) first", call. = FALSE)
  }
  env <- benchmark_environment()
  version_of <- function(pkg) {
    if (is.na(pkg)) NA_character_ else as.character(utils::packageVersion(pkg))
  }
  cases <- benchmark_cases(
    boinetc_prog = benchmark_unzip(boinetc_zip, "BOINETC-main/prog"),
    takeda2023_prog = benchmark_unzip(takeda2023_zip, "TITE-gBOIN-ET/prog")
  )
  rows <- list()
  for (k in seq_along(cases)) {
    cs <- cases[[k]]
    message(sprintf("[%d/%d] %s, %s trials", k, length(cases), cs$design,
                    format(cs$n_trials, big.mark = ",")))
    has_ref <- !is.null(cs$reference)
    is_code <- !is.null(cs$available)
    ref_ok <- has_ref && (if (is_code) isTRUE(cs$available) else
      requireNamespace(cs$package, quietly = TRUE))
    missing_status <- if (is_code) "not available" else "not installed"
    for (r in seq_len(reps)) {
      impl <- list(list(implementation = "simFastBOIN", fun = cs$fun,
                        package = "simFastBOIN", run = cs$ours))
      if (has_ref || is_code) {
        impl[[2]] <- list(implementation = "reference",
                          fun = cs$reference_fun, package = cs$package,
                          run = if (ref_ok) cs$reference else NULL,
                          code = is_code)
      }
      for (im in impl) {
        tm <- if (is.null(im$run)) {
          c(seconds = NA_real_, runs = 0)
        } else {
          benchmark_time(im$run, cs$n_trials, min_time)
        }
        rows[[length(rows) + 1L]] <- cbind(data.frame(
          design = cs$design, setting = cs$setting, n_trials = cs$n_trials,
          implementation = im$implementation, fun = im$fun,
          package = im$package,
          version = if (is.null(im$run) || isTRUE(im$code)) {
            NA_character_
          } else {
            version_of(im$package)
          },
          rep = r, runs = unname(tm[["runs"]]),
          seconds = unname(tm[["seconds"]]),
          status = if (is.null(im$run)) missing_status else "ok",
          stringsAsFactors = FALSE
        ), env)
      }
    }
    res <- do.call(rbind, rows)
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(res, out_file, row.names = FALSE)
  }
  invisible(res)
}
