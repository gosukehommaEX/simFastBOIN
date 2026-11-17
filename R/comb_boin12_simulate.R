#' Simulate Comb-BOIN12 Trials
#'
#' @description
#'   Simulate trials of the Comb-BOIN12 design of Lu et al. (2025), the
#'   extension of BOIN12 to two drugs given together, which assigns each
#'   cohort to the most desirable combination by a utility that weighs the
#'   response against the DLT, and return the patients, DLTs and responses at
#'   every combination. This is the simulation engine behind
#'   \code{\link{sim_comb_boin12}}.
#'
#' @param target
#'   Numeric scalar. Highest acceptable DLT probability, phi_T. It is the
#'   target of the BOIN boundaries \code{lambda_e} and \code{lambda_d}, of
#'   the admissibility of a combination and of the MTD combination that
#'   bounds the OBD combination.
#'
#' @param target_resp
#'   Numeric scalar. Lowest acceptable response probability, phi_E, as in
#'   \code{\link{sim_boin12}}.
#'
#' @param p_true
#'   Numeric matrix. True DLT probability of every combination, with one row
#'   per level of drug A and one column per level of drug B.
#'
#' @param p_resp
#'   Numeric matrix of the same size. True response probability of every
#'   combination.
#'
#' @param n_cohort
#'   Integer scalar. Number of cohorts in a trial.
#'
#' @param cohort_size
#'   Integer scalar or vector. Number of patients per cohort. A scalar is used
#'   for every cohort; a shorter vector is padded with its last element.
#'
#' @param n_trials
#'   Integer scalar. Number of trials to simulate. Defaults to 10000.
#'
#' @param start_dose
#'   Integer vector of length 2. Levels of drug A and drug B of the first
#'   cohort. Defaults to \code{c(1, 1)}.
#'
#' @param n_earlystop
#'   Integer scalar or \code{NULL}. Number of patients at the current
#'   combination that ends the trial when the next cohort would stay there
#'   (see \code{n_earlystop_rule}). Defaults to \code{NULL}, which runs every
#'   trial to the maximum sample size, as in Lu et al. (2025).
#'
#' @param n_earlystop_rule
#'   Character scalar. \code{"with_stay"} (the default) ends the trial at
#'   \code{n_earlystop} only when the next cohort would stay at the current
#'   combination; \code{"simple"} ends it in any case.
#'
#' @param p_saf,p_tox
#'   Numeric scalars or \code{NULL}. Highest DLT probability deemed
#'   subtherapeutic and lowest deemed overly toxic, which give the BOIN
#'   boundaries. Default to \code{0.6 * target} and \code{1.4 * target}.
#'
#' @param cutoff_eli
#'   Numeric scalar. A combination is inadmissible, together with every
#'   combination with at least the same levels of both drugs, when the
#'   posterior probability that its DLT rate exceeds \code{target} is at
#'   least this cutoff. Defaults to 0.95.
#'
#' @param cutoff_eli_resp
#'   Numeric scalar. A combination is inadmissible when the posterior
#'   probability that its response rate is below \code{target_resp} is at
#'   least this cutoff. Defaults to 0.90.
#'
#' @param extrasafe
#'   Logical scalar. Apply the cutoff \code{cutoff_eli - offset} to the
#'   lowest combination, as Lu et al. (2025) do with 0.85. Defaults to
#'   \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar. Amount by which \code{cutoff_eli} is lowered at the
#'   lowest combination when \code{extrasafe} is \code{TRUE}. Defaults to
#'   0.05, as in \code{\link{sim_comb_boin}}; Lu et al. (2025) use 0.10.
#'
#' @param psi00,psi11
#'   Numeric scalars from 0 to 100. Utility of a patient without DLT and
#'   without response, and with both, as in \code{\link{sim_boin12}}. Default
#'   to 40 and 60. The utility \code{pE - wT * pT} of Lu et al. (2025) is
#'   \code{psi00 = 100 * wT / (1 + wT)} and \code{psi11 = 100 / (1 + wT)};
#'   their \code{wT = 0.11} is \code{psi00 = 9.90991} and
#'   \code{psi11 = 90.09009}.
#'
#' @param benchmark_weight
#'   Numeric scalar above 0 and at most 1. Weight of the lowest acceptable
#'   utility in the utility benchmark, as in \code{\link{sim_boin12}}; it is
#'   the weight \code{w_u} of Lu et al. (2025), who use 0.6. Defaults to
#'   0.5.
#'
#' @param utility_prior
#'   Numeric vector of two positive numbers. Beta prior of the utility of a
#'   combination divided by 100. Defaults to \code{c(1, 1)}.
#'
#' @param n_star
#'   Integer scalar. With at least this many patients at the current
#'   combination and a DLT rate between \code{lambda_e} and \code{lambda_d},
#'   the next combination is chosen from the current one and the two below
#'   it. Defaults to 6.
#'
#' @param tie_break
#'   Character scalar. Choice among the combinations of equal desirability:
#'   \code{"higher"} (the default) keeps those with the highest sum of the
#'   two levels and draws one of them at random when several remain;
#'   \code{"random"} draws one of all the tied combinations. See Details.
#'
#' @param no_candidate
#'   Character scalar. What happens when none of the combinations allowed
#'   for the next cohort is admissible: \code{"nearest"} (the default) goes
#'   to the admissible combination nearest to the current one;
#'   \code{"stop"} ends the trial without an OBD combination. See Details.
#'
#' @param resp_cor
#'   Numeric scalar between -1 and 1. Correlation between the DLT and the
#'   response of a patient, see \code{resp_cor_type}. Defaults to 0.
#'
#' @param resp_cor_type
#'   Character scalar. \code{"copula"} (the default) or \code{"binary"}, as
#'   in \code{\link{sim_boin12}}.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. The state of the calling
#'   session is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{comb_boin12_trials}, which is a list with
#'   components
#'   \item{n_pts, n_tox, n_resp, n_tox_resp}{Integer arrays with dimensions
#'     trial by level of drug A by level of drug B, of the patients treated,
#'     the DLTs, the responses and the patients with both a DLT and a
#'     response.}
#'   \item{cohorts_used}{Integer vector of the cohorts enrolled.}
#'   \item{stop_reason}{Character vector: \code{"no_admissible_dose"},
#'     \code{"n_earlystop"}, \code{"max_cohorts"} or
#'     \code{"no_admissible_candidate"} (only with
#'     \code{no_candidate = "stop"}).}
#'   \item{early_stop}{Logical vector, \code{TRUE} for the trials stopped
#'     because no combination was admissible or, with
#'     \code{no_candidate = "stop"}, because no combination allowed for the
#'     next cohort was admissible. These trials select no OBD combination.}
#'   \item{design}{\code{"Comb-BOIN12"}.}
#'   \item{boundary}{List with \code{lambda_e}, \code{lambda_d}, \code{u_low}
#'     and \code{u_benchmark}.}
#'   \item{settings}{List of the design parameters.}
#'
#' @details
#'   After each cohort, with \code{phat} the DLT rate at the current
#'   combination (j, k): the trial stops when no combination is admissible;
#'   when \code{phat} is at least \code{lambda_d}, the next combination is
#'   the most desirable of the admissible combinations below (j, k), that is
#'   with levels at most j and k, that have the highest sum of the two
#'   levels: (j - 1, k) or (j, k - 1) when either is admissible, as in Lu et
#'   al. (2025), and lower ones otherwise; the current one when there is
#'   none. When \code{phat} is above \code{lambda_e} and at least
#'   \code{n_star} patients have been treated, the next combination is the
#'   most desirable of (j, k), (j - 1, k) and (j, k - 1); otherwise of (j, k)
#'   and its four neighbors (j - 1, k), (j, k - 1), (j + 1, k) and
#'   (j, k + 1). Only admissible combinations inside the matrix are
#'   considered. When none of these is admissible, the admissible combination
#'   nearest to the current one (by the sum of the differences of the two
#'   levels) is chosen, then the one with the lowest sum of the levels, then
#'   the most desirable. The desirability is that of BOIN12
#'   (\code{\link{boin12_decision_table}}). Untried combinations have a DLT
#'   rate of 0 and are admissible unless excluded by the DLT data of a lower
#'   combination.
#'
#'   With \code{no_candidate = "stop"}, the trial stops instead, without an
#'   OBD combination and with the stopping reason
#'   \code{"no_admissible_candidate"}, when none of these combinations is
#'   admissible and when \code{phat} is at least \code{lambda_d} and no
#'   combination below the current one is admissible. At the lowest
#'   combination, a \code{phat} of at least \code{lambda_d} still keeps the
#'   next cohort there when it is admissible.
#'
#'   Lu et al. (2025) do not state how ties of desirability are broken nor
#'   what happens when no candidate is admissible. With
#'   \code{tie_break = "higher"}, among tied combinations those with the
#'   highest sum of the levels are kept, so that the design escalates when an
#'   untried combination is as desirable as the current one, and a tie that
#'   remains (between (j + 1, k) and (j, k + 1), for example) is broken at
#'   random. The allocation percentages published by Lu et al. (2025) agree
#'   with this rule, and not with a random choice among all the tied
#'   combinations nor with a fixed preference for one drug (see
#'   \code{inst/reproduce/reproduce_Lu2025.R}). The two settings of
#'   \code{no_candidate} could not be told apart from the published results.
#'   With one row or one column, the trials are those of
#'   \code{\link{boin12_simulate}} without its exploration rule; Lu et al.
#'   (2025) have no such rule.
#'
#'   When \code{psi00 + psi11} is 100, as for the utility of Lu et al.
#'   (2025), the desirability does not depend on the number of patients with
#'   both a DLT and a response, and combinations with the same numbers of
#'   patients, DLTs and responses tie exactly, whatever the rounding of
#'   \code{psi00} and \code{psi11}.
#'
#'   Each patient draws two uniform variates in enrollment order, the first
#'   for the DLT and the second for the response. A decision between tied
#'   combinations draws one more variate when a tie remains after the rule
#'   of \code{tie_break}.
#'
#' @references
#'   Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
#'   utility-based Bayesian optimal interval design for dose optimization in
#'   cancer drug-combination trials. Statistics in Biopharmaceutical
#'   Research, 17(2), 266-276. \doi{10.1080/19466315.2024.2370403}
#'
#'   Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
#'   optimal interval phase I/II trial design for utility-based dose finding
#'   in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
#'   1393-1402. \doi{10.1200/PO.20.00257}
#'
#' @examples
#' # Scenario 16 of Lu et al. (2025), with their utility (wT = 0.11)
#' p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
#' p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
#' trials <- comb_boin12_simulate(
#'   target = 0.33, target_resp = 0.2, p_true = p_true, p_resp = p_resp,
#'   n_cohort = 14, cohort_size = 3, n_trials = 200,
#'   psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11, benchmark_weight = 0.6,
#'   extrasafe = TRUE, offset = 0.10, seed = 123
#' )
#' trials
#'
#' @seealso \code{\link{sim_comb_boin12}}, \code{\link{comb_boin12_select_obd}},
#'   \code{\link{boin12_simulate}}
#'
#' @export
comb_boin12_simulate <- function(target, target_resp, p_true, p_resp,
                                 n_cohort, cohort_size, n_trials = 10000,
                                 start_dose = c(1, 1), n_earlystop = NULL,
                                 n_earlystop_rule = c("with_stay", "simple"),
                                 p_saf = NULL, p_tox = NULL,
                                 cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
                                 extrasafe = FALSE, offset = 0.05,
                                 psi00 = 40, psi11 = 60,
                                 benchmark_weight = 0.5,
                                 utility_prior = c(1, 1), n_star = 6,
                                 tie_break = c("higher", "random"),
                                 no_candidate = c("nearest", "stop"),
                                 resp_cor = 0,
                                 resp_cor_type = c("copula", "binary"),
                                 seed = 123) {

  n_earlystop_rule <- match.arg(n_earlystop_rule)
  tie_break <- match.arg(tie_break)
  no_candidate <- match.arg(no_candidate)
  resp_cor_type <- match.arg(resp_cor_type)

  check_p_true_matrix(p_true)
  check_p_resp_matrix(p_resp, p_true)
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_start_combination(start_dose, p_true)
  check_count(n_star, "n_star", 1L)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  check_flag(extrasafe, "extrasafe")
  check_comb_offset(offset, cutoff_eli)
  check_scalar_prob(resp_cor, "resp_cor", lower = -1, upper = 1)
  des <- boin12_check_design(target, target_resp, p_saf = p_saf,
                             p_tox = p_tox, psi00 = psi00, psi11 = psi11,
                             benchmark_weight = benchmark_weight,
                             utility_prior = utility_prior)

  n_rows <- nrow(p_true)
  n_cols <- ncol(p_true)
  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)
  if (!is.null(n_earlystop)) check_count(n_earlystop, "n_earlystop", 1L)

  adm <- comb_boin12_tables(target, target_resp, cutoff_eli, cutoff_eli_resp,
                            extrasafe, offset, max_total_pts)
  if (resp_cor_type == "binary") {
    link <- boin_et_resp_link(as.vector(p_true), as.vector(p_resp), resp_cor)
  } else {
    link <- list(q0 = rep(0, n_rows * n_cols), q1 = rep(0, n_rows * n_cols))
  }

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- comb_boin12_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    p_resp = as.numeric(p_resp),
    q0 = as.numeric(link$q0),
    q1 = as.numeric(link$q1),
    cor_type = if (resp_cor_type == "binary") 1L else 0L,
    resp_cor = as.numeric(resp_cor),
    n_rows = as.integer(n_rows),
    n_cols = as.integer(n_cols),
    cohort_size = cohort_size_vec,
    start_row = as.integer(start_dose[1L]),
    start_col = as.integer(start_dose[2L]),
    n_earlystop = if (is.null(n_earlystop)) 0L else as.integer(n_earlystop),
    early_stop_simple = identical(n_earlystop_rule, "simple"),
    n_star = as.integer(n_star),
    tie_break = if (tie_break == "random") 1L else 0L,
    no_candidate = if (no_candidate == "stop") 1L else 0L,
    lambda_e = des$lambda_e,
    lambda_d = des$lambda_d,
    psi00 = as.numeric(psi00),
    psi11 = as.numeric(psi11),
    alpha = des$utility_prior[1L],
    beta = des$utility_prior[2L],
    u_benchmark = des$u_benchmark / 100,
    max_n = as.integer(max_total_pts),
    tox_ok = adm$tox_ok,
    tox_ok_lowest = adm$tox_ok_lowest,
    resp_ok = adm$resp_ok
  )

  trial_names <- list(NULL, paste0("A", seq_len(n_rows)),
                      paste0("B", seq_len(n_cols)))
  for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
    dimnames(res[[nm]]) <- trial_names
  }

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$n_tox,
      n_resp = res$n_resp,
      n_tox_resp = res$n_tox_resp,
      cohorts_used = res$cohorts_used,
      stop_reason = boin12_stop_reasons()[res$stop_code + 1L],
      early_stop = res$stop_code %in% c(0L, 3L),
      design = "Comb-BOIN12",
      boundary = list(lambda_e = des$lambda_e, lambda_d = des$lambda_d,
                      u_low = des$u_low, u_benchmark = des$u_benchmark),
      settings = list(
        target = target,
        target_resp = target_resp,
        p_true = p_true,
        p_resp = p_resp,
        n_cohort = n_cohort,
        cohort_size = cohort_size_vec,
        n_trials = n_trials,
        start_dose = as.integer(start_dose),
        n_earlystop = n_earlystop,
        n_earlystop_rule = n_earlystop_rule,
        p_saf = des$p_saf,
        p_tox = des$p_tox,
        cutoff_eli = cutoff_eli,
        cutoff_eli_resp = cutoff_eli_resp,
        extrasafe = extrasafe,
        offset = offset,
        psi00 = psi00,
        psi11 = psi11,
        benchmark_weight = benchmark_weight,
        utility_prior = des$utility_prior,
        n_star = n_star,
        tie_break = tie_break,
        no_candidate = no_candidate,
        resp_cor = resp_cor,
        resp_cor_type = resp_cor_type,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "comb_boin12_trials"
  )
}
