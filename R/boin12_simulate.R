#' Simulate BOIN12 Trials
#'
#' @description
#'   Simulate trials of the BOIN12 design of Lin et al. (2020), which assigns
#'   each cohort to the most desirable dose by a utility that weighs the
#'   response against the DLT, and return the patients, DLTs and responses at
#'   every dose. This is the simulation engine behind \code{\link{sim_boin12}}.
#'
#' @param target
#'   Numeric scalar. Highest acceptable DLT probability, called phi_T by Lin
#'   et al. (2020). It is the target of the BOIN boundaries
#'   \code{lambda_e} and \code{lambda_d}, of the admissibility of a dose and
#'   of the MTD that bounds the OBD. Lin et al. recommend a value slightly
#'   above the target of a design that only considers toxicity, for example
#'   0.35 instead of 0.30.
#'
#' @param target_resp
#'   Numeric scalar. Lowest acceptable response probability, phi_E. Note that
#'   this is a lower limit, as in \code{\link{sim_bf_boin_et}}, whereas the
#'   \code{target_resp} of \code{\link{sim_boin_et}} is a target.
#'
#' @param p_true
#'   Numeric vector. True DLT probability at each dose.
#'
#' @param p_resp
#'   Numeric vector. True response probability at each dose.
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
#'   Integer scalar. Dose level of the first cohort. Defaults to 1.
#'
#' @param n_earlystop
#'   Integer scalar or \code{NULL}. Number of patients at the current dose
#'   that ends the trial when the next cohort would stay there (see
#'   \code{n_earlystop_rule}). Defaults to \code{NULL}, which runs every
#'   trial to the maximum sample size, as in Lin et al. (2020).
#'
#' @param n_earlystop_rule
#'   Character scalar. \code{"with_stay"} (the default) ends the trial at
#'   \code{n_earlystop} only when the next cohort would stay at the current
#'   dose, as in the example of Zhou et al. (2022); \code{"simple"} ends it in
#'   any case.
#'
#' @param p_saf,p_tox
#'   Numeric scalars or \code{NULL}. Highest DLT probability deemed
#'   subtherapeutic and lowest deemed overly toxic, which give the BOIN
#'   boundaries. Default to \code{0.6 * target} and \code{1.4 * target}.
#'
#' @param cutoff_eli
#'   Numeric scalar. A dose is inadmissible, together with every higher dose,
#'   when the posterior probability that its DLT rate exceeds \code{target}
#'   is at least this cutoff. Defaults to 0.95.
#'
#' @param cutoff_eli_resp
#'   Numeric scalar. A dose is inadmissible when the posterior probability
#'   that its response rate is below \code{target_resp} is at least this
#'   cutoff. Defaults to 0.90.
#'
#' @param psi00,psi11
#'   Numeric scalars from 0 to 100. Utility of a patient without DLT and
#'   without response, and with both. The utility of a response without DLT
#'   is 100 and that of a DLT without response is 0. Default to 40 and 60, the
#'   example of Lin et al. (2020).
#'
#' @param benchmark_weight
#'   Numeric scalar above 0 and at most 1. The utility benchmark is
#'   \code{benchmark_weight * u_low + (1 - benchmark_weight) * 100}, where
#'   \code{u_low} is the expected utility of a dose whose DLT and response
#'   probabilities are \code{target} and \code{target_resp}, the two
#'   outcomes being independent. Defaults to 0.5, the benchmark of Lin et al.
#'   (2020); Lu et al. (2025) use 0.6.
#'
#' @param utility_prior
#'   Numeric vector of two positive numbers. Beta prior of the utility of a
#'   dose divided by 100. Defaults to \code{c(1, 1)}.
#'
#' @param n_star
#'   Integer scalar. With at least this many patients at the current dose and
#'   a DLT rate between \code{lambda_e} and \code{lambda_d}, the next dose is
#'   chosen from the current dose and the dose below. Defaults to 6.
#'
#' @param n_explore
#'   Integer scalar. With at least this many patients at the current dose, a
#'   DLT rate below \code{lambda_d} and an untried admissible dose above, the
#'   next cohort goes to that dose. Defaults to 9.
#'
#' @param tie_break
#'   Character scalar. Choice among the doses of equal desirability:
#'   \code{"higher"} (the default) takes the highest of them, the rule of the
#'   BOIN12 software according to the 'escalation' package; \code{"random"}
#'   takes one of them at random, as stated in Table S2 of the Supplement of
#'   Lin et al. (2020). See Details.
#'
#' @param no_candidate
#'   Character scalar. What happens when none of the doses allowed for the
#'   next cohort is admissible: \code{"nearest"} (the default) goes to the
#'   admissible dose nearest to the current one; \code{"stop"} ends the trial
#'   without an OBD. See Details.
#'
#' @param explore_from
#'   Character scalar. Doses from which the exploration rule (see
#'   \code{n_explore}) escalates: \code{"any"} (the default) from any dose, as
#'   in the 'escalation' package; \code{"admissible"} only from an admissible
#'   dose, as the results published by Lin et al. (2020) indicate. See
#'   Details.
#'
#' @param resp_cor
#'   Numeric scalar between -1 and 1. Correlation between the DLT and the
#'   response of a patient, see \code{resp_cor_type}. Defaults to 0.
#'
#' @param resp_cor_type
#'   Character scalar. \code{"copula"} (the default) links the two outcomes
#'   by a Gaussian copula with correlation \code{resp_cor}, the bivariate
#'   normal model of Lin et al. (2020); \code{"binary"} sets the correlation
#'   of the two binary outcomes, as in \code{\link{boin_et_simulate}}.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. The state of the calling
#'   session is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{boin12_trials}, which is a list with components
#'   \item{n_pts, n_tox, n_resp, n_tox_resp}{Integer matrices, one row per
#'     trial and one column per dose, of the patients treated, the DLTs, the
#'     responses and the patients with both a DLT and a response.}
#'   \item{cohorts_used}{Integer vector of the cohorts enrolled.}
#'   \item{stop_reason}{Character vector: \code{"no_admissible_dose"},
#'     \code{"n_earlystop"}, \code{"max_cohorts"} or
#'     \code{"no_admissible_candidate"} (only with
#'     \code{no_candidate = "stop"}).}
#'   \item{early_stop}{Logical vector, \code{TRUE} for the trials stopped
#'     because no dose was admissible or, with \code{no_candidate = "stop"},
#'     because no dose allowed for the next cohort was admissible. These
#'     trials select no OBD.}
#'   \item{design}{\code{"BOIN12"}.}
#'   \item{boundary}{List with \code{lambda_e}, \code{lambda_d}, \code{u_low}
#'     and \code{u_benchmark}.}
#'   \item{settings}{List of the design parameters.}
#'
#' @details
#'   After each cohort, with \code{phat} the DLT rate at the current dose:
#'   the trial stops when no dose is admissible; when \code{n_explore}
#'   patients have been treated at the current dose, \code{phat} is below
#'   \code{lambda_d} and the dose above is untried, the dose is escalated
#'   (with \code{explore_from = "admissible"}, only when the current dose is
#'   admissible);
#'   otherwise, when \code{phat} is at least \code{lambda_d}, the next dose is
#'   the highest admissible dose below the current one (the current one if
#'   there is none); when \code{phat} is above \code{lambda_e} and at least
#'   \code{n_star} patients have been treated, the next dose is the more
#'   desirable of the current dose and the dose below; otherwise it is the
#'   most desirable of the dose below, the current dose and the dose above.
#'   Only admissible doses are considered. Ties go to the higher dose, as in
#'   the BOIN12 software according to the 'escalation' package, or with
#'   \code{tie_break = "random"} to one of the tied doses drawn at random.
#'   When none of these doses is admissible, the admissible dose nearest to
#'   the current one is chosen, ties to the lower dose. The desirability is
#'   that of \code{\link{boin12_decision_table}}.
#'
#'   With \code{no_candidate = "stop"}, the trial stops instead, without an
#'   OBD and with the stopping reason \code{"no_admissible_candidate"}, when
#'   none of these doses is admissible and when \code{phat} is at least
#'   \code{lambda_d} and no dose below the current one is admissible. At the
#'   lowest dose, a \code{phat} of at least \code{lambda_d} still keeps the
#'   next cohort there when that dose is admissible.
#'
#'   Lin et al. (2020) leave three points open, decided here as follows. A
#'   dose made inadmissible by its DLT data also makes every higher dose
#'   inadmissible, as in BOIN and in Lu et al. (2025), whereas the
#'   'escalation' package only checks the doses next to the current one. A
#'   de-escalation never goes to an inadmissible dose; 'escalation' moves to
#'   the dose below whatever its data. The trial stops when no dose is
#'   admissible; 'escalation' goes on at the lowest dose when its DLT rate is
#'   high. The article breaks ties at random; the software chooses the higher
#'   dose. \code{tie_break = "random"}, \code{no_candidate = "stop"} and
#'   \code{explore_from = "admissible"} are the settings compared with the
#'   results published by Lin et al. (2020) in
#'   \code{inst/reproduce/reproduce_Lin2020.R}; the last two rules are not
#'   stated in the article and were inferred from those results.
#'
#'   Each patient draws two uniform variates in enrollment order, the first
#'   for the DLT and the second for the response. With
#'   \code{tie_break = "random"}, a decision among tied doses draws one more
#'   variate, so that the trials depart from those with \code{"higher"} after
#'   the first tie.
#'
#' @references
#'   Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
#'   optimal interval phase I/II trial design for utility-based dose finding
#'   in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
#'   1393-1402. \doi{10.1200/PO.20.00257}
#'
#'   Zhou, Y., Lin, R., Lee, J. J., Li, D., Wang, L., Li, R. and Yuan, Y.
#'   (2022). TITE-BOIN12: A Bayesian phase I/II trial design to find the
#'   optimal biological dose with late-onset toxicity and efficacy. Statistics
#'   in Medicine, 41(11), 1918-1931. \doi{10.1002/sim.9337}
#'
#'   Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
#'   utility-based Bayesian optimal interval design for dose optimization in
#'   cancer drug-combination trials. Statistics in Biopharmaceutical
#'   Research, 17(2), 266-276. \doi{10.1080/19466315.2024.2370403}
#'
#' @examples
#' # Scenario 1 of Lin et al. (2020)
#' trials <- boin12_simulate(
#'   target = 0.35, target_resp = 0.25,
#'   p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'   p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
#'   n_cohort = 12, cohort_size = 3, n_trials = 200, seed = 123
#' )
#' trials
#'
#' @seealso \code{\link{sim_boin12}}, \code{\link{boin12_select_obd}},
#'   \code{\link{boin12_decision_table}}
#'
#' @export
boin12_simulate <- function(target, target_resp, p_true, p_resp, n_cohort,
                            cohort_size, n_trials = 10000, start_dose = 1,
                            n_earlystop = NULL,
                            n_earlystop_rule = c("with_stay", "simple"),
                            p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                            cutoff_eli_resp = 0.90, psi00 = 40, psi11 = 60,
                            benchmark_weight = 0.5, utility_prior = c(1, 1),
                            n_star = 6, n_explore = 9,
                            tie_break = c("higher", "random"),
                            no_candidate = c("nearest", "stop"),
                            explore_from = c("any", "admissible"),
                            resp_cor = 0,
                            resp_cor_type = c("copula", "binary"),
                            seed = 123) {

  n_earlystop_rule <- match.arg(n_earlystop_rule)
  tie_break <- match.arg(tie_break)
  no_candidate <- match.arg(no_candidate)
  explore_from <- match.arg(explore_from)
  resp_cor_type <- match.arg(resp_cor_type)

  check_p_true(p_true)
  n_doses <- length(p_true)
  if (!is.numeric(p_resp) || length(p_resp) != n_doses ||
      any(!is.finite(p_resp)) || any(p_resp < 0) || any(p_resp > 1)) {
    stop("'p_resp' must hold one response probability between 0 and 1 ",
         "per dose in 'p_true'", call. = FALSE)
  }
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(start_dose, "start_dose", 1L)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of doses in 'p_true'",
         call. = FALSE)
  }
  check_count(n_star, "n_star", 1L)
  check_count(n_explore, "n_explore", 1L)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  check_scalar_prob(resp_cor, "resp_cor", lower = -1, upper = 1)
  des <- boin12_check_design(target, target_resp, p_saf = p_saf,
                             p_tox = p_tox, psi00 = psi00, psi11 = psi11,
                             benchmark_weight = benchmark_weight,
                             utility_prior = utility_prior)

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)
  if (!is.null(n_earlystop)) check_count(n_earlystop, "n_earlystop", 1L)

  adm <- boin12_admissible_table(target, target_resp, cutoff_eli,
                                 cutoff_eli_resp, max_total_pts)
  if (resp_cor_type == "binary") {
    link <- boin_et_resp_link(p_true, p_resp, resp_cor)
  } else {
    link <- list(q0 = rep(0, n_doses), q1 = rep(0, n_doses))
  }

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- boin12_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    p_resp = as.numeric(p_resp),
    q0 = as.numeric(link$q0),
    q1 = as.numeric(link$q1),
    cor_type = if (resp_cor_type == "binary") 1L else 0L,
    resp_cor = as.numeric(resp_cor),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    n_earlystop = if (is.null(n_earlystop)) 0L else as.integer(n_earlystop),
    early_stop_simple = identical(n_earlystop_rule, "simple"),
    n_star = as.integer(n_star),
    n_explore = as.integer(n_explore),
    tie_break = if (tie_break == "random") 1L else 0L,
    no_candidate = if (no_candidate == "stop") 1L else 0L,
    explore_from = if (explore_from == "admissible") 1L else 0L,
    lambda_e = des$lambda_e,
    lambda_d = des$lambda_d,
    psi00 = as.numeric(psi00),
    psi11 = as.numeric(psi11),
    alpha = des$utility_prior[1L],
    beta = des$utility_prior[2L],
    u_benchmark = des$u_benchmark / 100,
    max_n = as.integer(max_total_pts),
    tox_ok = adm$tox_ok,
    resp_ok = adm$resp_ok
  )

  dose_names <- paste0("DL", seq_len(n_doses))
  for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp")) {
    colnames(res[[nm]]) <- dose_names
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
      design = "BOIN12",
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
        start_dose = start_dose,
        n_earlystop = n_earlystop,
        n_earlystop_rule = n_earlystop_rule,
        p_saf = des$p_saf,
        p_tox = des$p_tox,
        cutoff_eli = cutoff_eli,
        cutoff_eli_resp = cutoff_eli_resp,
        psi00 = psi00,
        psi11 = psi11,
        benchmark_weight = benchmark_weight,
        utility_prior = des$utility_prior,
        n_star = n_star,
        n_explore = n_explore,
        tie_break = tie_break,
        no_candidate = no_candidate,
        explore_from = explore_from,
        resp_cor = resp_cor,
        resp_cor_type = resp_cor_type,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "boin12_trials"
  )
}
