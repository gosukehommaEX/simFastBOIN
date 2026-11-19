#' Simulate BOIN-ETC Trials
#'
#' @description
#'   Simulate trials of the BOIN-ETC design of Kakizume et al. (2024), which
#'   finds the optimal dose combination of two drugs given together by
#'   dividing the dose matrix into subtrials along which the DLT probability
#'   increases (the waterfall approach) and running the BOIN-ET design of
#'   Takeda et al. (2018) within each subtrial. Return the patients, DLTs
#'   and responses at every combination. This is the simulation engine behind
#'   \code{\link{sim_boin_etc}}.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, phi.
#'
#' @param target_resp
#'   Numeric scalar. Target response probability, delta, as in
#'   \code{\link{sim_boin_et}}.
#'
#' @param p_true
#'   Numeric matrix. True DLT probability of every combination, with one row
#'   per level of drug A and one column per level of drug B. At least two
#'   levels of each drug.
#'
#' @param p_resp
#'   Numeric matrix of the same size. True response probability of every
#'   combination.
#'
#' @param n_cohort
#'   Integer scalar. Number of cohorts in a trial; with \code{cohort_size}
#'   it gives the maximum sample size.
#'
#' @param cohort_size
#'   Integer scalar or vector. Number of patients per cohort. A scalar is used
#'   for every cohort; a shorter vector is padded with its last element.
#'
#' @param n_trials
#'   Integer scalar. Number of trials to simulate. Defaults to 10000.
#'
#' @param subtrials
#'   Character scalar. The division of the dose matrix into subtrials and
#'   their starting combinations: \code{"etc1"} (the default),
#'   \code{"etc2"} or \code{"etc3"}, BOIN-ETC1, 2 and 3 of the figures and
#'   tables of Kakizume et al. (2024). See Details.
#'
#' @param n_earlystop
#'   Integer scalar. A subtrial ends when the combination chosen for the next
#'   cohort has this many patients (n_stop of Kakizume et al., 2024).
#'   Defaults to 6.
#'
#' @param etc1_start
#'   Character scalar. Starting combination of a row subtrial of
#'   \code{"etc1"}: \code{"diagonal"} (the default) one level of drug B below
#'   the candidate of the previous subtrial, as the text and the figures of
#'   Kakizume et al. (2024) state; \code{"below"} the same level of drug B,
#'   the rule of the code of the authors that produced their published
#'   results. See Details.
#'
#' @param rule4
#'   Character scalar. Choice among the neighboring combinations when the DLT
#'   rate is between \code{lambda1} and \code{lambda2} and the response rate
#'   at most \code{eta1}: \code{"paper"} (the default) the rule of the text;
#'   \code{"authors_code"} the rule of the code of the authors, needed to
#'   reproduce their results for BOIN-ETC1. See Details.
#'
#' @param p_saf,p_tox,p_saf_resp,prior
#'   Design parameters of the BOIN-ET boundaries, as in
#'   \code{\link{boin_et_boundary}}. The defaults give
#'   \code{(lambda1, lambda2, eta1) = (0.14, 0.35, 0.48)} for a target DLT
#'   probability of 0.30 and a target response probability of 0.60.
#'
#' @param cutoff_eli
#'   Numeric scalar. A combination is excluded, with the higher combinations
#'   of its subtrial, when the posterior probability that its DLT rate
#'   exceeds \code{p_eli} is above this cutoff. Defaults to 0.95.
#'
#' @param p_eli
#'   Numeric scalar or \code{NULL}. DLT probability of the exclusion rule.
#'   Defaults to \code{NULL}, \code{target}.
#'
#' @param resp_cor
#'   Numeric scalar between -1 and 1. Correlation between the DLT and the
#'   response of a patient, see \code{resp_cor_type}. Defaults to 0, as in
#'   the simulations of Kakizume et al. (2024); the BOIN-ET functions default
#'   to 0.2.
#'
#' @param resp_cor_type
#'   Character scalar. \code{"copula"} (the default) or \code{"binary"}, as
#'   in \code{\link{sim_boin_et}}.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. The state of the calling
#'   session is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{boin_etc_trials}, which is a list with
#'   components
#'   \item{n_pts, n_tox, n_resp, n_tox_resp}{Integer arrays with dimensions
#'     trial by level of drug A by level of drug B, of the patients treated,
#'     the DLTs, the responses and the patients with both a DLT and a
#'     response.}
#'   \item{excluded}{Logical array of the same dimensions, \code{TRUE} for the
#'     combinations excluded by the DLT rule during the trial.}
#'   \item{cohorts_used, n_subtrials}{Integer vectors of the cohorts enrolled
#'     and of the subtrials run.}
#'   \item{stop_reason}{Character vector: \code{"lowest_excluded"} (the
#'     lowest combination was excluded), \code{"subtrials_completed"} or
#'     \code{"max_cohorts"} (the maximum sample size was reached).}
#'   \item{early_stop}{Logical vector, \code{TRUE} for the trials stopped
#'     because the lowest combination was excluded. These trials select no
#'     optimal dose combination.}
#'   \item{design}{\code{"BOIN-ETC1"}, \code{"BOIN-ETC2"} or
#'     \code{"BOIN-ETC3"}.}
#'   \item{boundary}{List with \code{lambda1}, \code{lambda2} and
#'     \code{eta1}.}
#'   \item{settings}{List of the design parameters.}
#'
#' @details
#'   With \code{L} levels of drug A and \code{M} of drug B (\code{AaBb} is the
#'   combination of level \code{a} of drug A and level \code{b} of drug B),
#'   the first subtrial is \code{A1B1, A2B1, ..., ALB1, ALB2, ..., ALBM},
#'   started at \code{A1B1}. Then:
#'   \describe{
#'     \item{\code{"etc1"}, \code{"etc2"}}{one subtrial per row
#'       \code{l = L - 1, ..., 1}, \code{AlB2, ..., AlBM}. When the candidate
#'       of the first subtrial is \code{AlB1}, the rows start at row
#'       \code{l - 1} (at row 1 when \code{l} is 1) and the rows above are
#'       skipped. \code{"etc2"} starts every row at its lowest combination.
#'       \code{"etc1"} starts the row below the candidate \code{AlBm} of the
#'       previous subtrial at level \code{m - 1} of drug B
#'       (\code{etc1_start = "diagonal"}) or at level \code{m}
#'       (\code{"below"}), and at level 2 when that level is 1, after a
#'       candidate in the first column and after a subtrial that ended
#'       without one.}
#'     \item{\code{"etc3"}}{L-shaped subtrials \code{j = 1, ...,
#'       min(L, M) - 1}: column \code{j + 1} of drug B from level 1 of drug A
#'       up to level \code{L - j}, then row \code{L - j} from column
#'       \code{j + 2} to \code{M}, each started at its lowest combination.
#'       With \code{L > M} the shapes run out of columns after \code{M - 1}
#'       of them, which still covers the matrix.}
#'   }
#'   Within a subtrial, after each cohort at a combination with DLT rate
#'   \code{p} and response rate \code{q}: the combination and the higher
#'   combinations of the subtrial are excluded when its DLT data fail the
#'   rule of \code{cutoff_eli}; the subtrial ends without a candidate when
#'   its lowest combination is excluded, and the trial ends when that
#'   combination is \code{A1B1}. Otherwise the next cohort escalates when
#'   \code{p <= lambda1} and \code{q <= eta1} (it stays at the highest
#'   combination of the subtrial and below an excluded one), de-escalates
#'   when \code{p >= lambda2} (it stays at the lowest), stays when
#'   \code{q > eta1}, and otherwise (rule 4) goes to the next combination up
#'   if it has not been treated, or else to the treated, not excluded
#'   combination of the current one and its two neighbors in the subtrial
#'   with the highest response rate, a tie being broken at random. The
#'   subtrial ends when the combination chosen for the next cohort has
#'   \code{n_earlystop} patients (that combination is the candidate of the
#'   subtrial), and the trial ends when all cohorts have been enrolled or
#'   all subtrials have been run.
#'
#'   Kakizume et al. (2024) do not state the rules in every case; the
#'   following were checked against their code and their published results
#'   (\code{inst/reproduce/reproduce_Kakizume2024.R}):
#'   \itemize{
#'     \item The text of Section 2.2 describes BOIN-ETC2 with the L-shaped
#'       subtrials and BOIN-ETC3 with the rows started at their lowest
#'       combination. Their Figure 2 and Supplementary Figures S1 and S2,
#'       their tables and their code have it the other way round, and so has
#'       this function.
#'     \item The text and the figures start a row of BOIN-ETC1 one level of
#'       drug B below the previous candidate (\code{"diagonal"}); the code
#'       that produced the tables starts it at the same level
#'       (\code{"below"}), and only that rule reproduces the published
#'       allocation.
#'     \item In rule 4, the code of the authors drops the response rate of an
#'       untreated combination below the current one but not the
#'       combination, so that the choice moves one combination down (from
#'       the current combination to the untreated one, or from the one above
#'       to the current one); it also breaks a tie among the first eligible
#'       combinations rather than among the tied ones.
#'       \code{rule4 = "authors_code"} applies these rules. An untreated
#'       combination below the current one only happens in BOIN-ETC1, and
#'       its published allocation is reproduced only with these rules.
#'     \item The code excludes only the current combination and may treat an
#'       excluded lowest combination again; here the text is followed (the
#'       published results do not depend on the difference). The text does
#'       not say what happens when \code{A1B1} is excluded.
#'   }
#'
#'   The response of a patient does not exclude a combination, as in
#'   Kakizume et al. (2024). Each patient draws two uniform variates in
#'   enrollment order, the first for the DLT and the second for the response;
#'   a tie under rule 4 draws one more after the cohort.
#'
#' @references
#'   Kakizume, T., Takeda, K., Taguri, M. and Morita, S. (2024). BOIN-ETC: A
#'   Bayesian optimal interval design considering efficacy and toxicity to
#'   identify the optimal dose combinations. Statistical Methods in Medical
#'   Research, 33(4), 716-727. \doi{10.1177/09622802241236936}
#'
#'   Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
#'   interval design for dose finding based on both efficacy and toxicity
#'   outcomes. Pharmaceutical Statistics, 17(4), 383-395.
#'   \doi{10.1002/pst.1864}
#'
#' @examples
#' # Scenario 1 of Kakizume et al. (2024)
#' p_true <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
#' p_resp <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
#' trials <- boin_etc_simulate(
#'   target = 0.30, target_resp = 0.60, p_true = p_true, p_resp = p_resp,
#'   n_cohort = 8, cohort_size = 3, n_trials = 200, seed = 123
#' )
#' trials
#'
#' @seealso \code{\link{sim_boin_etc}}, \code{\link{boin_etc_select_obd}},
#'   \code{\link{boin_et_simulate}}
#'
#' @export
boin_etc_simulate <- function(target, target_resp, p_true, p_resp, n_cohort,
                              cohort_size, n_trials = 10000,
                              subtrials = c("etc1", "etc2", "etc3"),
                              n_earlystop = 6,
                              etc1_start = c("diagonal", "below"),
                              rule4 = c("paper", "authors_code"),
                              p_saf = NULL, p_tox = NULL, p_saf_resp = NULL,
                              prior = rep(1 / 6, 6), cutoff_eli = 0.95,
                              p_eli = NULL, resp_cor = 0,
                              resp_cor_type = c("copula", "binary"),
                              seed = 123) {

  subtrials <- match.arg(subtrials)
  etc1_start <- match.arg(etc1_start)
  rule4 <- match.arg(rule4)
  resp_cor_type <- match.arg(resp_cor_type)

  check_p_true_matrix(p_true)
  check_p_resp_matrix(p_resp, p_true)
  if (nrow(p_true) < 2L || ncol(p_true) < 2L) {
    stop("'p_true' must have at least two levels of each drug", call. = FALSE)
  }
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(n_earlystop, "n_earlystop", 1L)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(resp_cor, "resp_cor", lower = -1, upper = 1)
  bd <- boin_et_boundary(target, target_resp, p_saf = p_saf, p_tox = p_tox,
                         p_saf_resp = p_saf_resp, prior = prior)
  if (is.null(p_eli)) p_eli <- target
  check_scalar_prob(p_eli, "p_eli")

  n_rows <- nrow(p_true)
  n_cols <- ncol(p_true)
  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)

  dec <- as.integer(boin_et_decision_codes(bd$lambda1, bd$lambda2, bd$eta1,
                                           max_total_pts))
  tox_ok <- boin_etc_tox_table(p_eli, cutoff_eli, max_total_pts)
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

  res <- boin_etc_simulate_cpp(
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
    subtrials = match(subtrials, c("etc1", "etc2", "etc3")) - 1L,
    start_rule = if (etc1_start == "below") 1L else 0L,
    rule4 = if (rule4 == "authors_code") 1L else 0L,
    n_stop = as.integer(n_earlystop),
    max_n = as.integer(max_total_pts),
    dec = dec,
    tox_ok = tox_ok
  )

  trial_names <- list(NULL, paste0("A", seq_len(n_rows)),
                      paste0("B", seq_len(n_cols)))
  for (nm in c("n_pts", "n_tox", "n_resp", "n_tox_resp", "excluded")) {
    dimnames(res[[nm]]) <- trial_names
  }
  excluded <- res$excluded != 0L
  dim(excluded) <- dim(res$excluded)
  dimnames(excluded) <- trial_names
  stop_reason <- boin_etc_stop_reasons()[res$stop_code + 1L]

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$n_tox,
      n_resp = res$n_resp,
      n_tox_resp = res$n_tox_resp,
      excluded = excluded,
      cohorts_used = res$cohorts_used,
      n_subtrials = res$n_subtrials,
      stop_reason = stop_reason,
      early_stop = stop_reason == "lowest_excluded",
      design = c(etc1 = "BOIN-ETC1", etc2 = "BOIN-ETC2",
                 etc3 = "BOIN-ETC3")[[subtrials]],
      boundary = list(lambda1 = bd$lambda1, lambda2 = bd$lambda2,
                      eta1 = bd$eta1),
      settings = list(
        target = target,
        target_resp = target_resp,
        p_true = p_true,
        p_resp = p_resp,
        n_cohort = n_cohort,
        cohort_size = cohort_size_vec,
        n_trials = n_trials,
        subtrials = subtrials,
        n_earlystop = n_earlystop,
        etc1_start = etc1_start,
        rule4 = rule4,
        p_saf = bd$p_saf,
        p_tox = bd$p_tox,
        p_saf_resp = bd$p_saf_resp,
        prior = prior,
        cutoff_eli = cutoff_eli,
        p_eli = p_eli,
        resp_cor = resp_cor,
        resp_cor_type = resp_cor_type,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "boin_etc_trials"
  )
}
