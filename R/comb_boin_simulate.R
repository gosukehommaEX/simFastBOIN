#' Simulate BOIN Drug-Combination Trials
#'
#' @description
#'   Run the BOIN drug-combination design of Lin and Yin (2017) over many
#'   simulated trials of two drugs given together, and return the patient and
#'   DLT counts at every dose combination. This is the simulation engine behind
#'   \code{\link{sim_comb_boin}}, exposed for users who want the raw trial data.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param p_true
#'   Numeric matrix. True DLT probability of each dose combination, with one row
#'   per dose level of drug A and one column per dose level of drug B, both in
#'   increasing dose order.
#'
#' @param n_cohort
#'   Integer scalar. Number of cohorts in a trial.
#'
#' @param cohort_size
#'   Integer scalar or vector. Number of patients per cohort. A scalar is used for
#'   every cohort. A vector shorter than \code{n_cohort} is padded with its last
#'   element and a longer one is truncated.
#'
#' @param n_trials
#'   Integer scalar. Number of trials to simulate. Defaults to 10000.
#'
#' @param start_dose
#'   Integer vector of length two. Dose levels of drug A and drug B for the first
#'   cohort. Defaults to \code{c(1, 1)}. Titration, when requested, also starts
#'   here.
#'
#' @param n_earlystop
#'   Integer scalar. The trial stops once this many patients have been treated at
#'   the current combination and the design would neither escalate nor
#'   de-escalate, or cannot. Defaults to 18. Set it to a value above the maximum
#'   sample size to switch this rule off.
#'
#' @param p_saf
#'   Numeric scalar. Highest DLT probability deemed subtherapeutic.
#'   Defaults to \code{0.6 * target}.
#'
#' @param p_tox
#'   Numeric scalar. Lowest DLT probability deemed overly toxic.
#'   Defaults to \code{1.4 * target}.
#'
#' @param cutoff_eli
#'   Numeric scalar. Posterior probability cutoff for eliminating a combination.
#'   Defaults to 0.95.
#'
#' @param extrasafe
#'   Logical scalar. Apply the stricter safety stopping rule at the lowest
#'   combination. Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is
#'   relaxed for the safety stopping rule. Defaults to 0.05.
#'
#' @param titration
#'   Logical scalar. Treat one patient per combination, moving up one level of
#'   either drug, until the first DLT, before switching to full cohorts. Ignored
#'   when the first cohort size is one. Defaults to \code{FALSE}.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. When \code{NULL} the current
#'   state of the random number generator is used. The state of the calling
#'   session is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{comb_boin_trials}, which is a list with components
#'   \item{n_pts}{Integer array of patients treated, with dimensions trial by level of drug A by level of drug B.}
#'   \item{n_tox}{Integer array of the same shape, of DLTs observed.}
#'   \item{eliminated}{Logical array of the same shape, marking combinations eliminated during the trial.}
#'   \item{cohorts_used}{Integer vector of cohorts completed in each trial.}
#'   \item{stop_reason}{Character vector giving why each trial ended.}
#'   \item{boundary}{The \code{boin_boundary} object used.}
#'   \item{settings}{List of the design parameters.}
#'
#' @details
#'   The decision boundaries are those of the single-agent design, from
#'   \code{\link{boin_boundary}}. After each cohort the observed DLT rate at the
#'   current combination is compared with them. When the design escalates, the
#'   next cohort goes one level up of drug A or of drug B, whichever has the
#'   larger posterior probability, under a Beta(0.5, 0.5) prior, that its DLT
#'   probability lies between \code{lambda_e} and \code{lambda_d}; de-escalation
#'   chooses one level down in the same way. Eliminated combinations are not
#'   escalated to. When a combination is eliminated, so is every combination
#'   with at least the same level of both drugs, and the trial stops when the
#'   lowest combination is eliminated.
#'
#'   The rules, and the order in which they are applied, are those of
#'   \code{BOIN::get.oc.comb()} with \code{mtd.contour = FALSE}. That
#'   implementation adds two things to the article: a candidate's posterior
#'   probability is increased by 0.0005 for every patient already treated
#'   there, so that ties go to the better informed combination, and the trial
#'   can stop at \code{n_earlystop}. Random variates are
#'   drawn in the same order as there: one per patient, one for the direction
#'   of every titration step that can go either way, and one for every
#'   escalation or de-escalation that has a candidate, even when there is only
#'   one. With the same seed and matching arguments a run therefore reproduces
#'   the reference implementation trial by trial.
#'
#'   During titration the patients treated before the first DLT are added to
#'   the planned cohorts, as in the reference implementation: the first full
#'   cohort only tops up the combination reached. A trial can therefore treat a
#'   few more patients than \code{n_cohort} cohorts hold.
#'
#'   Possible values of \code{stop_reason} are \code{"lowest_dose_eliminated"},
#'   \code{"lowest_dose_too_toxic"}, \code{"n_earlystop"} and
#'   \code{"max_cohorts"}.
#'
#' @references
#'   Lin, R. and Yin, G. (2017). Bayesian Optimal Interval Design for Dose
#'   Finding in Drug-Combination Trials. Statistical Methods in Medical
#'   Research, 26(5), 2155-2167.
#'
#'   Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020). BOIN: An R
#'   Package for Designing Single-Agent and Drug-Combination Dose-Finding Trials Using
#'   Bayesian Optimal Interval Designs. Journal of Statistical Software, 94(13), 1-32.
#'
#' @examples
#' p_true <- matrix(c(0.02, 0.04, 0.08, 0.14,
#'                    0.08, 0.25, 0.42, 0.48,
#'                    0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE)
#'
#' trials <- comb_boin_simulate(
#'   target = 0.25,
#'   p_true = p_true,
#'   n_cohort = 16,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#' trials
#'
#' # Patients treated at each combination in the first trial
#' trials$n_pts[1, , ]
#'
#' @seealso \code{\link{sim_comb_boin}}, \code{\link{comb_boin_select_mtd}}
#'
#' @export
comb_boin_simulate <- function(target, p_true, n_cohort, cohort_size,
                               n_trials = 10000, start_dose = c(1, 1),
                               n_earlystop = 18, p_saf = NULL, p_tox = NULL,
                               cutoff_eli = 0.95, extrasafe = FALSE,
                               offset = 0.05, titration = FALSE, seed = 123) {

  check_p_true_matrix(p_true)
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(n_earlystop, "n_earlystop", 1L)
  check_flag(titration, "titration")
  check_start_combination(start_dose, p_true)
  if (n_earlystop <= 6) {
    warning("'n_earlystop' is low; values between 9 and 18 are recommended",
            call. = FALSE)
  }

  n_rows <- nrow(p_true)
  n_cols <- ncol(p_true)

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)

  # Titration needs a cohort size above one, otherwise it is the design itself.
  if (titration && cohort_size_vec[1L] == 1L) titration <- FALSE

  # No combination can receive more patients than the planned cohorts hold,
  # titration included, because titration treats one patient per combination
  # and the first cohort only tops the last one up.
  max_total_pts <- sum(cohort_size_vec)

  bound <- boin_boundary(target, max_total_pts, p_saf = p_saf, p_tox = p_tox,
                         cutoff_eli = cutoff_eli, extrasafe = extrasafe,
                         offset = offset)

  # The engine encodes a missing elimination boundary as zero.
  b_elim <- bound$b_elim
  b_elim[is.na(b_elim)] <- 0L

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- comb_boin_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    n_rows = as.integer(n_rows),
    n_cols = as.integer(n_cols),
    cohort_size = cohort_size_vec,
    start_row = as.integer(start_dose[1L]),
    start_col = as.integer(start_dose[2L]),
    n_earlystop = as.integer(n_earlystop),
    titration = titration,
    extrasafe = extrasafe,
    target = as.numeric(target),
    cutoff_eli = as.numeric(cutoff_eli),
    offset = as.numeric(offset),
    lambda_e = bound$lambda_e,
    lambda_d = bound$lambda_d,
    b_esc = bound$b_esc,
    b_deesc = bound$b_deesc,
    b_elim = as.integer(b_elim)
  )

  reasons <- c("lowest_dose_eliminated", "lowest_dose_too_toxic",
               "n_earlystop", "max_cohorts")

  trial_names <- list(NULL, paste0("A", seq_len(n_rows)),
                      paste0("B", seq_len(n_cols)))
  dimnames(res$n_pts) <- trial_names
  dimnames(res$n_tox) <- trial_names
  dimnames(res$eliminated) <- trial_names

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$n_tox,
      eliminated = res$eliminated,
      cohorts_used = res$cohorts_used,
      stop_reason = reasons[res$stop_code + 1L],
      boundary = bound,
      settings = list(
        target = target,
        p_true = p_true,
        n_cohort = n_cohort,
        cohort_size = cohort_size_vec,
        n_trials = n_trials,
        start_dose = as.integer(start_dose),
        n_earlystop = n_earlystop,
        p_saf = bound$p_saf,
        p_tox = bound$p_tox,
        cutoff_eli = cutoff_eli,
        extrasafe = extrasafe,
        offset = offset,
        titration = titration,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "comb_boin_trials"
  )
}
