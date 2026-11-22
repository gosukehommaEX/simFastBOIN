#' Simulate iBOIN Trials
#'
#' @description
#'   Run the iBOIN dose-finding algorithm of Zhou et al. (2021) over many
#'   simulated trials and return the patient and DLT counts at every dose.
#'   This is the simulation engine behind \code{\link{sim_iboin}}, exposed
#'   for users who want the raw trial data.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param p_true
#'   Numeric vector. True DLT probability at each dose level, in increasing
#'   dose order.
#'
#' @param skeleton
#'   Numeric vector, one value per dose. Prior DLT probability at each dose,
#'   strictly between 0 and 1.
#'
#' @param prior_ess
#'   Whole number or vector of whole numbers, one per dose. Prior effective
#'   sample size. See \code{\link{iboin_boundary}}.
#'
#' @param n_cohort
#'   Integer scalar. Number of cohorts in a trial.
#'
#' @param cohort_size
#'   Integer scalar or vector. Number of patients per cohort. A scalar is
#'   used for every cohort. A vector shorter than \code{n_cohort} is padded
#'   with its last element and a longer one is truncated.
#'
#' @param n_trials
#'   Integer scalar. Number of trials to simulate. Defaults to 10000.
#'
#' @param start_dose
#'   Integer scalar. Dose level for the first cohort. Defaults to 1.
#'
#' @param n_earlystop
#'   Integer scalar. The trial stops once this many patients have been
#'   treated at the current dose and the design would stay there. Defaults
#'   to 18. The simulations of the article do not stop early; a value above
#'   the maximum sample size switches the rule off.
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
#'   Numeric scalar. Posterior probability cutoff for dose elimination.
#'   Defaults to 0.95.
#'
#' @param extrasafe
#'   Logical scalar. Apply the stricter safety stopping rule at the lowest
#'   dose. Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is
#'   relaxed for the safety stopping rule. Defaults to 0.05.
#'
#' @param robust
#'   Logical scalar. Use the robust prior, which drops the prior information
#'   above the prior MTD when the prior MTD lies in the upper half of the
#'   doses. Defaults to \code{FALSE}. See \code{\link{iboin_boundary}}.
#'
#' @param n_earlystop_rule
#'   Character scalar, either \code{"with_stay"} or \code{"simple"}. Under
#'   \code{"with_stay"} the trial stops at \code{n_earlystop} only when the
#'   design would remain at the current dose. Under \code{"simple"} it stops
#'   as soon as the threshold is reached.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. When \code{NULL} the current
#'   state of the random number generator is used. The state of the calling
#'   session is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{c("iboin_trials", "boin_trials")}, which is a
#'   list with components
#'   \item{n_pts}{Integer matrix, \code{n_trials} by number of doses, of
#'     patients treated.}
#'   \item{n_tox}{Integer matrix of the same shape, of DLTs observed.}
#'   \item{eliminated}{Logical matrix of the same shape, marking doses
#'     eliminated during the trial.}
#'   \item{cohorts_used}{Integer vector of cohorts completed in each trial.}
#'   \item{stop_reason}{Character vector giving why each trial ended.}
#'   \item{boundary}{The \code{iboin_boundary} object used.}
#'   \item{design}{The name of the design, \code{"iBOIN"}.}
#'   \item{settings}{List of the design parameters.}
#'
#' @details
#'   The trial runs as a BOIN trial, see \code{\link{boin_simulate}}, with
#'   the boundaries of \code{\link{iboin_boundary}}, which depend on the dose
#'   and on the number of patients treated there. Dose elimination and the
#'   safety stopping rule use a uniform prior, as in the BOIN design. The
#'   titration of the BOIN design is not available.
#'
#'   One uniform variate is drawn per patient, in enrollment order, and a
#'   cohort cut short by the maximum sample size still draws a full cohort's
#'   worth. With a prior effective sample size of 0 at every dose the trials
#'   are therefore, for the same seed, those of \code{\link{boin_simulate}}.
#'
#'   Possible values of \code{stop_reason} are \code{"lowest_dose_eliminated"},
#'   \code{"lowest_dose_too_toxic"}, \code{"n_earlystop"},
#'   \code{"max_sample_size"} and \code{"max_cohorts"}.
#'
#' @references
#'   Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
#'   Incorporating historical information to improve phase I clinical trials.
#'   Pharmaceutical Statistics, 20, 1017-1034. \doi{10.1002/pst.2121}
#'
#' @examples
#' trials <- iboin_simulate(
#'   target = 0.30,
#'   p_true = c(0.08, 0.15, 0.31, 0.45, 0.55),
#'   skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
#'   prior_ess = 3,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#' trials
#'
#' head(trials$n_pts)
#' table(trials$stop_reason)
#'
#' @seealso \code{\link{sim_iboin}}, \code{\link{iboin_boundary}},
#'   \code{\link{boin_select_mtd}}
#'
#' @export
iboin_simulate <- function(target, p_true, skeleton, prior_ess, n_cohort,
                           cohort_size, n_trials = 10000, start_dose = 1,
                           n_earlystop = 18, p_saf = NULL, p_tox = NULL,
                           cutoff_eli = 0.95, extrasafe = FALSE,
                           offset = 0.05, robust = FALSE,
                           n_earlystop_rule = c("with_stay", "simple"),
                           seed = 123) {

  n_earlystop_rule <- match.arg(n_earlystop_rule)

  check_p_true(p_true)
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(n_earlystop, "n_earlystop", 1L)
  check_count(start_dose, "start_dose", 1L)
  if (n_earlystop <= 6) {
    warning("'n_earlystop' is low; values between 9 and 18 are recommended",
            call. = FALSE)
  }

  n_doses <- length(p_true)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of doses in 'p_true'",
         call. = FALSE)
  }
  if (length(skeleton) != n_doses) {
    stop("'skeleton' must have one value for each dose of 'p_true'",
         call. = FALSE)
  }

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)

  bound <- iboin_boundary(target, skeleton, prior_ess, max_total_pts,
                          p_saf = p_saf, p_tox = p_tox,
                          cutoff_eli = cutoff_eli, extrasafe = extrasafe,
                          offset = offset, robust = robust)

  # The engine encodes a missing elimination boundary as zero.
  b_elim <- bound$b_elim
  b_elim[is.na(b_elim)] <- 0L

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- iboin_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    n_earlystop = as.integer(n_earlystop),
    early_stop_simple = identical(n_earlystop_rule, "simple"),
    extrasafe = extrasafe,
    target = as.numeric(target),
    cutoff_eli = as.numeric(cutoff_eli),
    offset = as.numeric(offset),
    b_esc = unname(bound$b_esc),
    b_deesc = unname(bound$b_deesc),
    b_elim = as.integer(b_elim),
    max_total_pts = as.integer(max_total_pts)
  )

  reasons <- c("lowest_dose_eliminated", "lowest_dose_too_toxic",
               "n_earlystop", "max_sample_size", "max_cohorts")

  dose_names <- paste0("DL", seq_len(n_doses))
  colnames(res$n_pts) <- dose_names
  colnames(res$n_tox) <- dose_names
  colnames(res$eliminated) <- dose_names

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$n_tox,
      eliminated = res$eliminated,
      cohorts_used = res$cohorts_used,
      stop_reason = reasons[res$stop_code + 1L],
      boundary = bound,
      design = "iBOIN",
      settings = list(
        target = target,
        p_true = p_true,
        skeleton = skeleton,
        prior_ess = bound$prior_ess,
        prior_ess_input = bound$prior_ess_input,
        robust = robust,
        n_cohort = n_cohort,
        cohort_size = cohort_size_vec,
        n_trials = n_trials,
        start_dose = start_dose,
        n_earlystop = n_earlystop,
        n_earlystop_rule = n_earlystop_rule,
        p_saf = bound$p_saf,
        p_tox = bound$p_tox,
        cutoff_eli = cutoff_eli,
        extrasafe = extrasafe,
        offset = offset,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = c("iboin_trials", "boin_trials")
  )
}
