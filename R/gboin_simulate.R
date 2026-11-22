#' Simulate gBOIN Trials
#'
#' @description
#'   Run the gBOIN dose-finding algorithm of Mu et al. (2019) over many
#'   simulated trials, for a toxicity outcome with several grades that each
#'   carry a score, and return the patients, their grades and the total
#'   score at every dose. This is the simulation engine behind
#'   \code{\link{sim_gboin}}, exposed for users who want the raw trial data.
#'
#' @param target
#'   Numeric scalar. Target value of the mean normalized toxicity score, for
#'   example 0.313 for a target equivalent toxicity score of 0.47 with a
#'   largest score of 1.5. See the details.
#'
#' @param p_grade
#'   Numeric matrix. True probability of each toxicity category (rows, from
#'   the lowest to the highest) at each dose (columns, in increasing dose
#'   order). Every column sums to 1.
#'
#' @param score
#'   Numeric vector, one value per row of \code{p_grade}. Score of each
#'   toxicity category, for example \code{c(0, 0.5, 1, 1.5)} for grades 0 or
#'   1, 2, 3 and 4. The scores are divided by their maximum.
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
#'   to 18. A value above the maximum sample size switches the rule off.
#'
#' @param p_saf
#'   Numeric scalar. Highest mean score deemed subtherapeutic. Defaults to
#'   \code{0.6 * target}.
#'
#' @param p_tox
#'   Numeric scalar. Lowest mean score deemed overly toxic. Defaults to
#'   \code{1.4 * target}.
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
#' @param elim_rule
#'   Character scalar, either \code{"posterior"} or \code{"count"}. How a
#'   dose is eliminated. Under \code{"posterior"}, the default and the rule
#'   of the article, the posterior probability is evaluated at the total
#'   score. Under \code{"count"} the dose is eliminated when the total score
#'   reaches the elimination boundary of the number of DLTs of
#'   \code{\link{boin_boundary}}, a whole number, as in the
#'   'UnifiedDoseFinding' package. See the details.
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
#'   An object of class \code{c("gboin_trials", "boin_trials")}, which is a
#'   list with components
#'   \item{n_pts}{Integer matrix, \code{n_trials} by number of doses, of
#'     patients treated.}
#'   \item{n_tox}{Numeric matrix of the same shape: the total of the
#'     normalized scores of the patients treated at each dose (the total of
#'     the scores as given, divided by the largest score), which takes the
#'     place of the number of DLTs.}
#'   \item{n_grade}{Integer array, \code{n_trials} by doses by toxicity
#'     categories, of the patients in each category. The categories take
#'     the row names of \code{p_grade} when every row has one, and
#'     \code{G1}, \code{G2}, ... otherwise.}
#'   \item{eliminated}{Logical matrix marking doses eliminated during the
#'     trial.}
#'   \item{cohorts_used}{Integer vector of cohorts completed in each trial.}
#'   \item{stop_reason}{Character vector giving why each trial ended.}
#'   \item{lambda_e, lambda_d}{The escalation and de-escalation boundaries
#'     of the mean score.}
#'   \item{design}{The name of the design, \code{"gBOIN"}.}
#'   \item{settings}{List of the design parameters, with the true mean
#'     normalized score of each dose in \code{p_true}.}
#'
#' @details
#'   The gBOIN design generalizes the BOIN design to toxicity outcomes from
#'   the exponential family. This function covers the quasi-Bernoulli
#'   outcome: every toxicity category carries a score, such as the
#'   equivalent toxicity score (ETS) of Yuan et al. (2007), and the scores
#'   divided by their maximum, which lie between 0 and 1, are treated like
#'   fractional DLTs. With \code{x} the total normalized score of the
#'   \code{n} patients treated at the current dose, the design
#'   \itemize{
#'     \item escalates when \code{x / n <= lambda_e},
#'     \item de-escalates when \code{x / n >= lambda_d},
#'     \item stays otherwise,
#'   }
#'   where \code{lambda_e} and \code{lambda_d} are the BOIN boundaries of
#'   \code{\link{boin_lambda}} for \code{target}, \code{p_saf} and
#'   \code{p_tox}; they are the optimal boundaries for the quasi-Bernoulli
#'   outcome as well (Theorem 1 of the article). A dose with at least three
#'   patients is eliminated, together with every higher dose, when the
#'   posterior probability that its mean score exceeds \code{target}, from
#'   the distribution Beta(\code{x + 1}, \code{n - x + 1}), is above
#'   \code{cutoff_eli}; the next cohort is then treated at the dose below,
#'   and the trial stops when the lowest dose is eliminated.
#'
#'   \code{target}, \code{p_saf} and \code{p_tox} are on the normalized
#'   scale. The article computes its boundaries for the ETS example (0.373
#'   and 0.566) from the unnormalized target of 0.47; on the normalized
#'   scale the target is 0.47 / 1.5 = 0.313 and the boundaries, multiplied
#'   by 1.5, are 0.371 and 0.562.
#'
#'   The reference implementation in the 'UnifiedDoseFinding' package
#'   applies the elimination boundary of the number of DLTs, a whole number,
#'   to the total score, which eliminates a dose a little later than the
#'   posterior probability evaluated at the total score itself, the rule the
#'   article describes. \code{elim_rule = "count"} gives the rule of the
#'   package; it is meant for comparisons with it. For whole numbers of DLTs
#'   the two rules are the same.
#'
#'   One uniform variate is drawn per patient, in enrollment order: the
#'   patient falls in the highest category for which the probability of
#'   that category or a higher one exceeds the variate. A cohort cut short
#'   by the maximum sample size still draws a full cohort's worth. With two
#'   categories of scores 0 and 1 the trials are therefore, for the same
#'   seed, those of \code{\link{boin_simulate}} without titration.
#'
#'   The total normalized score of a dose is the total of the scores as
#'   given, divided once by the largest score, as in
#'   \code{\link{gboin_et_simulate}}. With scores that are binary fractions,
#'   such as 0, 0.5, 1 and 1.5, the total of the scores as given is exact, so
#'   that the same patients give the same total, and the same decision, on
#'   every platform.
#'
#'   Possible values of \code{stop_reason} are \code{"lowest_dose_eliminated"},
#'   \code{"lowest_dose_too_toxic"}, \code{"n_earlystop"},
#'   \code{"max_sample_size"} and \code{"max_cohorts"}.
#'
#' @references
#'   Mu, R., Yuan, Y., Xu, J., Mandrekar, S. J. and Yin, J. (2019). gBOIN: a
#'   unified model-assisted phase I trial design accounting for toxicity
#'   grades, and binary or continuous end points. Journal of the Royal
#'   Statistical Society: Series C, 68, 289-308. \doi{10.1111/rssc.12263}
#'
#'   Yuan, Z., Chappell, R. and Bailey, H. (2007). The continual reassessment
#'   method for multiple toxicity grades: a Bayesian quasi-likelihood
#'   approach. Biometrics, 63, 173-179.
#'
#' @examples
#' # Scenario 1 of Mu et al. (2019): grades 0 or 1, 2, 3 and 4 at six doses
#' p_grade <- rbind(
#'   c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19),
#'   c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11),
#'   c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11),
#'   c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)
#' )
#' trials <- gboin_simulate(
#'   target = 0.47 / 1.5,
#'   p_grade = p_grade,
#'   score = c(0, 0.5, 1, 1.5),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#' trials
#'
#' # Average number of patients of each grade at each dose
#' apply(trials$n_grade, c(3, 2), mean)
#'
#' @seealso \code{\link{sim_gboin}}, \code{\link{gboin_select_mtd}},
#'   \code{\link{boin_lambda}}
#'
#' @export
gboin_simulate <- function(target, p_grade, score, n_cohort, cohort_size,
                           n_trials = 10000, start_dose = 1,
                           n_earlystop = 18, p_saf = NULL, p_tox = NULL,
                           cutoff_eli = 0.95, extrasafe = FALSE,
                           offset = 0.05,
                           elim_rule = c("posterior", "count"),
                           n_earlystop_rule = c("with_stay", "simple"),
                           seed = 123) {

  elim_rule <- match.arg(elim_rule)
  n_earlystop_rule <- match.arg(n_earlystop_rule)

  grades <- gboin_check_grades(p_grade, score)
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(n_earlystop, "n_earlystop", 1L)
  check_count(start_dose, "start_dose", 1L)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_flag(extrasafe, "extrasafe")
  if (!is.numeric(offset) || length(offset) != 1L || !is.finite(offset) ||
      offset <= 0 || offset >= 0.5) {
    stop("'offset' must be a single number strictly between 0 and 0.5",
         call. = FALSE)
  }
  if (n_earlystop <= 6) {
    warning("'n_earlystop' is low; values between 9 and 18 are recommended",
            call. = FALSE)
  }

  n_doses <- ncol(p_grade)
  n_cat <- nrow(p_grade)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of doses in 'p_grade'",
         call. = FALSE)
  }

  if (is.null(p_saf)) p_saf <- 0.6 * target
  if (is.null(p_tox)) p_tox <- 1.4 * target
  lambda <- boin_lambda(target, p_saf = p_saf, p_tox = p_tox)

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)

  # Elimination boundary of the number of DLTs; a missing boundary is
  # encoded as zero.
  b_elim <- boin_boundary(target, max_total_pts, p_saf = p_saf,
                          p_tox = p_tox, cutoff_eli = cutoff_eli)$b_elim
  b_elim[is.na(b_elim)] <- 0L

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- gboin_simulate_cpp(
    n_trials = as.integer(n_trials),
    cum_prob = grades$cum_prob,
    score = as.numeric(score),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    n_earlystop = as.integer(n_earlystop),
    early_stop_simple = identical(n_earlystop_rule, "simple"),
    extrasafe = extrasafe,
    target = as.numeric(target),
    lambda_e = lambda$lambda_e,
    lambda_d = lambda$lambda_d,
    cutoff_eli = as.numeric(cutoff_eli),
    offset = as.numeric(offset),
    elim_by_count = identical(elim_rule, "count"),
    b_elim = as.integer(b_elim),
    max_total_pts = as.integer(max_total_pts)
  )

  reasons <- c("lowest_dose_eliminated", "lowest_dose_too_toxic",
               "n_earlystop", "max_sample_size", "max_cohorts")

  dose_names <- paste0("DL", seq_len(n_doses))
  # Row names are used when every row has one; rbind() of named and unnamed
  # vectors leaves some of them empty.
  grade_names <- rownames(p_grade)
  if (is.null(grade_names) || any(!nzchar(grade_names))) {
    grade_names <- paste0("G", seq_len(n_cat))
  }
  colnames(res$n_pts) <- dose_names
  colnames(res$score_sum) <- dose_names
  colnames(res$eliminated) <- dose_names
  dimnames(res$n_grade) <- list(NULL, dose_names, grade_names)

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$score_sum,
      n_grade = res$n_grade,
      eliminated = res$eliminated,
      cohorts_used = res$cohorts_used,
      stop_reason = reasons[res$stop_code + 1L],
      lambda_e = lambda$lambda_e,
      lambda_d = lambda$lambda_d,
      design = "gBOIN",
      settings = list(
        target = target,
        p_true = grades$p_true,
        p_grade = p_grade,
        score = score,
        score_norm = grades$score_norm,
        n_cohort = n_cohort,
        cohort_size = cohort_size_vec,
        n_trials = n_trials,
        start_dose = start_dose,
        n_earlystop = n_earlystop,
        n_earlystop_rule = n_earlystop_rule,
        p_saf = p_saf,
        p_tox = p_tox,
        cutoff_eli = cutoff_eli,
        extrasafe = extrasafe,
        offset = offset,
        elim_rule = elim_rule,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = c("gboin_trials", "boin_trials")
  )
}
