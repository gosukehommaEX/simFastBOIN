#' Simulate gBOIN-ET Trials
#'
#' @description
#'   Run the gBOIN-ET dose-finding algorithm of Takeda, Morita and Taguri
#'   (2022) over many simulated trials, for toxicity and efficacy outcomes
#'   with several categories that each carry a score, and return the
#'   patients, their categories and the total scores at every dose. This is
#'   the simulation engine behind \code{\link{sim_gboin_et}}, exposed for
#'   users who want the raw trial data.
#'
#' @inheritParams boin_et_simulate
#' @inheritParams boin_et_decision_table
#'
#' @param target
#'   Numeric scalar. Target value of the mean normalized toxicity score, for
#'   example \code{0.47 / 1.5} for a target equivalent toxicity score of 0.47
#'   with a largest score of 1.5.
#'
#' @param target_resp
#'   Numeric scalar. Target value of the mean normalized efficacy score, for
#'   example \code{0.95 / 1.5} for a target equivalent efficacy score of 0.95
#'   with a largest score of 1.5.
#'
#' @param p_grade
#'   Numeric matrix. True probability of each toxicity category (rows) at
#'   each dose (columns, in increasing dose order). Every column sums to 1.
#'   The first row is the category without a toxicity event.
#'
#' @param p_resp_grade
#'   Numeric matrix with one column per dose. True probability of each
#'   efficacy category (rows). Every column sums to 1. The first row is the
#'   category without an efficacy event.
#'
#' @param score
#'   Numeric vector, one value per row of \code{p_grade}. Score of each
#'   toxicity category, for example \code{c(0, 0.5, 1, 1.5)} for grades 0 or
#'   1, 2, 3 and 4. The scores are divided by their maximum.
#'
#' @param score_resp
#'   Numeric vector, one value per row of \code{p_resp_grade}. Score of each
#'   efficacy category, for example \code{c(0, 1, 1.5)} for stable or
#'   progressive disease, partial response and complete response. The scores
#'   are divided by their maximum.
#'
#' @param lambda1,lambda2,eta1
#'   Numeric scalars or \code{NULL}. Boundaries of the mean toxicity score
#'   and of the mean efficacy score. Defaults to \code{NULL}, the boundaries
#'   of \code{\link{boin_et_boundary}} for the targets. Give all three to
#'   use other values, such as those printed in an article.
#'
#' @param resp_cor
#'   Numeric scalar strictly between -1 and 1. Correlation of the Gaussian
#'   copula that links the toxicity event and the efficacy event of a
#'   patient. Defaults to 0.2, the default of the 'boinet' package.
#'
#' @return
#'   An object of class \code{boin_et_trials}, which is a list with components
#'   \item{n_pts}{Integer matrix, \code{n_trials} by number of doses, of
#'     patients treated.}
#'   \item{n_tox}{Numeric matrix of the same shape: the total of the
#'     normalized toxicity scores of the patients treated at each dose.}
#'   \item{n_resp}{Numeric matrix of the same shape: the total of the
#'     normalized efficacy scores.}
#'   \item{n_grade}{Integer array, \code{n_trials} by doses by toxicity
#'     categories, of the patients in each category.}
#'   \item{n_resp_grade}{Integer array, \code{n_trials} by doses by efficacy
#'     categories. The categories take the row names of \code{p_grade} and
#'     \code{p_resp_grade} when every row has one, and \code{G1}, \code{G2},
#'     ... and \code{R1}, \code{R2}, ... otherwise.}
#'   \item{cohorts_used}{Integer vector of cohorts enrolled in each trial.}
#'   \item{stop_reason}{Character vector giving why each trial ended.}
#'   \item{early_stop}{Logical vector, \code{TRUE} for the trials stopped
#'     early, which select no OBD.}
#'   \item{design}{\code{"gBOIN-ET"}.}
#'   \item{boundary}{The boundaries used, an object of class
#'     \code{boin_et_boundary}.}
#'   \item{settings}{List of the design parameters, with the true mean
#'     normalized scores of each dose in \code{p_true} and \code{p_resp}.}
#'
#' @details
#'   The gBOIN-ET design is the BOIN-ET design for ordinal outcomes. Every
#'   toxicity category and every efficacy category carries a score, such as
#'   the equivalent toxicity score and the equivalent efficacy score of the
#'   article. The scores divided by their maximum lie between 0 and 1 and
#'   are treated like fractional DLTs and responses (quasi-Bernoulli
#'   outcomes): the mean normalized scores of the patients at the current
#'   dose take the place of the observed DLT and response rates in the rules
#'   of \code{\link{boin_et_simulate}}. With \code{pt} and \code{pe} the mean
#'   toxicity and efficacy scores at the current dose, the design
#'   \itemize{
#'     \item escalates when \code{pt <= lambda1} and \code{pe <= eta1},
#'     \item stays when \code{pt < lambda2} and \code{pe > eta1},
#'     \item de-escalates when \code{pt >= lambda2},
#'     \item and otherwise moves to the dose with the highest mean efficacy
#'       score among the current dose and its neighbors, an untried higher
#'       dose being preferred and ties being broken at random.
#'   }
#'   A dose is excluded when the posterior probability that its mean toxicity
#'   score is below \code{p_eli} falls below \code{1 - cutoff_eli}, or when
#'   the posterior probability that its mean efficacy score is above
#'   \code{p_eli_resp} falls below \code{1 - cutoff_eli_resp}; the posterior
#'   distribution of a mean score with total \code{x} among \code{n} patients
#'   is Beta(\code{1 + x}, \code{1 + n - x}). The defaults are the rules of
#'   'boinet': \code{p_eli = target} with \code{cutoff_eli = 0.95}, the
#'   toxicity rule of Section 2.2 of Takeda et al. (2022), and
#'   \code{p_eli_resp = p_saf_resp} (by default \code{0.6 * target_resp})
#'   with \code{cutoff_eli_resp = 0.99}. Section 2.2 of the article excludes
#'   a dose for its efficacy when the posterior probability that its mean
#'   efficacy score is below half the target exceeds 0.50, which
#'   \code{p_eli_resp = target_resp / 2} and \code{cutoff_eli_resp = 0.50}
#'   give. How the next dose is found
#'   among the admissible doses, when the trial ends, and how the rules
#'   differ from those of 'boinet' unless \code{boinet_compatible} is
#'   \code{TRUE} is described in \code{\link{boin_et_simulate}}; with
#'   \code{boinet_compatible = TRUE} the rules are those of
#'   \code{boinet::gboinet()}.
#'
#'   \code{target}, \code{target_resp} and the other design parameters are on
#'   the normalized scale. The boundaries of \code{\link{boin_et_boundary}}
#'   come from a grid in steps of 0.01 that starts at \code{p_saf},
#'   \code{target} and \code{p_saf_resp}, so that they depend on the decimals
#'   of the targets: the targets \code{0.47 / 1.5} and \code{0.95 / 1.5} give
#'   0.1413, 0.3733 and 0.50, the values used by Takeda et al. (2022), who
#'   print 0.14, 0.37 and 0.50. The mean scores are not multiples of
#'   \code{1 / n}, so boundaries rounded to two decimals can give other
#'   decisions; \code{lambda1}, \code{lambda2} and \code{eta1} take the
#'   boundaries of an article directly.
#'
#'   A patient has a toxicity event when the toxicity category is not the
#'   first one, and an efficacy event likewise. The two events are linked by
#'   a Gaussian copula with correlation \code{resp_cor}, and the category of
#'   a patient with an event is drawn from the probabilities of the
#'   categories above the first, independently of everything else, which is
#'   the joint distribution of \code{boinet::gboinet()}. Each patient draws
#'   one uniform variate for the toxicity event, one for the efficacy event,
#'   and one more for each outcome with more than two categories; a tie
#'   under the rule for the neighboring doses draws one more. With two
#'   categories of scores 0 and 1 for both outcomes the trials are, for the
#'   same seed, those of \code{\link{boin_et_simulate}}. The 'boinet' package
#'   draws its random numbers differently, so the trials agree with it in
#'   distribution, not one by one.
#'
#'   The mean score of a dose is computed as the total of the scores as
#'   given, divided by the largest score times the number of patients. With
#'   scores that are binary fractions, such as 0, 0.5, 1 and 1.5, the total
#'   is exact, so that two doses whose mean scores are equal as fractions
#'   tie, and the trials are the same on every platform. Scores such as 0.1
#'   are not exact in floating point, and doses with equal mean scores may
#'   then fail to tie.
#'
#'   Possible values of \code{stop_reason} are those of
#'   \code{\link{boin_et_simulate}}.
#'
#' @references
#'   Takeda, K., Morita, S. and Taguri, M. (2022). gBOIN-ET: The generalized
#'   Bayesian optimal interval design for optimal dose-finding accounting for
#'   ordinal graded efficacy and toxicity in early clinical trials.
#'   Biometrical Journal, 64(7), 1178-1191. \doi{10.1002/bimj.202100263}
#'
#'   Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
#'   Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
#'   and Toxicity Outcomes. The R Journal, 17(2), 15-34.
#'   \doi{10.32614/RJ-2025-012}
#'
#' @examples
#' # Toxicity grades 0 or 1, 2, 3 and 4, and stable or progressive disease,
#' # partial response and complete response, at four doses
#' p_grade <- rbind(
#'   c(0.83, 0.62, 0.51, 0.34),
#'   c(0.12, 0.18, 0.19, 0.16),
#'   c(0.04, 0.11, 0.14, 0.15),
#'   c(0.01, 0.09, 0.16, 0.35)
#' )
#' p_resp_grade <- rbind(
#'   c(0.90, 0.50, 0.30, 0.25),
#'   c(0.05, 0.35, 0.20, 0.20),
#'   c(0.05, 0.15, 0.50, 0.55)
#' )
#' trials <- gboin_et_simulate(
#'   target = 0.47 / 1.5, target_resp = 0.95 / 1.5,
#'   p_grade = p_grade, p_resp_grade = p_resp_grade,
#'   score = c(0, 0.5, 1, 1.5), score_resp = c(0, 1, 1.5),
#'   n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 123
#' )
#' trials
#'
#' # Patients of the first trial by dose and toxicity category
#' trials$n_grade[1, , ]
#'
#' @seealso \code{\link{sim_gboin_et}}, \code{\link{gboin_et_select_obd}},
#'   \code{\link{boin_et_simulate}}
#'
#' @export
gboin_et_simulate <- function(target, target_resp, p_grade, p_resp_grade,
                              score, score_resp, n_cohort, cohort_size,
                              n_trials = 10000, start_dose = 1,
                              n_earlystop = NULL, p_saf = NULL, p_tox = NULL,
                              p_saf_resp = NULL, lambda1 = NULL,
                              lambda2 = NULL, eta1 = NULL, cutoff_eli = 0.95,
                              cutoff_eli_resp = 0.99, p_eli = NULL,
                              p_eli_resp = NULL, resp_cor = 0.2,
                              boinet_compatible = FALSE,
                              prior = rep(1 / 6, 6), seed = 123) {

  check_flag(boinet_compatible, "boinet_compatible")

  tox <- gboin_check_grades(p_grade, score)
  n_doses <- ncol(p_grade)
  resp <- gboin_et_check_resp(p_resp_grade, score_resp, n_doses)
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(start_dose, "start_dose", 1L)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of doses in 'p_grade'",
         call. = FALSE)
  }
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  check_scalar_prob(resp_cor, "resp_cor", lower = -1, upper = 1)

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)
  if (is.null(n_earlystop)) n_earlystop <- max_total_pts
  check_count(n_earlystop, "n_earlystop", 1L)

  bd <- gboin_et_boundary_arg(target, target_resp, p_saf, p_tox, p_saf_resp,
                              prior, lambda1, lambda2, eta1)
  eli <- boin_et_check_eli(p_eli, p_eli_resp, bd)

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  res <- gboin_et_simulate_cpp(
    n_trials = as.integer(n_trials),
    cum_tox = tox$cum_prob,
    cum_resp = resp$cum_prob,
    score_tox = as.numeric(score),
    score_resp = as.numeric(score_resp),
    resp_cor = as.numeric(resp_cor),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    n_earlystop = as.integer(n_earlystop),
    lambda1 = as.numeric(bd$lambda1),
    lambda2 = as.numeric(bd$lambda2),
    eta1 = as.numeric(bd$eta1),
    p_eli = as.numeric(eli$p_eli),
    p_eli_resp = as.numeric(eli$p_eli_resp),
    cutoff_eli = as.numeric(cutoff_eli),
    cutoff_eli_resp = as.numeric(cutoff_eli_resp),
    boinet_compatible = as.integer(boinet_compatible)
  )

  reasons <- boin_et_stop_reasons()

  dose_names <- paste0("DL", seq_len(n_doses))
  # Row names are used when every row has one; rbind() of named and unnamed
  # vectors leaves some of them empty.
  grade_names <- rownames(p_grade)
  if (is.null(grade_names) || any(!nzchar(grade_names))) {
    grade_names <- paste0("G", seq_len(nrow(p_grade)))
  }
  resp_names <- rownames(p_resp_grade)
  if (is.null(resp_names) || any(!nzchar(resp_names))) {
    resp_names <- paste0("R", seq_len(nrow(p_resp_grade)))
  }
  colnames(res$n_pts) <- dose_names
  colnames(res$tox_sum) <- dose_names
  colnames(res$resp_sum) <- dose_names
  dimnames(res$n_tox_cat) <- list(NULL, dose_names, grade_names)
  dimnames(res$n_resp_cat) <- list(NULL, dose_names, resp_names)

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$tox_sum,
      n_resp = res$resp_sum,
      n_grade = res$n_tox_cat,
      n_resp_grade = res$n_resp_cat,
      cohorts_used = res$cohorts_used,
      stop_reason = reasons[res$stop_code + 1L],
      early_stop = res$stop_code <= 1L,
      design = "gBOIN-ET",
      boundary = bd,
      settings = list(
        target = bd$target,
        target_resp = bd$target_resp,
        p_true = tox$p_true,
        p_resp = resp$p_resp,
        p_grade = p_grade,
        p_resp_grade = p_resp_grade,
        score = score,
        score_resp = score_resp,
        score_norm = tox$score_norm,
        score_resp_norm = resp$score_norm,
        n_cohort = n_cohort,
        cohort_size = cohort_size_vec,
        n_trials = n_trials,
        start_dose = start_dose,
        n_earlystop = n_earlystop,
        p_saf = bd$p_saf,
        p_tox = bd$p_tox,
        p_saf_resp = bd$p_saf_resp,
        cutoff_eli = cutoff_eli,
        cutoff_eli_resp = cutoff_eli_resp,
        p_eli = eli$p_eli,
        p_eli_resp = eli$p_eli_resp,
        resp_cor = resp_cor,
        boinet_compatible = boinet_compatible,
        prior = prior,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "boin_et_trials"
  )
}
