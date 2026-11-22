#' Operating Characteristics of the gBOIN Design
#'
#' @description
#'   Simulate trials of the gBOIN design of Mu et al. (2019) for a toxicity
#'   outcome with several grades that each carry a score, and summarize how
#'   often each dose is selected as the MTD and how the patients are
#'   distributed over the doses.
#'
#' @inheritParams gboin_simulate
#'
#' @param min_mtd_sample
#'   Integer scalar. Smallest number of patients a dose must have been given
#'   to be selected as the MTD. Defaults to 1.
#'
#' @param overdose_cutoff
#'   Numeric scalar or \code{NULL}. Doses with a true mean normalized score
#'   above this value count as overdoses in the summary. Defaults to
#'   \code{NULL}, the target.
#'
#' @param keep_trials
#'   Logical scalar. Keep the trial level data in the result. Defaults to
#'   \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Report the progress in messages. Defaults to
#'   \code{FALSE}.
#'
#' @return
#'   An object of class \code{c("gboin_oc", "boin_oc")}, a list with the
#'   components of the result of \code{\link{sim_boin}}, in which
#'   \item{p_true}{is the true mean normalized score of each dose,}
#'   \item{n_tox_dose, total_n_tox}{are averages of the total normalized
#'     score, in the place of the number of DLTs,}
#'   and in addition
#'   \item{design}{The name of the design, \code{"gBOIN"}.}
#'   \item{p_grade, score}{The true grade probabilities and the scores.}
#'   \item{n_grade_dose}{Matrix, toxicity categories by doses, of the
#'     average number of patients in each category.}
#'   With \code{keep_trials = TRUE} the component \code{trials} also holds
#'   the array \code{n_grade} of \code{\link{gboin_simulate}}.
#'
#' @details
#'   The trials are simulated by \code{\link{gboin_simulate}}, which
#'   describes the design, and the MTD is selected by
#'   \code{\link{gboin_select_mtd}}. A trial stopped because the lowest dose
#'   is eliminated selects no dose.
#'
#'   With two toxicity categories of scores 0 and 1 the design is the BOIN
#'   design, and the result is, for the same seed, that of
#'   \code{\link{sim_boin}} without titration.
#'
#'   The result can be passed to \code{\link{oc_results}} and
#'   \code{\link{oc_compare}}, where it appears under the name gBOIN.
#'
#' @references
#'   Mu, R., Yuan, Y., Xu, J., Mandrekar, S. J. and Yin, J. (2019). gBOIN: a
#'   unified model-assisted phase I trial design accounting for toxicity
#'   grades, and binary or continuous end points. Journal of the Royal
#'   Statistical Society: Series C, 68, 289-308. \doi{10.1111/rssc.12263}
#'
#' @examples
#' # Scenario 1 of Mu et al. (2019): grades 0 or 1, 2, 3 and 4 at six doses,
#' # equivalent toxicity scores 0, 0.5, 1 and 1.5, target score 0.47
#' p_grade <- rbind(
#'   c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19),
#'   c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11),
#'   c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11),
#'   c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)
#' )
#' oc <- sim_gboin(
#'   target = 0.47 / 1.5,
#'   p_grade = p_grade,
#'   score = c(0, 0.5, 1, 1.5),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_earlystop = 100,
#'   n_trials = 500,
#'   seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{gboin_simulate}}, \code{\link{gboin_select_mtd}},
#'   \code{\link{sim_gboin_multi}}, \code{\link{sim_boin}}
#'
#' @export
sim_gboin <- function(target, p_grade, score, n_cohort, cohort_size,
                      n_trials = 10000, start_dose = 1, n_earlystop = 18,
                      p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                      extrasafe = FALSE, offset = 0.05,
                      elim_rule = c("posterior", "count"),
                      min_mtd_sample = 1, overdose_cutoff = NULL,
                      n_earlystop_rule = c("with_stay", "simple"),
                      keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  elim_rule <- match.arg(elim_rule)
  n_earlystop_rule <- match.arg(n_earlystop_rule)

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- gboin_simulate(
    target = target, p_grade = p_grade, score = score, n_cohort = n_cohort,
    cohort_size = cohort_size, n_trials = n_trials, start_dose = start_dose,
    n_earlystop = n_earlystop, p_saf = p_saf, p_tox = p_tox,
    cutoff_eli = cutoff_eli, extrasafe = extrasafe, offset = offset,
    elim_rule = elim_rule, n_earlystop_rule = n_earlystop_rule, seed = seed
  )

  if (verbose) message("Selecting the MTD ...")

  selection <- gboin_select_mtd(
    n_pts = trials$n_pts, score_sum = trials$n_tox, target = target,
    cutoff_eli = cutoff_eli, extrasafe = extrasafe, offset = offset,
    min_mtd_sample = min_mtd_sample
  )

  out <- mtd_oc_summary(trials, selection, trials$settings$p_true, target,
                        overdose_cutoff = overdose_cutoff,
                        keep_trials = keep_trials)

  if (verbose) message("Done.")

  out$lambda_e <- trials$lambda_e
  out$lambda_d <- trials$lambda_d
  out$design <- "gBOIN"
  out$p_grade <- p_grade
  out$score <- score
  n_grade_dose <- t(apply(trials$n_grade, c(2L, 3L), mean))
  dimnames(n_grade_dose) <- rev(dimnames(trials$n_grade)[2:3])
  out$n_grade_dose <- n_grade_dose
  if (keep_trials) out$trials$n_grade <- trials$n_grade
  out$call <- call_expr

  structure(out, class = c("gboin_oc", "boin_oc"))
}
