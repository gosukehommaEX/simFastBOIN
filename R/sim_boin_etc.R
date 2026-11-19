#' Operating Characteristics of the BOIN-ETC Design
#'
#' @description
#'   Simulate trials of the BOIN-ETC design of Kakizume et al. (2024) for two
#'   drugs given together with \code{\link{boin_etc_simulate}}, select the
#'   optimal dose combination (ODC) of every trial with
#'   \code{\link{boin_etc_select_obd}}, and summarize the selection
#'   percentages and the numbers of patients, DLTs and responses at every
#'   combination.
#'
#' @inheritParams boin_etc_simulate
#' @inheritParams boin_etc_select_obd
#'
#' @param odc_tie
#'   Character scalar. \code{"lowest"} (the default) counts one ODC per
#'   trial, the tied candidate with the lowest sum of the levels;
#'   \code{"all"} counts every tied candidate as selected, the definition of
#'   the published percentages of Kakizume et al. (2024), so that
#'   \code{sel_percent} can add up to more than 100.
#'
#' @param keep_trials
#'   Logical scalar. Keep the simulated trials (\code{trials}) and the
#'   selection (\code{obd}, \code{tied}) in the result. Defaults to
#'   \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Print progress messages. Defaults to \code{FALSE}.
#'
#' @return
#'   An object of class \code{boin_etc_oc}, which is a list with components
#'   \item{design}{\code{"BOIN-ETC1"}, \code{"BOIN-ETC2"} or
#'     \code{"BOIN-ETC3"}.}
#'   \item{target, target_resp, p_true, p_resp}{The targets and the true
#'     probabilities.}
#'   \item{n_rows, n_cols, n_trials}{Numbers of levels of drug A and drug B,
#'     and of trials.}
#'   \item{lambda1, lambda2, eta1}{The BOIN-ET boundaries.}
#'   \item{sel_percent}{Matrix of the percentage of trials selecting each
#'     combination as the ODC.}
#'   \item{percent_no_obd, percent_early_stop}{Percentages of trials
#'     selecting no ODC and stopped early.}
#'   \item{n_pts_dose, n_tox_dose, n_resp_dose}{Matrices of the average
#'     numbers of patients, DLTs and responses at each combination.}
#'   \item{total_n_pts, total_n_tox, total_n_resp}{Their totals.}
#'   \item{mean_subtrials}{Average number of subtrials run.}
#'   \item{true_obd}{Two-column matrix of the true optimal combinations: the
#'     combinations with the highest true utility among those with a DLT
#'     probability at most \code{target}.}
#'   \item{pct_true_obd_selected, pct_patients_true_obd}{Percentage of
#'     trials selecting a true optimal combination (with
#'     \code{odc_tie = "all"}, the sum of \code{sel_percent} over them) and
#'     percentage of the patients treated there.}
#'   \item{overdose}{List with the \code{cutoff} (\code{target}), the
#'     \code{combinations} above it, and the percentages of trials selecting
#'     one of them (\code{pct_selected}) and of patients treated there
#'     (\code{pct_patients}).}
#'   \item{stop_reason_percent}{Percentage of trials by stopping reason.}
#'   \item{odc_tie, settings, call}{The settings used.}
#'   \item{trials, obd, tied}{With \code{keep_trials = TRUE}, the
#'     \code{boin_etc_trials} object, a two-column matrix of the selected
#'     combinations and the logical array of the tied candidates.}
#'
#' @details
#'   The percentages of patients are those of the average numbers of
#'   patients. Kakizume et al. (2024) report the average over trials of the
#'   percentage of the patients of a trial treated at the optimal
#'   combinations (%nODC); with \code{keep_trials = TRUE} it can be computed
#'   from \code{trials}, as in \code{inst/reproduce/reproduce_Kakizume2024.R}.
#'   Their optimal combinations (bold in their Table 2) are not always those
#'   with the highest true utility; \code{\link{oc_compare}} accepts them
#'   through its argument \code{truth}.
#'
#' @references
#'   Kakizume, T., Takeda, K., Taguri, M. and Morita, S. (2024). BOIN-ETC: A
#'   Bayesian optimal interval design considering efficacy and toxicity to
#'   identify the optimal dose combinations. Statistical Methods in Medical
#'   Research, 33(4), 716-727. \doi{10.1177/09622802241236936}
#'
#' @examples
#' # Scenario 1 of Kakizume et al. (2024), BOIN-ETC1 with n_stop = 6
#' p_true <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
#' p_resp <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
#' oc <- sim_boin_etc(
#'   target = 0.30, target_resp = 0.60, p_true = p_true, p_resp = p_resp,
#'   n_cohort = 8, cohort_size = 3, n_trials = 200, seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_boin_etc_multi}}, \code{\link{boin_etc_simulate}},
#'   \code{\link{boin_etc_select_obd}}, \code{\link{print.boin_etc_oc}}
#'
#' @export
sim_boin_etc <- function(target, target_resp, p_true, p_resp, n_cohort,
                         cohort_size, n_trials = 10000,
                         subtrials = c("etc1", "etc2", "etc3"),
                         n_earlystop = 6,
                         etc1_start = c("diagonal", "below"),
                         rule4 = c("paper", "authors_code"),
                         p_saf = NULL, p_tox = NULL, p_saf_resp = NULL,
                         prior = rep(1 / 6, 6), cutoff_eli = 0.95,
                         p_eli = NULL, resp_cor = 0,
                         resp_cor_type = c("copula", "binary"),
                         psi00 = 40, psi11 = 60, utility_prior = c(1, 1),
                         odc_tie = c("lowest", "all"), keep_trials = FALSE,
                         verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  subtrials <- match.arg(subtrials)
  etc1_start <- match.arg(etc1_start)
  rule4 <- match.arg(rule4)
  resp_cor_type <- match.arg(resp_cor_type)
  odc_tie <- match.arg(odc_tie)
  boin_etc_check_utility(psi00, psi11, utility_prior)
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- boin_etc_simulate(
    target = target, target_resp = target_resp, p_true = p_true,
    p_resp = p_resp, n_cohort = n_cohort, cohort_size = cohort_size,
    n_trials = n_trials, subtrials = subtrials, n_earlystop = n_earlystop,
    etc1_start = etc1_start, rule4 = rule4, p_saf = p_saf, p_tox = p_tox,
    p_saf_resp = p_saf_resp, prior = prior, cutoff_eli = cutoff_eli,
    p_eli = p_eli, resp_cor = resp_cor, resp_cor_type = resp_cor_type,
    seed = seed
  )

  if (verbose) message("Selecting the optimal dose combination ...")

  selection <- boin_etc_select_obd(
    n_pts = trials$n_pts, n_tox = trials$n_tox, n_resp = trials$n_resp,
    n_tox_resp = trials$n_tox_resp, target = target, psi00 = psi00,
    psi11 = psi11, utility_prior = utility_prior, odc_tie = odc_tie,
    excluded = trials$excluded, early_stop = trials$early_stop
  )

  boin_etc_summarize(trials, selection, psi00, psi11, utility_prior, odc_tie,
                     keep_trials, call_expr)
}
