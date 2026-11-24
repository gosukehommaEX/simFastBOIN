#' Operating Characteristics of the Waterfall Drug-Combination Design
#'
#' @description
#'   Simulate a trial of two drugs given together, run with the waterfall design
#'   of Zhang and Yuan (2016), many times under one matrix of true DLT
#'   probabilities, and summarize how often each dose combination is selected
#'   as an MTD of the contour, how often the whole contour is selected, how
#'   many patients are treated at each combination and how many DLTs are
#'   observed.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param p_true
#'   Numeric matrix. True DLT probability of each dose combination, with one row
#'   per dose level of drug A and one column per dose level of drug B, both in
#'   increasing dose order. It needs at least two rows and no more rows than
#'   columns.
#'
#' @param cohort_size
#'   Integer scalar. Number of patients per cohort.
#'
#' @param n_cohort
#'   Integer vector with one element per dose level of drug A, or \code{NULL}.
#'   Maximum number of cohorts of the first, second, ... subtrial. Defaults to
#'   \code{NULL}, which gives every subtrial about four patients per dose of
#'   the subtrial (see Details).
#'
#' @param n_trials
#'   Integer scalar. Number of trials to simulate. Defaults to 10000.
#'
#' @param n_earlystop
#'   Integer scalar. A subtrial ends once this many patients have been treated
#'   at the dose chosen for its next cohort. Defaults to 12.
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
#'   Logical scalar. Apply the stricter safety stopping rule at the lowest dose
#'   of every subtrial. Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is
#'   relaxed for the safety stopping rule. Defaults to 0.05.
#'
#' @param titration
#'   Logical scalar. Start the first subtrial with one patient per dose, up to
#'   the first DLT. Ignored when \code{cohort_size} is one.
#'   Defaults to \code{FALSE}.
#'
#' @param bound_mtd
#'   Logical scalar. Require the estimate at the candidate MTD of a subtrial to
#'   be below the de-escalation boundary. Defaults to \code{FALSE}.
#'
#' @param boin_compatible
#'   Logical scalar. Drop the data of a subtrial that ends without a candidate
#'   MTD, as \code{BOIN::get.oc.comb()} does (see Details).
#'   Defaults to \code{FALSE}.
#'
#' @param keep_trials
#'   Logical scalar. Keep the full trial by trial data in the result.
#'   Defaults to \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Report progress while the simulation runs.
#'   Defaults to \code{FALSE}.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. The state of the calling session
#'   is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{waterfall_oc}, which is a list with components
#'   \item{sel_percent}{Matrix of the percentage of trials selecting each combination as the MTD of its level of drug A. A trial can select one combination per level of drug A, so the percentages can add up to more than 100.}
#'   \item{percent_no_mtd}{Percentage of trials selecting no MTD at any level of drug A.}
#'   \item{n_pts_dose}{Matrix of the average number of patients treated at each combination.}
#'   \item{n_tox_dose}{Matrix of the average number of DLTs observed at each combination.}
#'   \item{total_n_pts}{Average total number of patients per trial.}
#'   \item{total_n_tox}{Average total number of DLTs per trial.}
#'   \item{total_n_enrolled}{Average number of patients treated per trial, including those of subtrials whose data were dropped; it differs from \code{total_n_pts} only with \code{boin_compatible = TRUE}.}
#'   \item{true_mtd}{Two-column matrix of the dose levels of drug A and drug B of the true MTD contour.}
#'   \item{pct_contour_selected}{Percentage of trials selecting exactly the true MTD contour; \code{NA} when no level of drug A has a true MTD.}
#'   \item{pct_patients_contour}{Percentage of all patients treated on the true MTD contour.}
#'   \item{pct_patients_above}{Percentage of all patients treated above it.}
#'   \item{pct_patients_below}{Percentage of all patients treated below it.}
#'   \item{avg_n_subtrials}{Average number of subtrials per trial.}
#'   \item{stop_reason_percent}{Percentage of trials by reason for ending.}
#'   \item{trials}{Trial level data when \code{keep_trials} is \code{TRUE}, otherwise \code{NULL}; \code{mtd} is a matrix with one row per trial and one column per level of drug A, holding the level of drug B of the MTD.}
#'   together with the design parameters and the call.
#'
#' @details
#'   The trials are simulated with \code{\link{waterfall_simulate}} and the MTD
#'   contour is selected with \code{\link{waterfall_select_mtd}}; the rules,
#'   and the points in which they differ from the article and from the
#'   reference implementation, are described there. With
#'   \code{boin_compatible = TRUE}, the same seed and matching arguments, the
#'   results equal those of \code{BOIN::get.oc.comb()} with
#'   \code{mtd.contour = TRUE}; with a seed of 6, 1000 trials and the default
#'   \code{n_cohort}, they are the results of the article.
#'
#'   The true MTD of a level of drug A is the combination of that row whose
#'   DLT probability is closest to \code{target}, provided it does not exceed
#'   \code{target + 0.05}. A combination is below the contour when both of its
#'   dose levels are at most those of some true MTD, and above it when both are
#'   at least those of some true MTD; every other combination, including the
#'   true MTDs, counts as on the contour. These are the definitions of the
#'   reference implementation. The three percentages of patients are \code{NA}
#'   when no level of drug A has a true MTD.
#'
#'   The phase I/II extension of the article is not provided.
#'
#' @references
#'   Zhang, L. and Yuan, Y. (2016). A Practical Bayesian Design to Identify the
#'   Maximum Tolerated Dose Contour for Drug Combination Trials. Statistics in
#'   Medicine, 35(27), 4924-4936.
#'
#'   Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020). BOIN: An R
#'   Package for Designing Single-Agent and Drug-Combination Dose-Finding Trials Using
#'   Bayesian Optimal Interval Designs. Journal of Statistical Software, 94(13), 1-32.
#'
#' @examples
#' # Scenario 1 of Zhang and Yuan (2016): the MTD contour is A1B3 and A2B2
#' p_true <- matrix(c(0.03, 0.10, 0.28,
#'                    0.10, 0.30, 0.50), nrow = 2, byrow = TRUE)
#'
#' oc <- sim_waterfall(
#'   target = 0.30,
#'   p_true = p_true,
#'   cohort_size = 3,
#'   n_trials = 500,
#'   seed = 123
#' )
#' oc
#'
#' \donttest{
#' # The settings of the article
#' oc_article <- sim_waterfall(
#'   target = 0.30,
#'   p_true = p_true,
#'   cohort_size = 3,
#'   n_trials = 1000,
#'   boin_compatible = TRUE,
#'   seed = 6
#' )
#' oc_article$sel_percent
#' }
#'
#' @seealso \code{\link{sim_waterfall_multi}}, \code{\link{waterfall_simulate}},
#'   \code{\link{waterfall_select_mtd}}, \code{\link{sim_comb_boin}}
#'
#' @export
sim_waterfall <- function(target, p_true, cohort_size, n_cohort = NULL,
                          n_trials = 10000, n_earlystop = 12, p_saf = NULL,
                          p_tox = NULL, cutoff_eli = 0.95, extrasafe = FALSE,
                          offset = 0.05, titration = FALSE, bound_mtd = FALSE,
                          boin_compatible = FALSE, keep_trials = FALSE,
                          verbose = FALSE, seed = 123) {

  call_expr <- match.call()

  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- waterfall_simulate(
    target = target, p_true = p_true, cohort_size = cohort_size,
    n_cohort = n_cohort, n_trials = n_trials, n_earlystop = n_earlystop,
    p_saf = p_saf, p_tox = p_tox, cutoff_eli = cutoff_eli,
    extrasafe = extrasafe, offset = offset, titration = titration,
    bound_mtd = bound_mtd, boin_compatible = boin_compatible, seed = seed
  )

  if (verbose) message("Selecting the MTD contour ...")

  selection <- waterfall_select_mtd(
    n_pts = trials$n_pts, n_tox = trials$n_tox,
    eliminated = trials$eliminated, target = target
  )

  n_rows <- nrow(p_true)
  n_cols <- ncol(p_true)
  labels <- list(paste0("A", seq_len(n_rows)), paste0("B", seq_len(n_cols)))

  # One row per trial, one column per level of drug A.
  mtd <- matrix(selection$dose_b, nrow = n_trials, ncol = n_rows,
                byrow = TRUE)
  colnames(mtd) <- labels[[1L]]

  sel_count <- matrix(0, n_rows, n_cols, dimnames = labels)
  for (a in seq_len(n_rows)) {
    sel_count[a, ] <- tabulate(mtd[, a], nbins = n_cols)
  }
  # The summaries are computed as the reference implementation computes
  # them, from sums over the trials and with its order of operations, so
  # that the two agree to the last digit.
  sel_percent <- sel_count * 100 / n_trials
  percent_no_mtd <- 100 * sum(rowSums(!is.na(mtd)) == 0) / n_trials

  n_pts_num <- array(as.numeric(trials$n_pts), dim(trials$n_pts))
  n_tox_num <- array(as.numeric(trials$n_tox), dim(trials$n_tox))
  n_pts_sum <- matrix(apply(n_pts_num, c(2L, 3L), sum), n_rows, n_cols,
                      dimnames = labels)
  n_tox_sum <- matrix(apply(n_tox_num, c(2L, 3L), sum), n_rows, n_cols,
                      dimnames = labels)
  n_pts_dose <- n_pts_sum / n_trials
  n_tox_dose <- n_tox_sum / n_trials

  contour <- waterfall_contour(p_true, target)
  has_contour <- nrow(contour$true_mtd) > 0L
  if (has_contour) {
    true_col <- matrix(contour$true_col, nrow = n_trials, ncol = n_rows,
                       byrow = TRUE)
    both_na <- is.na(mtd) & is.na(true_col)
    equal <- !is.na(mtd) & !is.na(true_col) & mtd == true_col
    correct <- rowSums(both_na | equal) == n_rows
    pct_contour_selected <- 100 * sum(correct) / n_trials
    total_pts <- sum(n_pts_sum)
    pct_patients_contour <- 100 * sum(n_pts_sum[contour$at]) / total_pts
    pct_patients_above <- 100 * sum(n_pts_sum[contour$above]) / total_pts
    pct_patients_below <- 100 * sum(n_pts_sum[contour$below]) / total_pts
  } else {
    pct_contour_selected <- NA_real_
    pct_patients_contour <- NA_real_
    pct_patients_above <- NA_real_
    pct_patients_below <- NA_real_
  }

  stop_reason_percent <- 100 * table(trials$stop_reason) / n_trials

  if (verbose) message("Done.")

  structure(
    list(
      sel_percent = sel_percent,
      percent_no_mtd = percent_no_mtd,
      n_pts_dose = n_pts_dose,
      n_tox_dose = n_tox_dose,
      total_n_pts = sum(n_pts_dose),
      total_n_tox = sum(n_tox_dose),
      total_n_enrolled = sum(trials$n_enrolled) / n_trials,
      true_mtd = contour$true_mtd,
      pct_contour_selected = pct_contour_selected,
      pct_patients_contour = pct_patients_contour,
      pct_patients_above = pct_patients_above,
      pct_patients_below = pct_patients_below,
      avg_n_subtrials = mean(trials$n_subtrials),
      stop_reason_percent = stop_reason_percent,
      target = target,
      p_true = p_true,
      p_saf = trials$settings$p_saf,
      p_tox = trials$settings$p_tox,
      lambda_e = trials$boundary$lambda_e,
      lambda_d = trials$boundary$lambda_d,
      n_rows = n_rows,
      n_cols = n_cols,
      n_trials = n_trials,
      settings = trials$settings,
      trials = if (keep_trials) {
        list(
          n_pts = trials$n_pts,
          n_tox = trials$n_tox,
          eliminated = trials$eliminated,
          n_subtrials = trials$n_subtrials,
          n_enrolled = trials$n_enrolled,
          stop_reason = trials$stop_reason,
          mtd = mtd
        )
      } else {
        NULL
      },
      call = call_expr
    ),
    class = "waterfall_oc"
  )
}
