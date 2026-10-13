#' Operating Characteristics of the BOIN Drug-Combination Design
#'
#' @description
#'   Simulate a trial of two drugs given together, run with the BOIN
#'   drug-combination design of Lin and Yin (2017), many times under one matrix
#'   of true DLT probabilities, and summarize how often each dose combination is
#'   selected as the MTD, how many patients are treated at each combination and
#'   how many DLTs are observed.
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
#'   cohort. Defaults to \code{c(1, 1)}.
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
#'   either drug, until the first DLT. Ignored when the first cohort size is one.
#'   Defaults to \code{FALSE}.
#'
#' @param bound_mtd
#'   Logical scalar. Require the isotonic estimate at the selected combination
#'   to be below the de-escalation boundary. Defaults to \code{FALSE}.
#'
#' @param overdose_cutoff
#'   Numeric scalar or \code{NULL}. Combinations whose true DLT probability
#'   exceeds this value count as overdoses in the \code{overdose} component of
#'   the result. Defaults to \code{NULL}, which uses \code{target}.
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
#'   An object of class \code{comb_boin_oc}, which is a list with components
#'   \item{sel_percent}{Matrix of the percentage of trials selecting each combination as the MTD.}
#'   \item{percent_no_mtd}{Percentage of trials ending without an MTD.}
#'   \item{n_pts_dose}{Matrix of the average number of patients treated at each combination.}
#'   \item{n_tox_dose}{Matrix of the average number of DLTs observed at each combination.}
#'   \item{total_n_pts}{Average total number of patients per trial.}
#'   \item{total_n_tox}{Average total number of DLTs per trial.}
#'   \item{true_mtd}{Two-column matrix of the dose levels of drug A and drug B of the combinations whose true DLT probability is closest to \code{target}.}
#'   \item{pct_true_mtd_selected}{Percentage of trials selecting one of these combinations.}
#'   \item{pct_patients_true_mtd}{Percentage of all patients treated at one of them.}
#'   \item{overdose}{List describing the exposure to combinations above \code{overdose_cutoff}, with the components of the same list returned by \code{\link{sim_boin}}; \code{combinations} gives the dose levels of the combinations above the cutoff.}
#'   \item{stop_reason_percent}{Percentage of trials by reason for stopping.}
#'   \item{trials}{Trial level data when \code{keep_trials} is \code{TRUE}, otherwise \code{NULL}.}
#'   together with the design parameters and the call.
#'
#' @details
#'   The trials are simulated with \code{\link{comb_boin_simulate}} and the MTD
#'   is selected with \code{\link{comb_boin_select_mtd}}. A trial stopped for
#'   safety at the lowest combination selects none. Both follow
#'   \code{BOIN::get.oc.comb()} with \code{mtd.contour = FALSE}: with the same
#'   seed and matching arguments the two implementations agree trial by trial.
#'   \code{pct_true_mtd_selected} and \code{pct_patients_true_mtd} are the
#'   quantities that the reference implementation reports as the selection
#'   percentage of the MTD and the percentage of patients treated at the MTD.
#'
#'   Note that \code{n_earlystop} defaults to 18 here, whereas the reference
#'   implementation defaults to 100, which in practice switches the rule off.
#'
#'   Only the design that searches for a single MTD combination is provided; the
#'   waterfall design of Zhang and Yuan (2016), which searches for the MTD
#'   contour, is not.
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
#'   Zhang, L. and Yuan, Y. (2016). A Practical Bayesian Design to Identify the
#'   Maximum Tolerated Dose Contour for Drug Combination Trials. Statistics in
#'   Medicine, 35(27), 4924-4936.
#'
#' @examples
#' # The example of Yan and colleagues (2020): three levels of drug A, four of
#' # drug B and a target DLT rate of 0.25
#' p_true <- matrix(c(0.02, 0.04, 0.08, 0.14,
#'                    0.08, 0.25, 0.42, 0.48,
#'                    0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE)
#'
#' oc <- sim_comb_boin(
#'   target = 0.25,
#'   p_true = p_true,
#'   n_cohort = 16,
#'   cohort_size = 3,
#'   n_trials = 500,
#'   seed = 123
#' )
#' oc
#'
#' \donttest{
#' # The settings of the article, which turn early stopping off
#' oc_jss <- sim_comb_boin(
#'   target = 0.25,
#'   p_true = p_true,
#'   n_cohort = 16,
#'   cohort_size = 3,
#'   n_earlystop = 100,
#'   n_trials = 1000,
#'   seed = 6
#' )
#' oc_jss$pct_true_mtd_selected
#' }
#'
#' @seealso \code{\link{sim_comb_boin_multi}}, \code{\link{comb_boin_simulate}},
#'   \code{\link{comb_boin_select_mtd}}
#'
#' @export
sim_comb_boin <- function(target, p_true, n_cohort, cohort_size,
                          n_trials = 10000, start_dose = c(1, 1),
                          n_earlystop = 18, p_saf = NULL, p_tox = NULL,
                          cutoff_eli = 0.95, extrasafe = FALSE, offset = 0.05,
                          titration = FALSE, bound_mtd = FALSE,
                          overdose_cutoff = NULL, keep_trials = FALSE,
                          verbose = FALSE, seed = 123) {

  call_expr <- match.call()

  # A cutoff given explicitly is checked before the simulation runs; the
  # default is the target, which comb_boin_simulate() validates.
  if (!is.null(overdose_cutoff)) {
    check_scalar_prob(overdose_cutoff, "overdose_cutoff")
  }
  check_flag(bound_mtd, "bound_mtd")

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- comb_boin_simulate(
    target = target, p_true = p_true, n_cohort = n_cohort,
    cohort_size = cohort_size, n_trials = n_trials, start_dose = start_dose,
    n_earlystop = n_earlystop, p_saf = p_saf, p_tox = p_tox,
    cutoff_eli = cutoff_eli, extrasafe = extrasafe, offset = offset,
    titration = titration, seed = seed
  )

  if (verbose) message("Selecting the MTD ...")

  # A trial stopped for safety selects no combination, and the reference
  # implementation does not run the selection for it.
  safety_stop <- trials$stop_reason %in%
    c("lowest_dose_eliminated", "lowest_dose_too_toxic")
  mtd_a <- rep(NA_integer_, n_trials)
  mtd_b <- rep(NA_integer_, n_trials)
  selection_reason <- trials$stop_reason
  keep <- which(!safety_stop)
  if (length(keep) > 0L) {
    selection <- comb_boin_select_mtd(
      n_pts = trials$n_pts[keep, , , drop = FALSE],
      n_tox = trials$n_tox[keep, , , drop = FALSE],
      target = target, cutoff_eli = cutoff_eli, extrasafe = extrasafe,
      offset = offset, bound_mtd = bound_mtd, p_tox = trials$settings$p_tox
    )
    mtd_a[keep] <- selection$mtd_a
    mtd_b[keep] <- selection$mtd_b
    selection_reason[keep] <- selection$reason
  }

  n_rows <- nrow(p_true)
  n_cols <- ncol(p_true)
  labels <- list(paste0("A", seq_len(n_rows)), paste0("B", seq_len(n_cols)))

  selected <- !is.na(mtd_a)
  cell <- mtd_a[selected] + n_rows * (mtd_b[selected] - 1L)
  sel_count <- tabulate(cell, nbins = n_rows * n_cols)
  sel_percent <- matrix(sel_count / n_trials * 100, n_rows, n_cols,
                        dimnames = labels)
  percent_no_mtd <- mean(!selected) * 100

  # Averages over trials are taken as the reference implementation takes them,
  # from double precision arrays with mean(), so that they agree to the bit.
  n_pts_num <- array(as.numeric(trials$n_pts), dim(trials$n_pts))
  n_tox_num <- array(as.numeric(trials$n_tox), dim(trials$n_tox))
  n_pts_dose <- matrix(apply(n_pts_num, c(2L, 3L), mean), n_rows, n_cols,
                       dimnames = labels)
  n_tox_dose <- matrix(apply(n_tox_num, c(2L, 3L), mean), n_rows, n_cols,
                       dimnames = labels)

  closest <- abs(p_true - target) == min(abs(p_true - target))
  true_mtd <- which(closest, arr.ind = TRUE)
  dimnames(true_mtd) <- list(NULL, c("dose_a", "dose_b"))

  if (is.null(overdose_cutoff)) overdose_cutoff <- target
  above <- p_true > overdose_cutoff
  n_by_cell <- matrix(trials$n_pts, nrow = n_trials)
  n_per_trial <- rowSums(n_by_cell)
  max_pts <- trials$settings$max_total_pts
  n_above <- if (any(above)) {
    rowSums(n_by_cell[, which(above), drop = FALSE])
  } else {
    rep(0, n_trials)
  }

  # Whether the combination a trial recommended is itself above the cutoff.
  mtd_above <- rep(FALSE, n_trials)
  mtd_above[selected] <- above[cbind(mtd_a[selected], mtd_b[selected])]

  above_cells <- which(above, arr.ind = TRUE)
  dimnames(above_cells) <- list(NULL, c("dose_a", "dose_b"))

  overdose <- list(
    cutoff = overdose_cutoff,
    combinations = above_cells,
    pct_patients = sum(n_above) / sum(n_per_trial) * 100,
    pct_patients_by_trial = mean(n_above / n_per_trial) * 100,
    avg_n_patients = mean(n_above),
    pct_trials_any = mean(n_above > 0) * 100,
    pct_trials_over_60 = mean(n_above > 0.6 * max_pts) * 100,
    pct_trials_over_80 = mean(n_above > 0.8 * max_pts) * 100,
    pct_trials_mtd_above = mean(mtd_above) * 100,
    pct_mtd_above_when_selected = if (any(selected)) {
      mean(mtd_above[selected]) * 100
    } else {
      NA_real_
    }
  )

  stop_reason_percent <- 100 * table(trials$stop_reason) / n_trials

  if (verbose) message("Done.")

  structure(
    list(
      sel_percent = sel_percent,
      percent_no_mtd = percent_no_mtd,
      n_pts_dose = n_pts_dose,
      n_tox_dose = n_tox_dose,
      total_n_pts = sum(trials$n_pts) / n_trials,
      total_n_tox = sum(trials$n_tox) / n_trials,
      true_mtd = true_mtd,
      pct_true_mtd_selected = sum(sel_percent[closest]),
      pct_patients_true_mtd = sum(n_pts_dose[closest]) / sum(n_pts_dose) * 100,
      overdose = overdose,
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
          cohorts_used = trials$cohorts_used,
          stop_reason = trials$stop_reason,
          mtd_a = mtd_a,
          mtd_b = mtd_b,
          selection_reason = selection_reason
        )
      } else {
        NULL
      },
      call = call_expr
    ),
    class = "comb_boin_oc"
  )
}
