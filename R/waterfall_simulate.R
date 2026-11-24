#' Simulate Waterfall Drug-Combination Trials
#'
#' @description
#'   Simulate trials of two drugs given together, run with the waterfall design
#'   of Zhang and Yuan (2016), which looks for the MTD contour, that is, for one
#'   MTD at every dose level of drug A. This is the simulation engine behind
#'   \code{\link{sim_waterfall}}, exposed for users who want the raw trial
#'   data.
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
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. The state of the calling session
#'   is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{waterfall_trials}, which is a list with components
#'   \item{n_pts}{Integer array with dimensions trial by level of drug A by level of drug B: number of patients treated at each combination.}
#'   \item{n_tox}{Integer array of the same dimensions: number of DLTs observed.}
#'   \item{eliminated}{Logical array of the same dimensions: whether the combination is excluded from the selection of the MTD contour.}
#'   \item{n_subtrials}{Integer vector: number of subtrials run in each trial.}
#'   \item{n_enrolled}{Integer vector: number of patients treated in each trial, including those of a subtrial whose data were dropped.}
#'   \item{stop_reason}{Character vector: why each trial ended.}
#'   \item{boundary}{The \code{boin_boundary} object used for the decisions.}
#'   \item{settings}{List of the design parameters.}
#'
#' @details
#'   The dose matrix is divided into one-dimensional subtrials, each run with
#'   the BOIN design for a single agent. The first subtrial starts at the lowest
#'   combination, moves up the first column (the lowest level of drug B) to the
#'   highest level of drug A and then along that row. When a subtrial ends, the
#'   dose whose isotonic estimate of the DLT probability is closest to
#'   \code{target} is its candidate MTD. The next subtrial runs along the row
#'   below the row of the candidate, over the second to the last level of drug
#'   B, and starts one level of drug B above the candidate (at the last level
#'   when the candidate is there). The trial ends after the subtrial of the
#'   lowest level of drug A.
#'
#'   Within a subtrial, the dose escalates, stays or de-escalates by the
#'   boundaries of \code{\link{boin_boundary}}; a dose with too many DLTs is
#'   eliminated together with the doses above it, and the subtrial ends without
#'   a candidate when its lowest dose is eliminated. Otherwise a subtrial ends
#'   when \code{n_earlystop} patients have been treated at the dose chosen for
#'   the next cohort, or when its maximum number of cohorts, \code{n_cohort},
#'   is reached. The default of \code{n_cohort} is
#'   \code{ceiling(round(4 / cohort_size, 2) * m)} for a subtrial of \code{m}
#'   doses, the default of the reference implementation, which follows the
#'   article's suggestion of four patients per dose.
#'
#'   The rules are those of \code{BOIN::get.oc.comb()} with
#'   \code{mtd.contour = TRUE} in version 2.7.2 of the 'BOIN' package, the
#'   implementation by the authors of the article. Random variates are drawn
#'   one per patient in enrollment order, as that function draws them, so with
#'   \code{boin_compatible = TRUE}, the same seed and matching arguments the
#'   two give the same trials. Three rules of that function are not in the
#'   article:
#'   \itemize{
#'     \item When the candidate MTD of the first subtrial lies in the first
#'       column below the highest level of drug A, the rows above it are
#'       excluded, and when the data at the candidate would lead to an
#'       escalation, an additional subtrial explores the rest of its row.
#'     \item The combinations to the right of a candidate MTD in its row are
#'       excluded from the selection of the MTD contour.
#'     \item The trial ends when a subtrial ends without a candidate MTD; the
#'       rows below are then not explored.
#'   }
#'   With \code{extrasafe = TRUE}, that function evaluates the safety rule with
#'   a \code{Beta(y + 0.05, n - y + 0.1)} distribution during a subtrial and
#'   with a \code{Beta(y + 1, n - y + 1)} distribution when it chooses the
#'   candidate; the package does the same.
#'
#'   The package differs from that function in two points. When a subtrial
#'   other than the additional one ends without a candidate, that function
#'   drops its data, so that its patients and DLTs are missing from the counts
#'   and a trial stopped in the first subtrial has no patients. By default the
#'   package keeps these data and marks every dose of that subtrial as
#'   eliminated; \code{boin_compatible = TRUE} drops them. And that function
#'   takes the escalation and de-escalation boundaries from the default values
#'   of \code{p_saf} and \code{p_tox} whatever values are given; the package
#'   uses the values given.
#'
#'   Values of \code{stop_reason} are \code{"completed"},
#'   \code{"first_subtrial_no_mtd"}, \code{"later_subtrial_no_mtd"} and
#'   \code{"max_sample_size"}.
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
#' p_true <- matrix(c(0.03, 0.10, 0.28,
#'                    0.10, 0.30, 0.50), nrow = 2, byrow = TRUE)
#'
#' trials <- waterfall_simulate(
#'   target = 0.30,
#'   p_true = p_true,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#' trials
#'
#' @seealso \code{\link{sim_waterfall}}, \code{\link{waterfall_select_mtd}}
#'
#' @export
waterfall_simulate <- function(target, p_true, cohort_size, n_cohort = NULL,
                               n_trials = 10000, n_earlystop = 12,
                               p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                               extrasafe = FALSE, offset = 0.05,
                               titration = FALSE, bound_mtd = FALSE,
                               boin_compatible = FALSE, seed = 123) {

  check_p_true_matrix(p_true)
  n_rows <- nrow(p_true)
  n_cols <- ncol(p_true)
  if (n_rows < 2L) {
    stop("'p_true' must have at least two rows (levels of drug A)",
         call. = FALSE)
  }
  if (n_rows > n_cols) {
    stop("'p_true' must not have more rows than columns; transpose it so ",
         "that the drug with fewer levels is drug A", call. = FALSE)
  }
  check_count(cohort_size, "cohort_size", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(n_earlystop, "n_earlystop", 1L)
  check_flag(titration, "titration")
  check_flag(bound_mtd, "bound_mtd")
  check_flag(boin_compatible, "boin_compatible")
  if (n_earlystop <= 6) {
    warning("'n_earlystop' is low; values between 9 and 18 are recommended",
            call. = FALSE)
  }

  # One subtrial per level of drug A: the first holds the first column and
  # the last row, the others the columns 2 to K of one row.
  n_doses <- c(n_rows + n_cols - 1L, rep(n_cols - 1L, n_rows - 1L))
  if (is.null(n_cohort)) {
    n_cohort <- ceiling(round(4 / cohort_size, 2) * n_doses)
  }
  if (!is.numeric(n_cohort) || length(n_cohort) != n_rows ||
      any(!is.finite(n_cohort)) || any(n_cohort != round(n_cohort)) ||
      any(n_cohort < 1)) {
    stop("'n_cohort' must hold one whole number of at least 1 per level of ",
         "drug A, the maximum number of cohorts of each subtrial",
         call. = FALSE)
  }
  n_cohort <- as.integer(n_cohort)

  # Titration needs a cohort size above one, otherwise it is the design itself.
  if (titration && cohort_size == 1L) titration <- FALSE

  # No dose can receive more patients than the cohorts of its subtrial hold.
  max_dose_pts <- max(n_cohort) * as.integer(cohort_size)
  max_total_pts <- sum(n_cohort) * as.integer(cohort_size)

  bound <- boin_boundary(target, max_dose_pts, p_saf = p_saf, p_tox = p_tox,
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

  res <- waterfall_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    n_rows = as.integer(n_rows),
    n_cols = as.integer(n_cols),
    n_cohort = n_cohort,
    cohort_size = as.integer(cohort_size),
    n_earlystop = as.integer(n_earlystop),
    titration = titration,
    extrasafe = extrasafe,
    bound_mtd = bound_mtd,
    boin_compatible = boin_compatible,
    target = as.numeric(target),
    cutoff_eli = as.numeric(cutoff_eli),
    offset = as.numeric(offset),
    lambda_d = bound$lambda_d,
    b_esc = bound$b_esc,
    b_deesc = bound$b_deesc,
    b_elim = as.integer(b_elim)
  )

  reasons <- c("completed", "first_subtrial_no_mtd", "later_subtrial_no_mtd",
               "max_sample_size")

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
      n_subtrials = res$n_subtrials,
      n_enrolled = res$n_enrolled,
      stop_reason = reasons[res$stop_code + 1L],
      boundary = bound,
      settings = list(
        target = target,
        p_true = p_true,
        n_cohort = n_cohort,
        cohort_size = as.integer(cohort_size),
        n_trials = n_trials,
        n_earlystop = n_earlystop,
        p_saf = bound$p_saf,
        p_tox = bound$p_tox,
        cutoff_eli = cutoff_eli,
        extrasafe = extrasafe,
        offset = offset,
        titration = titration,
        bound_mtd = bound_mtd,
        boin_compatible = boin_compatible,
        max_total_pts = max_total_pts,
        seed = seed
      )
    ),
    class = "waterfall_trials"
  )
}
