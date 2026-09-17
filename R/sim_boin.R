#' Operating Characteristics of a BOIN Design
#'
#' @description
#'   Simulate a BOIN dose-finding trial many times under one dose-toxicity
#'   scenario and summarise how often each dose is selected as the MTD, how many
#'   patients are treated at each dose and how many DLTs are observed.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param p_true
#'   Numeric vector. True DLT probability at each dose level, in increasing dose
#'   order.
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
#'   Integer scalar. Dose level for the first cohort. Defaults to 1.
#'
#' @param n_earlystop
#'   Integer scalar. The trial stops once this many patients have been treated at
#'   the current dose and the design would stay there. Defaults to 18. Set it to a
#'   value above the maximum sample size to switch this rule off.
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
#'   Logical scalar. Apply the stricter safety stopping rule at the lowest dose.
#'   Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is
#'   relaxed for the safety stopping rule. Defaults to 0.05.
#'
#' @param titration
#'   Logical scalar. Start with single patient cohorts until the first DLT is
#'   seen. Ignored when the first cohort size is one. Defaults to \code{FALSE}.
#'
#' @param bound_mtd
#'   Logical scalar. Require the isotonic estimate at the selected dose to be at
#'   or below the de-escalation boundary. Defaults to \code{FALSE}.
#'
#' @param min_mtd_sample
#'   Integer scalar. Smallest number of patients a dose must have received to be
#'   eligible as the MTD. Defaults to 1.
#'
#' @param n_earlystop_rule
#'   Character scalar, either \code{"with_stay"} or \code{"simple"}. See
#'   \code{\link{boin_simulate}}.
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
#'   An object of class \code{boin_oc}, which is a list with components
#'   \item{sel_percent}{Percentage of trials selecting each dose as the MTD.}
#'   \item{percent_no_mtd}{Percentage of trials ending without an MTD.}
#'   \item{n_pts_dose}{Average number of patients treated at each dose.}
#'   \item{n_tox_dose}{Average number of DLTs observed at each dose.}
#'   \item{total_n_pts}{Average total number of patients per trial.}
#'   \item{total_n_tox}{Average total number of DLTs per trial.}
#'   \item{overdose60}{Percentage of trials treating more than 60 percent of patients above the target, \code{NA} when no dose is above the target.}
#'   \item{overdose80}{The same at 80 percent.}
#'   \item{stop_reason_percent}{Percentage of trials by reason for stopping.}
#'   \item{trials}{Trial level data when \code{keep_trials} is \code{TRUE}, otherwise \code{NULL}.}
#'   together with the design parameters and the call.
#'
#' @details
#'   The engine consumes one uniform random variate per patient, in enrollment
#'   order, and applies the decision rules in the order used by
#'   \code{BOIN::get.oc()}. With the same seed and matching arguments the two
#'   implementations therefore agree trial by trial, not merely on average.
#'
#'   Note that \code{n_earlystop} defaults to 18 here, whereas the reference
#'   implementation defaults to 100, which in practice switches the rule off.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' oc <- sim_boin(
#'   target = 0.30,
#'   p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_trials = 500,
#'   seed = 123
#' )
#' oc
#'
#' \donttest{
#' # A larger run with the safety options switched on
#' oc_safe <- sim_boin(
#'   target = 0.30,
#'   p_true = c(0.30, 0.40, 0.50, 0.60, 0.70),
#'   n_cohort = 20,
#'   cohort_size = 3,
#'   n_trials = 10000,
#'   extrasafe = TRUE,
#'   bound_mtd = TRUE,
#'   titration = TRUE,
#'   seed = 123
#' )
#' oc_safe
#' }
#'
#' @seealso \code{\link{sim_boin_multi}}, \code{\link{boin_simulate}}
#'
#' @export
sim_boin <- function(target, p_true, n_cohort, cohort_size,
                     n_trials = 10000, start_dose = 1, n_earlystop = 18,
                     p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                     extrasafe = FALSE, offset = 0.05, titration = FALSE,
                     bound_mtd = FALSE, min_mtd_sample = 1,
                     n_earlystop_rule = c("with_stay", "simple"),
                     keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  n_earlystop_rule <- match.arg(n_earlystop_rule)

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- boin_simulate(
    target = target, p_true = p_true, n_cohort = n_cohort,
    cohort_size = cohort_size, n_trials = n_trials, start_dose = start_dose,
    n_earlystop = n_earlystop, p_saf = p_saf, p_tox = p_tox,
    cutoff_eli = cutoff_eli, extrasafe = extrasafe, offset = offset,
    titration = titration, n_earlystop_rule = n_earlystop_rule, seed = seed
  )

  if (verbose) message("Selecting the MTD ...")

  selection <- boin_select_mtd(
    n_pts = trials$n_pts, n_tox = trials$n_tox, target = target,
    cutoff_eli = cutoff_eli, extrasafe = extrasafe, offset = offset,
    bound_mtd = bound_mtd, p_tox = trials$settings$p_tox,
    min_mtd_sample = min_mtd_sample
  )

  # A trial stopped for safety selects no dose, whatever the final data show.
  safety_stop <- trials$stop_reason %in%
    c("lowest_dose_eliminated", "lowest_dose_too_toxic")
  selection$mtd[safety_stop] <- NA_integer_
  selection$reason[safety_stop] <- trials$stop_reason[safety_stop]

  n_doses <- length(p_true)
  dose_names <- paste0("DL", seq_len(n_doses))

  selected <- factor(selection$mtd, levels = seq_len(n_doses))
  sel_percent <- as.numeric(table(selected)) / n_trials * 100
  names(sel_percent) <- dose_names
  percent_no_mtd <- mean(is.na(selection$mtd)) * 100

  n_pts_dose <- colMeans(trials$n_pts)
  n_tox_dose <- colMeans(trials$n_tox)
  names(n_pts_dose) <- dose_names
  names(n_tox_dose) <- dose_names

  above_target <- p_true > target
  if (any(above_target)) {
    n_above <- rowSums(trials$n_pts[, above_target, drop = FALSE])
    max_pts <- trials$settings$max_total_pts
    overdose60 <- mean(n_above > 0.6 * max_pts) * 100
    overdose80 <- mean(n_above > 0.8 * max_pts) * 100
  } else {
    overdose60 <- NA_real_
    overdose80 <- NA_real_
  }

  stop_reason_percent <- 100 * table(trials$stop_reason) / n_trials

  if (verbose) message("Done.")

  structure(
    list(
      sel_percent = sel_percent,
      percent_no_mtd = percent_no_mtd,
      n_pts_dose = n_pts_dose,
      n_tox_dose = n_tox_dose,
      total_n_pts = mean(rowSums(trials$n_pts)),
      total_n_tox = mean(rowSums(trials$n_tox)),
      overdose60 = overdose60,
      overdose80 = overdose80,
      stop_reason_percent = stop_reason_percent,
      target = target,
      p_true = p_true,
      p_saf = trials$settings$p_saf,
      p_tox = trials$settings$p_tox,
      lambda_e = trials$boundary$lambda_e,
      lambda_d = trials$boundary$lambda_d,
      n_doses = n_doses,
      n_trials = n_trials,
      settings = trials$settings,
      trials = if (keep_trials) {
        list(
          n_pts = trials$n_pts,
          n_tox = trials$n_tox,
          eliminated = trials$eliminated,
          cohorts_used = trials$cohorts_used,
          stop_reason = trials$stop_reason,
          mtd = selection$mtd,
          selection_reason = selection$reason
        )
      } else {
        NULL
      },
      call = call_expr
    ),
    class = "boin_oc"
  )
}
