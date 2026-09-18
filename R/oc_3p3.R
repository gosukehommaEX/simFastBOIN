#' Exact Operating Characteristics of the 3+3 Design
#'
#' @description
#'   Compute the operating characteristics of the traditional 3+3 dose-finding
#'   design exactly. The probability of every path the design can take is
#'   evaluated in closed form, so the result carries no Monte Carlo error and no
#'   seed is involved. Use it as a comparator for the BOIN results returned by
#'   \code{\link{sim_boin}}.
#'
#' @param p_true
#'   Numeric vector. True DLT probability at each dose level, in increasing dose
#'   order.
#'
#' @param mtd_rule
#'   Character scalar, either \code{"previous"} or \code{"expand"}. The two
#'   definitions of the MTD in common use. See Details. Defaults to
#'   \code{"previous"}.
#'
#' @param start_dose
#'   Integer scalar. Dose level for the first cohort. Defaults to 1.
#'
#' @param overdose_cutoff
#'   Numeric scalar or \code{NULL}. Doses whose true DLT probability exceeds this
#'   value count as overdoses in the \code{overdose} component of the result.
#'   Defaults to \code{NULL}, which uses one third.
#'
#' @return
#'   An object of class \code{oc_3p3}, which is a list with components
#'   \item{sel_percent}{Percentage of trials selecting each dose as the MTD.}
#'   \item{percent_no_mtd}{Percentage of trials ending without an MTD.}
#'   \item{n_pts_dose}{Average number of patients treated at each dose.}
#'   \item{n_tox_dose}{Average number of DLTs observed at each dose.}
#'   \item{total_n_pts}{Average total number of patients per trial.}
#'   \item{total_n_tox}{Average total number of DLTs per trial.}
#'   \item{overdose}{List describing the design's relationship with doses above \code{overdose_cutoff}: \code{pct_patients} (the percentage of all patients treated at such a dose, that is the probability that a patient is dosed above the cutoff), \code{avg_n_patients}, \code{pct_trials_any} (the percentage of trials treating at least one patient there), \code{pct_trials_mtd_above} (the percentage of all trials whose selected MTD is above the cutoff) and \code{pct_mtd_above_when_selected} (the same among the trials that selected an MTD). The \code{cutoff} and the dose levels \code{doses} that exceed it are also returned.}
#'   together with \code{p_true}, \code{n_doses}, \code{mtd_rule},
#'   \code{start_dose}, \code{method} and the call.
#'
#' @details
#'   The escalation rule is the traditional one. Three patients are treated at
#'   the current dose. No DLT leads to the next dose up. One DLT out of three
#'   leads to three more patients at the same dose, and then to the next dose up
#'   if no further DLT is seen and to stopping otherwise. Two or more DLTs out of
#'   three stop the trial at once. A trial that escalates past the highest dose
#'   ends there.
#'
#'   The two MTD rules differ only in what happens once escalation has stopped.
#'   \code{"previous"} names the dose below the one declared too toxic, and the
#'   highest dose when escalation ran out of doses. \code{"expand"} instead takes
#'   the MTD to be the highest dose at which at most one of six patients had a
#'   DLT: starting from the highest dose still eligible, a dose that has only
#'   three patients is expanded to six and assessed, and the search moves down a
#'   level whenever the expanded dose fails. Neither rule is universally the
#'   standard, so the comparator being matched should be checked before a rule is
#'   chosen. The two agree on which doses are visited during escalation, and so
#'   on the probability that a patient is exposed to an overly toxic dose; they
#'   differ in the sample size and in the selection percentages.
#'
#'   The 3+3 design has no fixed maximum sample size, so the within-trial
#'   exposure percentage and the 60 percent and 80 percent thresholds reported
#'   for BOIN in \code{\link{sim_boin}} have no counterpart here.
#'
#' @examples
#' oc_3p3(p_true = c(0.30, 0.48, 0.67))
#'
#' # The other definition of the MTD
#' oc_3p3(p_true = c(0.30, 0.48, 0.67), mtd_rule = "expand")
#'
#' # Exposure to doses above a stated DLT rate
#' oc_3p3(p_true = c(0.20, 0.40, 0.60), overdose_cutoff = 1 / 3)$overdose
#'
#' @seealso \code{\link{sim_3p3}}, \code{\link{sim_boin}}
#'
#' @export
oc_3p3 <- function(p_true, mtd_rule = c("previous", "expand"), start_dose = 1,
                   overdose_cutoff = NULL) {

  call_expr <- match.call()
  check_p_true(p_true)
  mtd_rule <- match.arg(mtd_rule)
  check_count(start_dose, "start_dose")
  start_dose <- as.integer(start_dose)
  n_doses <- length(p_true)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of dose levels", call. = FALSE)
  }
  if (is.null(overdose_cutoff)) overdose_cutoff <- 1 / 3
  check_scalar_prob(overdose_cutoff, "overdose_cutoff")

  q <- p_true
  no_dlt_3 <- (1 - q)^3                 # no DLT among three patients
  one_dlt_3 <- 3 * q * (1 - q)^2        # exactly one DLT among three patients
  one_of_six <- one_dlt_3 * no_dlt_3    # one DLT, then three more without one
  escalate <- no_dlt_3 + one_of_six

  levels_used <- start_dose:n_doses

  reach <- numeric(n_doses)
  running <- 1
  for (d in levels_used) {
    reach[d] <- running
    running <- running * escalate[d]
  }

  stop_toxic <- numeric(n_doses)
  stop_toxic[levels_used] <- reach[levels_used] * (1 - escalate[levels_used])
  exhausted <- reach[n_doses] * escalate[n_doses]

  n_pts_dose <- numeric(n_doses)
  n_tox_dose <- numeric(n_doses)
  n_pts_dose[levels_used] <- reach[levels_used] * (3 + 3 * one_dlt_3[levels_used])
  n_tox_dose[levels_used] <- reach[levels_used] * 3 * q[levels_used] *
    (1 + one_dlt_3[levels_used])

  sel <- numeric(n_doses)

  if (mtd_rule == "previous") {

    if (start_dose < n_doses) {
      sel[start_dose:(n_doses - 1L)] <- stop_toxic[(start_dose + 1L):n_doses]
    }
    sel[n_doses] <- exhausted
    no_mtd <- stop_toxic[start_dose]

  } else {

    # A dose that was escalated from had either no DLT out of three or one out of
    # six. Only the first kind needs three more patients before it can be
    # assessed, and it passes if at most one DLT is seen in six.
    only_three <- ifelse(escalate > 0, no_dlt_3 / escalate, 0)
    already_six <- ifelse(escalate > 0, one_of_six / escalate, 0)
    pass_expansion <- no_dlt_3 + one_dlt_3
    accept <- only_three * pass_expansion + already_six
    reject <- only_three * (1 - pass_expansion)

    examined <- numeric(n_doses)
    for (k in rev(levels_used)) {
      entry <- if (k == n_doses) exhausted else stop_toxic[k + 1L]
      carried <- if (k < n_doses) examined[k + 1L] * reject[k + 1L] else 0
      examined[k] <- entry + carried
    }

    sel[levels_used] <- examined[levels_used] * accept[levels_used]
    expanded <- examined[levels_used] * only_three[levels_used] * 3
    n_pts_dose[levels_used] <- n_pts_dose[levels_used] + expanded
    n_tox_dose[levels_used] <- n_tox_dose[levels_used] + expanded * q[levels_used]
    no_mtd <- stop_toxic[start_dose] + examined[start_dose] * reject[start_dose]
  }

  above <- p_true > overdose_cutoff
  total_n <- sum(n_pts_dose)
  n_above <- sum(n_pts_dose[above])
  lowest_above <- if (any(above)) min(which(above)) else NA_integer_
  pct_trials_any <- if (!is.na(lowest_above) && lowest_above >= start_dose) {
    reach[lowest_above] * 100
  } else {
    0
  }
  pct_trials_mtd_above <- sum(sel[above]) * 100

  overdose <- list(
    cutoff = overdose_cutoff,
    doses = which(above),
    pct_patients = if (total_n > 0) n_above / total_n * 100 else 0,
    avg_n_patients = n_above,
    pct_trials_any = pct_trials_any,
    pct_trials_mtd_above = pct_trials_mtd_above,
    pct_mtd_above_when_selected = if (no_mtd < 1) {
      pct_trials_mtd_above / (1 - no_mtd)
    } else {
      NA_real_
    }
  )

  structure(
    list(
      sel_percent = sel * 100,
      percent_no_mtd = no_mtd * 100,
      n_pts_dose = n_pts_dose,
      n_tox_dose = n_tox_dose,
      total_n_pts = total_n,
      total_n_tox = sum(n_tox_dose),
      overdose = overdose,
      p_true = p_true,
      n_doses = n_doses,
      mtd_rule = mtd_rule,
      start_dose = start_dose,
      method = "exact",
      n_trials = NA_integer_,
      call = call_expr
    ),
    class = "oc_3p3"
  )
}
