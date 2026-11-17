#' Compare Several Designs Against the True Doses
#'
#' @description
#'   Judge the results of several designs against the true doses of the
#'   scenarios: how often each design selects the true MTD (designs that
#'   select an MTD) or the true OBD (designs that select an OBD), how many
#'   patients it treats there, how often it selects an overdose and how many
#'   patients it treats at overdoses, how often it selects no dose, and its
#'   average sample size and duration. The result has \code{print()} and
#'   \code{plot()} methods.
#'
#' @param ...
#'   Design results, named by the design, as for \code{\link{oc_results}}.
#'
#' @param truth
#'   The result of \code{\link{oc_scenarios}} for the same scenarios, with the
#'   same scenario names. The DLT probabilities of every scenario must agree
#'   with those stored in the results. Designs that select an OBD need
#'   scenarios with response probabilities.
#'
#' @param overdose_cutoff
#'   Numeric scalar or \code{NULL}. A dose whose true DLT probability is above
#'   this value is an overdose. Defaults to \code{NULL}, the \code{tox_limit}
#'   of \code{truth}.
#'
#' @return
#'   A data frame of class \code{c("oc_compare", "data.frame")} with one row
#'   per design and scenario and the columns
#'   \item{design, target_type, scenario}{As in \code{\link{oc_results}}.}
#'   \item{pct_correct}{Percentage of trials selecting a true MTD or OBD;
#'     \code{NA} when the scenario has none.}
#'   \item{n_pts_correct, pct_patients_correct}{Average number of patients
#'     treated at the true MTDs or OBDs, and its percentage of all the
#'     patients.}
#'   \item{pct_overdose_selected}{Percentage of trials selecting an
#'     overdose.}
#'   \item{n_pts_overdose, pct_patients_overdose}{Average number of patients
#'     treated at overdoses, and its percentage of all the patients.}
#'   \item{pct_no_selection}{Percentage of trials selecting no dose.}
#'   \item{total_n_pts}{Average number of patients per trial.}
#'   \item{duration}{Average trial duration, \code{NA} for the designs
#'     without time.}
#'   with the attribute \code{overdose_cutoff}.
#'
#' @details
#'   Every quantity is computed from the per-dose averages of the results
#'   (\code{sel_percent} and \code{n_pts_dose}), so nothing is simulated
#'   again and the results need not keep the trials. A dose is an overdose
#'   when its DLT probability exceeds \code{overdose_cutoff} by more than
#'   \code{1e-12}.
#'
#' @examples
#' scen <- list(
#'   "Scenario 1" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                       p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
#'   "Scenario 2" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                       p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
#' )
#' truth <- oc_scenarios(scen, target = 0.35, target_resp = 0.25)
#' boin12 <- sim_boin12_multi(target = 0.35, target_resp = 0.25,
#'                            scenarios = scen, n_cohort = 12,
#'                            cohort_size = 3, n_trials = 200, seed = 1)
#' boin <- sim_boin_multi(target = 0.35, scenarios = scen, n_cohort = 12,
#'                        cohort_size = 3, n_trials = 200, seed = 1)
#' oc_compare(BOIN12 = boin12, BOIN = boin, truth = truth)
#'
#' @seealso \code{\link{oc_scenarios}}, \code{\link{oc_results}},
#'   \code{\link{print.oc_compare}}, \code{\link{plot.oc_compare}}
#'
#' @export
oc_compare <- function(..., truth, overdose_cutoff = NULL) {
  if (missing(truth) || !inherits(truth, "oc_scenarios") ||
      is.null(attr(truth, "type"))) {
    stop("'truth' must be the result of oc_scenarios()", call. = FALSE)
  }
  if (is.null(overdose_cutoff)) overdose_cutoff <- attr(truth, "tox_limit")
  check_scalar_prob(overdose_cutoff, "overdose_cutoff")
  sets <- oc_collect(list(...))
  type <- sets[[1L]]$type
  if (type != attr(truth, "type")) {
    stop("the designs and 'truth' must be both for a single agent or both ",
         "for a drug combination", call. = FALSE)
  }
  n_cells <- length(sets[[1L]]$scenarios[[1L]]$sel)
  if (n_cells != prod(attr(truth, "dims"))) {
    stop("the designs and 'truth' must have the same number of doses",
         call. = FALSE)
  }
  has_resp <- "true_obd" %in% names(truth)

  rows <- list()
  for (set in sets) {
    if (set$target_type == "OBD" && !has_resp) {
      stop("design '", set$design, "' selects an OBD, so 'truth' needs ",
           "response probabilities", call. = FALSE)
    }
    for (s in names(set$scenarios)) {
      if (!s %in% levels(truth$scenario)) {
        stop("scenario '", s, "' of design '", set$design, "' is not in ",
             "'truth'", call. = FALSE)
      }
      z <- set$scenarios[[s]]
      tr <- truth[truth$scenario == s, , drop = FALSE]
      p_true <- as.vector(z$p_true)
      if (length(p_true) != nrow(tr) ||
          any(abs(p_true - tr$p_true) > 1e-10)) {
        stop("the DLT probabilities of scenario '", s, "' differ between ",
             "design '", set$design, "' and 'truth'", call. = FALSE)
      }
      correct <- if (set$target_type == "MTD") tr$true_mtd else tr$true_obd
      over <- tr$p_true - overdose_cutoff > 1e-12
      sel <- as.vector(z$sel)
      n_pts <- as.vector(z$n_pts)
      n_all <- sum(n_pts)
      rows[[length(rows) + 1L]] <- data.frame(
        design = set$design,
        target_type = set$target_type,
        scenario = s,
        pct_correct = if (any(correct)) sum(sel[correct]) else NA_real_,
        n_pts_correct = sum(n_pts[correct]),
        pct_patients_correct = if (n_all > 0) {
          100 * sum(n_pts[correct]) / n_all
        } else {
          NA_real_
        },
        pct_overdose_selected = sum(sel[over]),
        n_pts_overdose = sum(n_pts[over]),
        pct_patients_overdose = if (n_all > 0) {
          100 * sum(n_pts[over]) / n_all
        } else {
          NA_real_
        },
        pct_no_selection = z$no_sel,
        total_n_pts = z$total_n_pts,
        duration = z$duration,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$design <- factor(out$design,
                       levels = vapply(sets, function(z) z$design,
                                       character(1)))
  out$scenario <- factor(out$scenario, levels = levels(truth$scenario))
  out$scenario <- droplevels(out$scenario)
  rownames(out) <- NULL
  structure(out, class = c("oc_compare", "data.frame"),
            overdose_cutoff = overdose_cutoff)
}
