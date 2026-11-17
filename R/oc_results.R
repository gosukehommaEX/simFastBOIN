#' Operating Characteristics of Several Designs in One Table
#'
#' @description
#'   Collect the results of the simulation functions of different designs
#'   (and of \code{\link{oc_3p3}}) into one long data frame, so that the
#'   designs can be tabulated and drawn side by side. The result has
#'   \code{print()} and \code{plot()} methods; \code{\link{oc_compare}}
#'   judges the same results against the true doses.
#'
#' @param ...
#'   Design results, named by the design, for example
#'   \code{oc_results(BOIN12 = res_a, "BOIN-ET" = res_b)}. Each is the result
#'   of a simulation function for one scenario (such as \code{sim_boin()} or
#'   \code{sim_boin12()}), of its \code{_multi} version, of
#'   \code{\link{oc_3p3}}, or a list of such single-scenario results named by
#'   the scenarios (the way to give \code{oc_3p3()} several scenarios).
#'   Unnamed results are named after their design. All results must be for a
#'   single agent or all for a drug combination, with the same number of
#'   doses.
#'
#' @return
#'   A data frame of class \code{c("oc_results", "data.frame")} with the
#'   columns
#'   \item{design}{Factor of the design names, in the order given.}
#'   \item{target_type}{\code{"MTD"} or \code{"OBD"}, what the design
#'     selects.}
#'   \item{scenario}{Factor of the scenario names. A single-scenario result
#'     is named "Scenario 1".}
#'   \item{dose}{Dose level, or \code{dose_a} and \code{dose_b} for a
#'     combination; \code{NA} for the quantities of the whole trial.}
#'   \item{metric}{Factor: per dose, \code{"sel_percent"} (percentage of
#'     trials selecting the dose), \code{"n_pts"}, \code{"n_tox"} and
#'     \code{"n_resp"} (average numbers of patients, DLTs and responses);
#'     per trial, \code{"percent_no_selection"} (percentage of trials
#'     selecting no dose), \code{"total_n_pts"}, \code{"total_n_tox"},
#'     \code{"total_n_resp"} and \code{"duration"} (average trial duration,
#'     for the designs with time). Quantities a design does not record are
#'     left out.}
#'   \item{value}{The value.}
#'   with the attributes \code{type} (\code{"single"} or
#'   \code{"combination"}) and \code{dims}.
#'
#' @details
#'   The values are those of the components \code{sel_percent},
#'   \code{n_pts_dose}, \code{n_tox_dose}, \code{n_resp_dose},
#'   \code{percent_no_mtd} or \code{percent_no_obd}, \code{total_n_pts},
#'   \code{total_n_tox}, \code{total_n_resp} and \code{duration_mean} or
#'   \code{duration} of the results; nothing is simulated again.
#'
#' @examples
#' scen <- list(
#'   "Scenario 1" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                       p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))
#' )
#' boin12 <- sim_boin12_multi(target = 0.35, target_resp = 0.25,
#'                            scenarios = scen, n_cohort = 12,
#'                            cohort_size = 3, n_trials = 200, seed = 1)
#' boin <- sim_boin_multi(target = 0.35, scenarios = scen, n_cohort = 12,
#'                        cohort_size = 3, n_trials = 200, seed = 1)
#' res <- oc_results(BOIN12 = boin12, BOIN = boin)
#' head(res)
#' res
#'
#' @seealso \code{\link{oc_scenarios}}, \code{\link{oc_compare}},
#'   \code{\link{print.oc_results}}, \code{\link{plot.oc_results}}
#'
#' @export
oc_results <- function(...) {
  sets <- oc_collect(list(...))
  type <- sets[[1L]]$type
  dims <- dim(as.matrix(sets[[1L]]$scenarios[[1L]]$sel))
  if (type == "single") dims <- dims[1L]
  designs <- vapply(sets, function(z) z$design, character(1))
  scenario_names <- unique(unlist(lapply(sets, function(z) {
    names(z$scenarios)
  })))
  dose_metrics <- c("sel_percent", "n_pts", "n_tox", "n_resp")
  trial_metrics <- c("percent_no_selection", "total_n_pts", "total_n_tox",
                     "total_n_resp", "duration")

  rows <- list()
  for (set in sets) {
    for (s in names(set$scenarios)) {
      z <- set$scenarios[[s]]
      per_dose <- list(sel_percent = z$sel, n_pts = z$n_pts,
                       n_tox = z$n_tox, n_resp = z$n_resp)
      for (m in dose_metrics) {
        v <- per_dose[[m]]
        if (is.null(v)) next
        if (type == "single") {
          pos <- data.frame(dose = seq_along(v))
        } else {
          pos <- data.frame(dose_a = as.vector(row(as.matrix(v))),
                            dose_b = as.vector(col(as.matrix(v))))
        }
        rows[[length(rows) + 1L]] <- data.frame(
          design = set$design, target_type = set$target_type, scenario = s,
          pos, metric = m, value = as.vector(unname(v)),
          stringsAsFactors = FALSE
        )
      }
      per_trial <- c(percent_no_selection = z$no_sel,
                     total_n_pts = z$total_n_pts,
                     total_n_tox = z$total_n_tox,
                     total_n_resp = z$total_n_resp,
                     duration = z$duration)
      per_trial <- per_trial[!is.na(per_trial)]
      if (type == "single") {
        pos <- data.frame(dose = rep(NA_integer_, length(per_trial)))
      } else {
        pos <- data.frame(dose_a = rep(NA_integer_, length(per_trial)),
                          dose_b = rep(NA_integer_, length(per_trial)))
      }
      rows[[length(rows) + 1L]] <- data.frame(
        design = rep(set$design, length(per_trial)),
        target_type = rep(set$target_type, length(per_trial)),
        scenario = rep(s, length(per_trial)), pos,
        metric = names(per_trial), value = unname(per_trial),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$design <- factor(out$design, levels = designs)
  out$scenario <- factor(out$scenario, levels = scenario_names)
  out$metric <- factor(out$metric, levels = c(dose_metrics, trial_metrics))
  out$metric <- droplevels(out$metric)
  rownames(out) <- NULL
  structure(out, class = c("oc_results", "data.frame"), type = type,
            dims = dims)
}
