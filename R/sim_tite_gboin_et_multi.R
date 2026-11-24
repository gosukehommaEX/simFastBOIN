#' Operating Characteristics of the TITE-gBOIN-ET Design in Several Scenarios
#'
#' @description
#'   Run \code{\link{sim_tite_gboin_et}} for several scenarios of
#'   probabilities of toxicity and efficacy categories with the same design,
#'   and collect the results in one table.
#'
#' @inheritParams sim_tite_gboin_et
#' @inheritParams sim_gboin_et_multi
#'
#' @return
#'   An object of class \code{boin_et_oc_multi}, as returned by
#'   \code{\link{sim_gboin_et_multi}}, whose \code{results} are the
#'   \code{boin_et_oc} objects of \code{\link{sim_tite_gboin_et}}.
#'
#' @details
#'   Every scenario is simulated with the same seed, as in
#'   \code{\link{sim_boin_et_multi}}.
#'
#' @references
#'   Takeda, K., Yamaguchi, Y., Taguri, M. and Morita, S. (2023).
#'   TITE-gBOIN-ET: Time-to-event generalized Bayesian optimal interval
#'   design to accelerate dose-finding accounting for ordinal graded efficacy
#'   and toxicity outcomes. Biometrical Journal, 65(7), 2200265.
#'   \doi{10.1002/bimj.202200265}
#'
#' @examples
#' # The same toxicity grades with an increasing and with a plateauing
#' # efficacy
#' p_grade <- rbind(
#'   c(0.94, 0.79, 0.62, 0.50),
#'   c(0.05, 0.15, 0.20, 0.20),
#'   c(0.01, 0.05, 0.15, 0.25),
#'   c(0.00, 0.01, 0.03, 0.05)
#' )
#' scenarios <- list(
#'   "Increasing" = list(
#'     p_grade = p_grade,
#'     p_resp_grade = rbind(
#'       c(0.64, 0.45, 0.20, 0.05),
#'       c(0.30, 0.40, 0.40, 0.15),
#'       c(0.05, 0.10, 0.20, 0.35),
#'       c(0.01, 0.05, 0.20, 0.45)
#'     )
#'   ),
#'   "Plateau" = list(
#'     p_grade = p_grade,
#'     p_resp_grade = rbind(
#'       c(0.44, 0.10, 0.05, 0.05),
#'       c(0.40, 0.30, 0.20, 0.20),
#'       c(0.15, 0.35, 0.25, 0.25),
#'       c(0.01, 0.25, 0.50, 0.50)
#'     )
#'   )
#' )
#' res <- sim_tite_gboin_et_multi(
#'   target = 0.47 / 1.5, target_resp = 1.75 / 3, scenarios = scenarios,
#'   score = c(0, 0.5, 1, 1.5), score_resp = c(0, 0.25, 1, 3),
#'   n_cohort = 12, cohort_size = 3, window = 30, resp_window = 60,
#'   accrual_rate = 1 / 10, n_trials = 300, seed = 123
#' )
#' res
#'
#' @seealso \code{\link{sim_tite_gboin_et}},
#'   \code{\link{sim_gboin_et_multi}}
#'
#' @export
sim_tite_gboin_et_multi <- function(target, target_resp, scenarios, score,
                                    score_resp, n_cohort, cohort_size,
                                    window, accrual_rate,
                                    resp_window = window,
                                    accrual = c("uniform", "exponential",
                                                "fixed"),
                                    event_time = c("weibull", "uniform"),
                                    late_fraction = 0.5,
                                    resp_late_fraction = 0.5, resp_cor = 0.2,
                                    min_completed_share = 0.5,
                                    n_trials = 10000, start_dose = 1,
                                    n_earlystop = NULL, p_saf = NULL,
                                    p_tox = NULL, p_saf_resp = NULL,
                                    lambda1 = NULL, lambda2 = NULL,
                                    eta1 = NULL, cutoff_eli = 0.95,
                                    cutoff_eli_resp = 0.99, p_eli = NULL,
                                    p_eli_resp = NULL,
                                    boinet_compatible = FALSE,
                                    prior = rep(1 / 6, 6),
                                    est_method = c("observed", "fp_linear",
                                                   "fp_logistic"),
                                    obd_method = c("max_resp",
                                                   "utility_weighted",
                                                   "utility_truncated_linear",
                                                   "utility_scoring"),
                                    min_obd_sample = 1, w1 = 0.33, w2 = 1.09,
                                    t_low = NULL, t_upp = NULL, e_low = NULL,
                                    e_upp = NULL, psi00 = 40, psi11 = 60,
                                    keep_trials = FALSE, verbose = FALSE,
                                    seed = 123) {

  call_expr <- match.call()
  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  est_method <- match.arg(est_method)
  obd_method <- match.arg(obd_method)
  check_flag(verbose, "verbose")

  if (!is.list(scenarios) || length(scenarios) == 0L) {
    stop("'scenarios' must be a non-empty list", call. = FALSE)
  }
  n_scenarios <- length(scenarios)
  nms <- names(scenarios)
  if (is.null(nms)) nms <- rep("", n_scenarios)
  scenario_names <- character(n_scenarios)
  for (i in seq_len(n_scenarios)) {
    element <- scenarios[[i]]
    if (!is.list(element) || !is.matrix(element$p_grade) ||
        !is.matrix(element$p_resp_grade)) {
      stop("each scenario must be a list with the matrices 'p_grade' and ",
           "'p_resp_grade'", call. = FALSE)
    }
    name <- if (nzchar(nms[i])) nms[i] else paste("Scenario", i)
    if (!is.null(element$name)) name <- as.character(element$name)
    scenario_names[i] <- name
  }
  shapes <- vapply(scenarios, function(z) {
    c(dim(z$p_grade), dim(z$p_resp_grade))
  }, integer(4))
  if (any(shapes != shapes[, 1L])) {
    stop("all scenarios must have the same numbers of categories and of ",
         "doses", call. = FALSE)
  }
  if (anyDuplicated(scenario_names)) {
    stop("scenario names must be unique", call. = FALSE)
  }
  n_doses <- ncol(scenarios[[1L]]$p_grade)

  results <- vector("list", n_scenarios)
  names(results) <- scenario_names

  for (i in seq_len(n_scenarios)) {
    if (verbose) {
      message("Scenario ", i, " of ", n_scenarios, ": ", scenario_names[i])
    }
    results[[i]] <- sim_tite_gboin_et(
      target = target, target_resp = target_resp,
      p_grade = scenarios[[i]]$p_grade,
      p_resp_grade = scenarios[[i]]$p_resp_grade, score = score,
      score_resp = score_resp, n_cohort = n_cohort,
      cohort_size = cohort_size, window = window,
      accrual_rate = accrual_rate, resp_window = resp_window,
      accrual = accrual, event_time = event_time,
      late_fraction = late_fraction,
      resp_late_fraction = resp_late_fraction, resp_cor = resp_cor,
      min_completed_share = min_completed_share, n_trials = n_trials,
      start_dose = start_dose, n_earlystop = n_earlystop, p_saf = p_saf,
      p_tox = p_tox, p_saf_resp = p_saf_resp, lambda1 = lambda1,
      lambda2 = lambda2, eta1 = eta1, cutoff_eli = cutoff_eli,
      cutoff_eli_resp = cutoff_eli_resp, p_eli = p_eli,
      p_eli_resp = p_eli_resp, boinet_compatible = boinet_compatible,
      prior = prior, est_method = est_method, obd_method = obd_method,
      min_obd_sample = min_obd_sample, w1 = w1, w2 = w2, t_low = t_low,
      t_upp = t_upp, e_low = e_low, e_upp = e_upp, psi00 = psi00,
      psi11 = psi11, keep_trials = keep_trials, verbose = FALSE, seed = seed
    )
  }

  structure(
    list(
      results = results,
      summary_table = boin_et_multi_table(results, scenario_names, n_doses),
      scenario_names = scenario_names,
      n_doses = n_doses,
      call = call_expr
    ),
    class = "boin_et_oc_multi"
  )
}
