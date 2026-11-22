#' Operating Characteristics of the gBOIN-ET Design in Several Scenarios
#'
#' @description
#'   Run \code{\link{sim_gboin_et}} for several scenarios of probabilities of
#'   toxicity and efficacy categories with the same design, and collect the
#'   results in one table.
#'
#' @inheritParams sim_gboin_et
#'
#' @param scenarios
#'   A list of scenarios, each a list with the components \code{p_grade} and
#'   \code{p_resp_grade}, matrices as the arguments of the same names of
#'   \code{\link{sim_gboin_et}}, and optionally \code{name}. The names of the
#'   list are used as scenario names. Every scenario must have the same
#'   numbers of categories and of doses.
#'
#' @return
#'   An object of class \code{boin_et_oc_multi}, which is a list with
#'   components
#'   \item{results}{The result of \code{\link{sim_gboin_et}} for every
#'     scenario.}
#'   \item{summary_table}{Data frame with the true mean scores, the
#'     selection percentages and the average numbers of patients and total
#'     scores of every scenario.}
#'   \item{scenario_names}{The names of the scenarios.}
#'   \item{n_doses}{The number of doses.}
#'   \item{call}{The call.}
#'
#' @details
#'   Every scenario is simulated with the same seed, as in
#'   \code{\link{sim_boin_et_multi}}.
#'
#' @references
#'   Takeda, K., Morita, S. and Taguri, M. (2022). gBOIN-ET: The generalized
#'   Bayesian optimal interval design for optimal dose-finding accounting for
#'   ordinal graded efficacy and toxicity in early clinical trials.
#'   Biometrical Journal, 64(7), 1178-1191. \doi{10.1002/bimj.202100263}
#'
#' @examples
#' # The same toxicity grades with an increasing and with a plateauing
#' # efficacy
#' p_grade <- rbind(
#'   c(0.83, 0.62, 0.51, 0.34),
#'   c(0.12, 0.18, 0.19, 0.16),
#'   c(0.04, 0.11, 0.14, 0.15),
#'   c(0.01, 0.09, 0.16, 0.35)
#' )
#' scenarios <- list(
#'   "Increasing" = list(
#'     p_grade = p_grade,
#'     p_resp_grade = rbind(
#'       c(0.90, 0.50, 0.30, 0.25),
#'       c(0.05, 0.35, 0.20, 0.20),
#'       c(0.05, 0.15, 0.50, 0.55)
#'     )
#'   ),
#'   "Plateau" = list(
#'     p_grade = p_grade,
#'     p_resp_grade = rbind(
#'       c(0.60, 0.30, 0.30, 0.30),
#'       c(0.30, 0.20, 0.20, 0.20),
#'       c(0.10, 0.50, 0.50, 0.50)
#'     )
#'   )
#' )
#' res <- sim_gboin_et_multi(
#'   target = 0.47 / 1.5, target_resp = 0.95 / 1.5, scenarios = scenarios,
#'   score = c(0, 0.5, 1, 1.5), score_resp = c(0, 1, 1.5),
#'   n_cohort = 10, cohort_size = 3, n_trials = 300, seed = 123
#' )
#' res
#'
#' @seealso \code{\link{sim_gboin_et}}, \code{\link{sim_boin_et_multi}}
#'
#' @export
sim_gboin_et_multi <- function(target, target_resp, scenarios, score,
                               score_resp, n_cohort, cohort_size,
                               n_trials = 10000, start_dose = 1,
                               n_earlystop = NULL, p_saf = NULL, p_tox = NULL,
                               p_saf_resp = NULL, lambda1 = NULL,
                               lambda2 = NULL, eta1 = NULL, cutoff_eli = 0.95,
                               cutoff_eli_resp = 0.99, p_eli = NULL,
                               p_eli_resp = NULL, resp_cor = 0.2,
                               boinet_compatible = FALSE,
                               prior = rep(1 / 6, 6),
                               est_method = c("observed", "fp_linear",
                                              "fp_logistic"),
                               obd_method = c("max_resp", "utility_weighted",
                                              "utility_truncated_linear",
                                              "utility_scoring"),
                               min_obd_sample = 1, w1 = 0.33, w2 = 1.09,
                               t_low = NULL, t_upp = NULL, e_low = NULL,
                               e_upp = NULL, psi00 = 40, psi11 = 60,
                               keep_trials = FALSE, verbose = FALSE,
                               seed = 123) {

  call_expr <- match.call()
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
    results[[i]] <- sim_gboin_et(
      target = target, target_resp = target_resp,
      p_grade = scenarios[[i]]$p_grade,
      p_resp_grade = scenarios[[i]]$p_resp_grade, score = score,
      score_resp = score_resp, n_cohort = n_cohort,
      cohort_size = cohort_size, n_trials = n_trials,
      start_dose = start_dose, n_earlystop = n_earlystop, p_saf = p_saf,
      p_tox = p_tox, p_saf_resp = p_saf_resp, lambda1 = lambda1,
      lambda2 = lambda2, eta1 = eta1, cutoff_eli = cutoff_eli,
      cutoff_eli_resp = cutoff_eli_resp, p_eli = p_eli,
      p_eli_resp = p_eli_resp, resp_cor = resp_cor,
      boinet_compatible = boinet_compatible, prior = prior,
      est_method = est_method, obd_method = obd_method,
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
