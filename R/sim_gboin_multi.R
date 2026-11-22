#' Operating Characteristics of the gBOIN Design in Several Scenarios
#'
#' @description
#'   Run \code{\link{sim_gboin}} for several sets of true grade probabilities
#'   with the same design, and collect the results in one table.
#'
#' @inheritParams sim_gboin
#'
#' @param scenarios
#'   A list of scenarios, each a matrix of true grade probabilities as the
#'   argument \code{p_grade} of \code{\link{sim_gboin}}, or a list with the
#'   components \code{name} and \code{p_grade}. The names of the list are
#'   used as scenario names. Every scenario must have the same numbers of
#'   toxicity categories and of doses.
#'
#' @return
#'   An object of class \code{boin_oc_multi}, which is a list with components
#'   \item{results}{The result of \code{\link{sim_gboin}} for every
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
#'   \code{\link{sim_boin_multi}}.
#'
#' @references
#'   Mu, R., Yuan, Y., Xu, J., Mandrekar, S. J. and Yin, J. (2019). gBOIN: a
#'   unified model-assisted phase I trial design accounting for toxicity
#'   grades, and binary or continuous end points. Journal of the Royal
#'   Statistical Society: Series C, 68, 289-308. \doi{10.1111/rssc.12263}
#'
#' @examples
#' # Scenarios 8 and 10 of Mu et al. (2019) have the same mean scores but
#' # different grades
#' scenarios <- list(
#'   "Grade 3 only" = rbind(
#'     c(0.92, 0.76, 0.68, 0.57, 0.45, 0.25),
#'     c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
#'     c(0.08, 0.24, 0.32, 0.43, 0.55, 0.75),
#'     c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
#'   ),
#'   "Grade 2, then grade 3" = rbind(
#'     c(0.84, 0.52, 0.36, 0.14, 0.45, 0.25),
#'     c(0.16, 0.48, 0.64, 0.86, 0.00, 0.00),
#'     c(0.00, 0.00, 0.00, 0.00, 0.55, 0.75),
#'     c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
#'   )
#' )
#' res <- sim_gboin_multi(
#'   target = 0.47 / 1.5,
#'   scenarios = scenarios,
#'   score = c(0, 0.5, 1, 1.5),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_earlystop = 100,
#'   n_trials = 300,
#'   seed = 123
#' )
#' res
#'
#' @seealso \code{\link{sim_gboin}}, \code{\link{sim_boin_multi}}
#'
#' @export
sim_gboin_multi <- function(target, scenarios, score, n_cohort, cohort_size,
                            n_trials = 10000, start_dose = 1,
                            n_earlystop = 18, p_saf = NULL, p_tox = NULL,
                            cutoff_eli = 0.95, extrasafe = FALSE,
                            offset = 0.05,
                            elim_rule = c("posterior", "count"),
                            min_mtd_sample = 1, overdose_cutoff = NULL,
                            n_earlystop_rule = c("with_stay", "simple"),
                            keep_trials = FALSE, verbose = FALSE,
                            seed = 123) {

  call_expr <- match.call()
  elim_rule <- match.arg(elim_rule)
  n_earlystop_rule <- match.arg(n_earlystop_rule)

  if (!is.list(scenarios) || length(scenarios) == 0L) {
    stop("'scenarios' must be a non-empty list", call. = FALSE)
  }
  n_scenarios <- length(scenarios)
  nms <- names(scenarios)
  if (is.null(nms)) nms <- rep("", n_scenarios)
  scenario_names <- character(n_scenarios)
  grades <- vector("list", n_scenarios)
  for (i in seq_len(n_scenarios)) {
    element <- scenarios[[i]]
    name <- if (nzchar(nms[i])) nms[i] else paste("Scenario", i)
    if (is.matrix(element)) {
      grades[[i]] <- element
    } else if (is.list(element) && is.matrix(element$p_grade)) {
      grades[[i]] <- element$p_grade
      if (!is.null(element$name)) name <- as.character(element$name)
    } else {
      stop("each scenario must be a matrix of grade probabilities or a ",
           "list with a 'p_grade' element", call. = FALSE)
    }
    scenario_names[i] <- name
  }
  shapes <- vapply(grades, dim, integer(2))
  if (any(shapes != shapes[, 1L])) {
    stop("all scenarios must have the same numbers of toxicity categories ",
         "and of doses", call. = FALSE)
  }
  if (anyDuplicated(scenario_names)) {
    stop("scenario names must be unique", call. = FALSE)
  }
  n_doses <- ncol(grades[[1L]])

  results <- vector("list", n_scenarios)
  names(results) <- scenario_names

  for (i in seq_len(n_scenarios)) {
    if (verbose) {
      message("Scenario ", i, " of ", n_scenarios, ": ", scenario_names[i])
    }
    results[[i]] <- sim_gboin(
      target = target, p_grade = grades[[i]], score = score,
      n_cohort = n_cohort, cohort_size = cohort_size, n_trials = n_trials,
      start_dose = start_dose, n_earlystop = n_earlystop, p_saf = p_saf,
      p_tox = p_tox, cutoff_eli = cutoff_eli, extrasafe = extrasafe,
      offset = offset, elim_rule = elim_rule,
      min_mtd_sample = min_mtd_sample, overdose_cutoff = overdose_cutoff,
      n_earlystop_rule = n_earlystop_rule, keep_trials = keep_trials,
      verbose = FALSE, seed = seed
    )
  }

  summary_table <- oc_multi_table(results, scenario_names, n_doses)

  structure(
    list(
      results = results,
      summary_table = summary_table,
      scenario_names = scenario_names,
      n_doses = n_doses,
      call = call_expr
    ),
    class = "boin_oc_multi"
  )
}
