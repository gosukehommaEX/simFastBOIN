#' True Dose-Toxicity and Dose-Response Scenarios
#'
#' @description
#'   Tabulate a list of scenarios, as given to the \code{_multi} simulation
#'   functions, with the true utility of every dose and the doses that are
#'   the true maximum tolerated dose (MTD) and the true optimal biological
#'   dose (OBD). The result is the truth against which
#'   \code{\link{oc_compare}} judges the designs, and it has
#'   \code{print()} and \code{plot()} methods.
#'
#' @param scenarios
#'   List of scenarios. For a single agent, each element is a numeric vector
#'   of DLT probabilities or a list with the component \code{p_true}, and
#'   optionally \code{p_resp} (response probabilities) and \code{name}. For
#'   two drugs given together, \code{p_true} (and \code{p_resp}) are matrices
#'   with one row per level of drug A and one column per level of drug B.
#'   Either every scenario or none has \code{p_resp}. Names of the list are
#'   used when \code{name} is absent, and "Scenario 1", "Scenario 2", ...
#'   otherwise.
#'
#' @param target
#'   Numeric scalar. Target DLT probability; the true MTD is the dose whose
#'   DLT probability is closest to it.
#'
#' @param target_resp
#'   Numeric scalar or \code{NULL}. Lowest acceptable response probability,
#'   as \code{target_resp} of \code{\link{sim_boin12}}; a dose with a lower
#'   response probability cannot be the true OBD. Defaults to \code{NULL}, no
#'   such limit. The \code{target_resp} of \code{\link{sim_boin_et}}, the
#'   target response probability of BOIN-ET, is a different quantity.
#'
#' @param tox_limit
#'   Numeric scalar. Highest acceptable DLT probability; a dose with a higher
#'   DLT probability cannot be the true OBD. Defaults to \code{target}.
#'
#' @param psi00,psi11
#'   Numeric scalars. Utilities of a patient with neither a DLT nor a response
#'   and with both, on the scale where a response without a DLT scores 100
#'   and a DLT without a response scores 0. Default to 40 and 60, as in
#'   \code{\link{sim_boin12}}.
#'
#' @param resp_cor
#'   Numeric scalar. Correlation between the DLT and the response, used for
#'   the probabilities of the four outcomes. Defaults to 0, independence.
#'
#' @param resp_cor_type
#'   Character scalar. How \code{resp_cor} is read, as in
#'   \code{\link{sim_boin12}}: \code{"copula"} (the default), the correlation
#'   of the latent normal variables, or \code{"binary"}, the correlation of
#'   the two binary outcomes.
#'
#' @param true_mtd,true_obd
#'   \code{NULL} (the default) or a list with one element per scenario, in
#'   the order of \code{scenarios} or named by the scenarios, giving the true
#'   MTDs or OBDs directly when the definition of an article differs from the
#'   one used here: a vector of dose levels for a single agent, a two-column
#'   matrix of the levels of drug A and drug B for a combination, or
#'   \code{NULL} for the computed doses. An empty vector means no such dose.
#'
#' @return
#'   A data frame of class \code{c("oc_scenarios", "data.frame")} with one row
#'   per scenario and dose (or combination) and the columns
#'   \item{scenario}{Factor of the scenario names, in the order given.}
#'   \item{dose}{Dose level, or \code{dose_a} and \code{dose_b} for a
#'     combination.}
#'   \item{p_true}{True DLT probability.}
#'   \item{p_resp, utility}{True response probability and utility, when the
#'     scenarios have response probabilities.}
#'   \item{true_mtd}{Logical, the dose is a true MTD.}
#'   \item{true_obd}{Logical, the dose is a true OBD, when the scenarios have
#'     response probabilities.}
#'   The arguments are kept as attributes, with \code{type} (\code{"single"}
#'   or \code{"combination"}) and \code{dims} (the number of doses, or the
#'   numbers of levels of the two drugs).
#'
#' @details
#'   The utility of a dose is
#'   \code{100 * Pr(response, no DLT) + psi00 * Pr(neither) + psi11 *
#'   Pr(both)}. The true MTDs are all the doses whose DLT probability is
#'   closest to \code{target} (several when tied). The true OBDs are all the
#'   doses with the highest utility among the doses whose DLT probability is
#'   at most \code{tox_limit} and, when \code{target_resp} is given, whose
#'   response probability is at least \code{target_resp}; a scenario may have
#'   none. Ties are decided with a tolerance of \code{1e-10}.
#'
#'   Some articles use other definitions, for example Lin et al. (2020), who
#'   count every dose within 3 points of the highest utility as an OBD; give
#'   such doses through \code{true_obd}.
#'
#' @references
#'   Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
#'   optimal interval phase I/II trial design for utility-based dose finding
#'   in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
#'   1393-1402. \doi{10.1200/PO.20.00257}
#'
#' @examples
#' truth <- oc_scenarios(
#'   list(
#'     "Scenario 1" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                         p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
#'     "Scenario 2" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'                         p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
#'   ),
#'   target = 0.35, target_resp = 0.25
#' )
#' truth
#'
#' @seealso \code{\link{oc_results}}, \code{\link{oc_compare}},
#'   \code{\link{print.oc_scenarios}}, \code{\link{plot.oc_scenarios}}
#'
#' @export
oc_scenarios <- function(scenarios, target, target_resp = NULL,
                         tox_limit = target, psi00 = 40, psi11 = 60,
                         resp_cor = 0,
                         resp_cor_type = c("copula", "binary"),
                         true_mtd = NULL, true_obd = NULL) {

  resp_cor_type <- match.arg(resp_cor_type)
  check_scalar_prob(target, "target")
  check_scalar_prob(tox_limit, "tox_limit")
  if (!is.null(target_resp)) check_scalar_prob(target_resp, "target_resp")
  for (nm in c("psi00", "psi11")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0 ||
        v > 100) {
      stop("'", nm, "' must be a single number between 0 and 100",
           call. = FALSE)
    }
  }
  if (!is.numeric(resp_cor) || length(resp_cor) != 1L ||
      !is.finite(resp_cor) || resp_cor <= -1 || resp_cor >= 1) {
    stop("'resp_cor' must be a single number strictly between -1 and 1",
         call. = FALSE)
  }

  if (!is.list(scenarios) || length(scenarios) == 0L) {
    stop("'scenarios' must be a non-empty list", call. = FALSE)
  }
  first <- scenarios[[1L]]
  first_p <- if (is.list(first)) first[["p_true"]] else first
  combination <- is.matrix(first_p)
  has_resp <- vapply(scenarios, function(e) {
    is.list(e) && !is.null(e[["p_resp"]])
  }, logical(1))
  if (any(has_resp) && !all(has_resp)) {
    stop("either every scenario or none must have 'p_resp'", call. = FALSE)
  }
  has_resp <- all(has_resp)

  if (combination) {
    norm <- normalize_comb_scenarios(scenarios)
    dims <- dim(norm[[1L]]$p_true)
    if (has_resp) {
      for (i in seq_along(norm)) {
        p_resp <- scenarios[[i]][["p_resp"]]
        if (!is.numeric(p_resp) || !is.matrix(p_resp) ||
            !identical(dim(p_resp), dims) || any(!is.finite(p_resp)) ||
            any(p_resp < 0) || any(p_resp > 1)) {
          stop("'p_resp' of scenario '", norm[[i]]$name, "' must be a ",
               "matrix of response probabilities of the same size as ",
               "'p_true'", call. = FALSE)
        }
        norm[[i]]$p_resp <- p_resp
      }
    }
  } else {
    norm <- if (has_resp) normalize_backfill_scenarios(scenarios) else
      normalize_scenarios(scenarios)
    dims <- length(norm[[1L]]$p_true)
  }
  scenario_names <- vapply(norm, function(z) z$name, character(1))

  if (!has_resp && !is.null(true_obd)) {
    stop("'true_obd' needs scenarios with 'p_resp'", call. = FALSE)
  }
  given_mtd <- oc_true_override(true_mtd, "true_mtd", scenario_names, dims)
  given_obd <- oc_true_override(true_obd, "true_obd", scenario_names, dims)

  tol <- 1e-10
  rows <- vector("list", length(norm))
  for (i in seq_along(norm)) {
    p_true <- as.vector(norm[[i]]$p_true)
    gap <- abs(p_true - target)
    mtd <- gap <= min(gap) + tol
    if (!is.null(given_mtd[[i]])) mtd <- as.vector(given_mtd[[i]])
    if (combination) {
      df <- data.frame(dose_a = rep(seq_len(dims[1L]), times = dims[2L]),
                       dose_b = rep(seq_len(dims[2L]), each = dims[1L]))
    } else {
      df <- data.frame(dose = seq_len(dims))
    }
    df$p_true <- p_true
    if (has_resp) {
      p_resp <- as.vector(norm[[i]]$p_resp)
      pr <- oc_outcome_probs(p_true, p_resp, resp_cor, resp_cor_type)
      utility <- 100 * pr$resp_no_dlt + psi00 * pr$neither +
        psi11 * pr$both
      ok <- p_true <= tox_limit + tol
      if (!is.null(target_resp)) ok <- ok & p_resp >= target_resp - tol
      obd <- logical(length(p_true))
      if (any(ok)) obd <- ok & utility >= max(utility[ok]) - tol
      if (!is.null(given_obd[[i]])) obd <- as.vector(given_obd[[i]])
      df$p_resp <- p_resp
      df$utility <- utility
    }
    df$true_mtd <- mtd
    if (has_resp) df$true_obd <- obd
    rows[[i]] <- cbind(
      data.frame(scenario = rep(scenario_names[i], nrow(df)),
                 stringsAsFactors = FALSE),
      df
    )
  }
  out <- do.call(rbind, rows)
  out$scenario <- factor(out$scenario, levels = scenario_names)
  rownames(out) <- NULL

  class(out) <- c("oc_scenarios", "data.frame")
  attr(out, "type") <- if (combination) "combination" else "single"
  attr(out, "dims") <- dims
  attr(out, "target") <- target
  attr(out, "target_resp") <- target_resp
  attr(out, "tox_limit") <- tox_limit
  attr(out, "psi00") <- psi00
  attr(out, "psi11") <- psi11
  attr(out, "resp_cor") <- resp_cor
  attr(out, "resp_cor_type") <- resp_cor_type
  out
}
