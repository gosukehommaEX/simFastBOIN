#' Validate the Design Parameters of U-BOIN
#'
#' @description
#'   Check the limits, the stage I target, the utilities, the cutoffs and the
#'   priors of the U-BOIN design, and fill in the defaults. Internal.
#'
#' @param tox_upper
#'   Numeric scalar. Highest acceptable DLT probability.
#'
#' @param resp_lower
#'   Numeric scalar. Lowest acceptable response probability.
#'
#' @param target
#'   Numeric scalar or \code{NULL}. Target DLT probability of stage I.
#'
#' @param p_saf
#'   Numeric scalar or \code{NULL}. Highest DLT probability deemed
#'   subtherapeutic, for the boundaries of stage I.
#'
#' @param p_tox
#'   Numeric scalar or \code{NULL}. Lowest DLT probability deemed overly
#'   toxic, for the boundaries of stage I.
#'
#' @param utility
#'   Numeric vector of length four.
#'
#' @param cutoff_eli
#'   Numeric scalar. Cutoff of the toxicity rule.
#'
#' @param cutoff_eli_resp
#'   Numeric scalar. Cutoff of the futility rule.
#'
#' @param prior_tox
#'   Numeric vector of length two. Beta prior of the toxicity rule.
#'
#' @param prior_resp
#'   Numeric vector of length two. Beta prior of the futility rule.
#'
#' @param stage1
#'   Logical scalar. Check the target and derive the boundaries of stage I;
#'   \code{FALSE} for the selection of the OBD, which does not use them.
#'
#' @return
#'   A list with the checked values and \code{lambda_e}, \code{lambda_d}
#'   (\code{NA} when \code{stage1} is \code{FALSE}).
#'
#' @keywords internal
#' @noRd
uboin_check_design <- function(tox_upper, resp_lower, target, p_saf, p_tox,
                               utility, cutoff_eli, cutoff_eli_resp,
                               prior_tox, prior_resp, stage1 = TRUE) {
  check_scalar_prob(tox_upper, "tox_upper")
  check_scalar_prob(resp_lower, "resp_lower")
  lambda <- list(lambda_e = NA_real_, lambda_d = NA_real_)
  if (stage1) {
    if (is.null(target)) target <- tox_upper - 0.05
    check_scalar_prob(target, "target")
    if (target > tox_upper) {
      stop("'target' must not exceed 'tox_upper'", call. = FALSE)
    }
    if (is.null(p_saf)) p_saf <- 0.6 * target
    if (is.null(p_tox)) p_tox <- 1.4 * target
    lambda <- boin_lambda(target, p_saf = p_saf, p_tox = p_tox)
  }
  if (!is.numeric(utility) || length(utility) != 4L ||
      any(!is.finite(utility))) {
    stop("'utility' must hold four finite numbers, for the outcomes ",
         "(no response, DLT), (no response, no DLT), (response, DLT) and ",
         "(response, no DLT)", call. = FALSE)
  }
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  for (prior in list(list(prior_tox, "prior_tox"),
                     list(prior_resp, "prior_resp"))) {
    if (!is.numeric(prior[[1L]]) || length(prior[[1L]]) != 2L ||
        any(!is.finite(prior[[1L]])) || any(prior[[1L]] <= 0)) {
      stop("'", prior[[2L]], "' must hold the two positive parameters of a ",
           "beta distribution", call. = FALSE)
    }
  }
  list(tox_upper = tox_upper, resp_lower = resp_lower, target = target,
       p_saf = p_saf, p_tox = p_tox, lambda_e = lambda$lambda_e,
       lambda_d = lambda$lambda_d, utility = as.numeric(utility),
       cutoff_eli = cutoff_eli, cutoff_eli_resp = cutoff_eli_resp,
       prior_tox = as.numeric(prior_tox), prior_resp = as.numeric(prior_resp))
}
