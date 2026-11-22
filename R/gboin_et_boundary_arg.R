#' Boundaries of a gBOIN-ET Trial
#'
#' @description
#'   Return the boundaries of the gBOIN-ET functions: those of
#'   \code{\link{boin_et_boundary}} for the targets, or the three values
#'   given by the user. Internal.
#'
#' @param target,target_resp,p_saf,p_tox,p_saf_resp,prior
#'   Arguments of \code{\link{boin_et_boundary}}.
#'
#' @param lambda1,lambda2,eta1
#'   Numeric scalars or \code{NULL}. Boundaries given directly; all three or
#'   none.
#'
#' @return
#'   An object of class \code{boin_et_boundary}. When the boundaries are
#'   given, \code{prob_incorrect} is \code{NA}.
#'
#' @keywords internal
#' @noRd
gboin_et_boundary_arg <- function(target, target_resp, p_saf, p_tox,
                                  p_saf_resp, prior, lambda1, lambda2, eta1) {
  given <- c(!is.null(lambda1), !is.null(lambda2), !is.null(eta1))
  if (!any(given)) {
    return(boin_et_boundary(target, target_resp, p_saf = p_saf, p_tox = p_tox,
                            p_saf_resp = p_saf_resp, prior = prior))
  }
  if (!all(given)) {
    stop("give all of 'lambda1', 'lambda2' and 'eta1', or none of them",
         call. = FALSE)
  }
  check_scalar_prob(lambda1, "lambda1")
  check_scalar_prob(lambda2, "lambda2")
  check_scalar_prob(eta1, "eta1")
  if (!(lambda1 < lambda2)) {
    stop("'lambda1' must be below 'lambda2'", call. = FALSE)
  }
  par <- boin_et_check_design(target, target_resp, p_saf, p_tox, p_saf_resp)
  structure(
    list(
      lambda1 = lambda1,
      lambda2 = lambda2,
      eta1 = eta1,
      prob_incorrect = NA_real_,
      target = par$target,
      target_resp = par$target_resp,
      p_saf = par$p_saf,
      p_tox = par$p_tox,
      p_saf_resp = par$p_saf_resp,
      prior = prior
    ),
    class = "boin_et_boundary"
  )
}
