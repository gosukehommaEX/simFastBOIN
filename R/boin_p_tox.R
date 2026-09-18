#' Toxic Threshold Giving a Required De-escalation Boundary
#'
#' @description
#'   Find the value of \code{p_tox} for which the BOIN de-escalation boundary
#'   \code{lambda_d} equals a required value. This inverts
#'   \code{\link{boin_lambda}}, which maps \code{p_tox} to \code{lambda_d}.
#'
#' @param target
#'   Numeric scalar. Target DLT probability.
#'
#' @param lambda_d
#'   Numeric scalar. Required de-escalation boundary. Must lie strictly between
#'   \code{target} and 1.
#'
#' @return
#'   A numeric scalar, the value of \code{p_tox} that produces \code{lambda_d}.
#'
#' @details
#'   \code{lambda_d} increases strictly with \code{p_tox}, from \code{target} in
#'   the limit as \code{p_tox} approaches \code{target} to 1 in the limit as
#'   \code{p_tox} approaches 1, so the solution exists and is unique for any
#'   \code{lambda_d} strictly between those bounds. The root is found
#'   numerically.
#'
#'   A boundary close to the target needs a toxic threshold close to the target,
#'   and the rest of the package refuses a \code{p_tox} within ten percent of
#'   \code{target}, on the grounds that phase I sample sizes cannot tell the two
#'   rates apart. The value is still returned, since the inversion is well
#'   defined, but a warning names the smallest boundary that leaves a usable
#'   threshold. For a target of 0.30 that smallest boundary is about 0.315.
#'
#'   A design team that wants to tighten the de-escalation boundary usually
#'   states the boundary rather than the threshold behind it. This function turns
#'   that statement into the \code{p_tox} to pass to the other functions, so that
#'   the whole design follows from it.
#'
#' @examples
#' # A target of 0.30 gives a de-escalation boundary of 0.359 by default
#' boin_lambda(target = 0.30)$lambda_d
#'
#' # The toxic threshold that tightens it to 0.33
#' p_tox <- boin_p_tox(target = 0.30, lambda_d = 0.33)
#' p_tox
#'
#' # which is what it claims to be
#' boin_lambda(target = 0.30, p_tox = p_tox)$lambda_d
#'
#' # A boundary too close to the target leaves an unusable threshold, and says so
#' boin_p_tox(target = 0.30, lambda_d = 0.31)
#'
#' # and shifts the de-escalation boundary down by one DLT at every cohort end
#' as.data.frame(boin_boundary(0.30, 18))[seq(3, 18, 3), ]
#' as.data.frame(boin_boundary(0.30, 18, p_tox = p_tox))[seq(3, 18, 3), ]
#'
#' @seealso \code{\link{boin_lambda}}, \code{\link{boin_boundary}}
#'
#' @importFrom stats uniroot
#'
#' @export
boin_p_tox <- function(target, lambda_d) {

  check_scalar_prob(target, "target")
  check_scalar_prob(lambda_d, "lambda_d")
  if (lambda_d <= target) {
    stop("'lambda_d' must be greater than 'target'; the de-escalation boundary ",
         "always lies above the target rate", call. = FALSE)
  }

  boundary_at <- function(p_tox) {
    log((1 - target) / (1 - p_tox)) /
      log(p_tox * (1 - target) / (target * (1 - p_tox)))
  }

  lower <- target + 1e-6
  upper <- 1 - 1e-9
  if (lambda_d >= boundary_at(upper)) {
    stop("'lambda_d' is too close to 1 to be attained for this target",
         call. = FALSE)
  }

  p_tox <- uniroot(
    function(p_tox) boundary_at(p_tox) - lambda_d,
    interval = c(lower, upper),
    tol = .Machine$double.eps^0.5
  )$root

  # The other functions reject a threshold within ten percent of the target, so
  # say so here rather than letting the caller meet the error further on.
  if ((p_tox - target) < 0.1 * target) {
    warning("a de-escalation boundary of ", format(lambda_d), " needs 'p_tox' = ",
            format(signif(p_tox, 6)), ", which is too close to 'target' for the ",
            "other functions to accept. The smallest usable boundary for a ",
            "target of ", format(target), " is about ",
            format(signif(boundary_at(1.1 * target), 4)), ".", call. = FALSE)
  }

  p_tox
}
