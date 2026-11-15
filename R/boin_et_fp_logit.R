#' Fractional Polynomial Estimate of the Response Rates
#'
#' @description
#'   Estimate the response rates of the doses tried in a BOIN-ET trial with a
#'   logistic regression on a fractional polynomial of the dose level with two
#'   degrees of freedom, as in Takeda et al. (2018). This follows
#'   \code{boinet::fp.logit()} (boinet 1.6.0, MIT license, by Yusuke Yamaguchi
#'   and Kentaro Takeda): the observed rates are fitted by
#'   \code{mfp::mfp()} with a quasi-binomial family, one observation per dose.
#'   Internal.
#'
#' @param obs
#'   Integer vector. Responses at each dose tried, in dose order.
#'
#' @param n
#'   Integer vector. Patients at each dose tried.
#'
#' @param dose
#'   Integer vector. The dose levels tried.
#'
#' @return
#'   Numeric vector of the fitted response rates.
#'
#' @keywords internal
#' @noRd
boin_et_fp_logit <- function(obs, n, dose) {
  if (!requireNamespace("mfp", quietly = TRUE)) {
    stop("Package 'mfp' is needed for est_method = \"fp_logistic\"",
         call. = FALSE)
  }
  prob <- obs / n
  obse_df <- data.frame(prob = prob, dose = dose)
  # fp() is a special of mfp::mfp() that the model frame evaluates in the
  # environment of the formula, so the formula gets one where fp() is found
  # without attaching 'mfp'.
  env <- new.env(parent = environment())
  assign("fp", mfp::fp, envir = env)
  model <- stats::as.formula(
    "prob ~ fp(dose, df = 4, select = 0.99999, alpha = 0.99999)", env = env
  )
  fpfit <- suppressWarnings(mfp::mfp(
    model,
    family = stats::quasi(link = "logit", variance = "mu(1-mu)"),
    data = obse_df
  ))
  as.numeric(fpfit$fit$fitted.values[match(dose, obse_df$dose)])
}
