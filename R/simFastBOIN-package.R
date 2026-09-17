#' simFastBOIN: Fast Simulation of Bayesian Optimal Interval Designs
#'
#' @description
#'   Tools for designing and evaluating phase I dose-finding trials that use the
#'   Bayesian optimal interval (BOIN) design of Liu and Yuan (2015). The
#'   simulation engine is written in C++ and reproduces the reference
#'   implementation in the \pkg{BOIN} package exactly, including the order in
#'   which random variates are consumed, so that results obtained with the same
#'   seed agree to the last trial.
#'
#' @section Design tools:
#'   \describe{
#'     \item{\code{\link{boin_lambda}}}{Escalation and de-escalation interval boundaries.}
#'     \item{\code{\link{boin_boundary}}}{Integer decision boundaries as a function of sample size.}
#'     \item{\code{\link{boin_decision_table}}}{Decision table indexed by DLTs and patients.}
#'   }
#'
#' @section Simulation tools:
#'   \describe{
#'     \item{\code{\link{sim_boin}}}{Operating characteristics for one scenario.}
#'     \item{\code{\link{sim_boin_multi}}}{Operating characteristics for several scenarios.}
#'     \item{\code{\link{boin_simulate}}}{Raw trial data from the simulation engine.}
#'     \item{\code{\link{boin_isotonic}}}{Isotonic estimates of the DLT probability.}
#'     \item{\code{\link{boin_select_mtd}}}{MTD selection from completed trials.}
#'   }
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#'   Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020). BOIN: An R
#'   Package for Designing Single-Agent and Drug-Combination Dose-Finding Trials Using
#'   Bayesian Optimal Interval Designs. Journal of Statistical Software, 94(13), 1-32.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib simFastBOIN, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats pbeta
## usethis namespace: end
NULL
