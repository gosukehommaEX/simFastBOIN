#' Select the MTD of an iBOIN Trial
#'
#' @description
#'   Select the maximum tolerated dose from the final data of one or more
#'   trials of the iBOIN design of Zhou et al. (2021), with or without the
#'   prior information in the estimates of the DLT probabilities.
#'
#' @param n_pts
#'   Integer vector or matrix. Number of patients treated at each dose. A
#'   matrix holds one trial per row.
#'
#' @param n_tox
#'   Integer vector or matrix of the same shape. Number of DLTs observed at
#'   each dose.
#'
#' @param target
#'   Numeric scalar. Target DLT probability.
#'
#' @param skeleton
#'   Numeric vector, one value per dose. Prior DLT probability at each dose.
#'
#' @param prior_ess
#'   Whole number or vector of whole numbers, one per dose. Prior effective
#'   sample size.
#'
#' @param with_prior
#'   Logical scalar. When \code{TRUE}, the default, the prior information is
#'   added to the data before the isotonic regression. When \code{FALSE} the
#'   selection is that of \code{\link{boin_select_mtd}}.
#'
#' @param cutoff_eli
#'   Numeric scalar. Posterior probability cutoff for dose elimination.
#'   Defaults to 0.95.
#'
#' @param extrasafe
#'   Logical scalar. Apply the stricter safety stopping rule at the lowest
#'   dose. Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is
#'   relaxed for the safety stopping rule. Defaults to 0.05.
#'
#' @param robust
#'   Logical scalar. Use the robust prior, see \code{\link{iboin_boundary}}.
#'   Defaults to \code{FALSE}.
#'
#' @param min_mtd_sample
#'   Integer scalar. Smallest number of patients a dose must have been given
#'   to be selected. Defaults to 1.
#'
#' @return
#'   A data frame with one row per trial and the columns \code{trial},
#'   \code{mtd} (the selected dose, \code{NA} when none is selected) and
#'   \code{reason} (\code{"selected"}, \code{"lowest_dose_eliminated"} or
#'   \code{"no_admissible_dose"}).
#'
#' @details
#'   The article gives two estimates of the DLT probability of a dose for the
#'   selection of the MTD: the observed rate \code{n_tox / n_pts}, and the
#'   rate with the prior information,
#'   \code{(n_tox + prior_ess * skeleton) / (n_pts + prior_ess)}. Either is
#'   made monotone by isotonic regression, and the dose whose estimate is
#'   closest to the target is selected. With \code{with_prior = TRUE} the
#'   numbers \code{n_tox + prior_ess * skeleton} and \code{n_pts + prior_ess}
#'   take the place of the DLTs and the patients in the estimate of
#'   \code{\link{boin_select_mtd}}; the doses that can be selected are still
#'   those that have been given to patients and are not eliminated, and the
#'   elimination uses the data alone, under a uniform prior.
#'
#'   The article does not say which estimate its simulation study uses. Its
#'   operating characteristics (Table S.1 of its supporting information) are
#'   reproduced with the prior information and not without it, see
#'   \code{inst/reproduce/reproduce_Zhou2021.R}. The estimate with the prior
#'   selects the true MTD more often when the prior MTD is right and less
#'   often when it is wrong.
#'
#' @references
#'   Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
#'   Incorporating historical information to improve phase I clinical trials.
#'   Pharmaceutical Statistics, 20, 1017-1034. \doi{10.1002/pst.2121}
#'
#' @examples
#' skeleton <- c(0.10, 0.19, 0.30, 0.42, 0.54)
#' n_pts <- c(3, 6, 12, 9, 0)
#' n_tox <- c(0, 1, 3, 4, 0)
#'
#' iboin_select_mtd(n_pts, n_tox, target = 0.30, skeleton = skeleton,
#'                  prior_ess = 3)
#' iboin_select_mtd(n_pts, n_tox, target = 0.30, skeleton = skeleton,
#'                  prior_ess = 3, with_prior = FALSE)
#'
#' @seealso \code{\link{sim_iboin}}, \code{\link{boin_select_mtd}}
#'
#' @export
iboin_select_mtd <- function(n_pts, n_tox, target, skeleton, prior_ess,
                             with_prior = TRUE, cutoff_eli = 0.95,
                             extrasafe = FALSE, offset = 0.05,
                             robust = FALSE, min_mtd_sample = 1) {

  n_pts <- as_count_matrix(n_pts, "n_pts")
  n_tox <- as_count_matrix(n_tox, "n_tox")

  if (!identical(dim(n_pts), dim(n_tox))) {
    stop("'n_pts' and 'n_tox' must have the same dimensions", call. = FALSE)
  }
  if (any(n_tox > n_pts)) {
    stop("'n_tox' must not exceed 'n_pts' at any dose", call. = FALSE)
  }
  check_flag(with_prior, "with_prior")
  check_flag(extrasafe, "extrasafe")
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_count(min_mtd_sample, "min_mtd_sample", 1L)
  check_scalar_prob(target, "target")
  if (length(skeleton) != ncol(n_pts)) {
    stop("'skeleton' must have one value for each dose", call. = FALSE)
  }
  prior <- iboin_prior(target, skeleton, prior_ess, 0.6 * target,
                       1.4 * target, robust = robust)

  ess <- if (with_prior) prior$prior_ess else rep(0L, ncol(n_pts))
  est_y <- n_tox + matrix(ess * skeleton, nrow(n_pts), ncol(n_pts),
                          byrow = TRUE)
  est_n <- n_pts + matrix(ess, nrow(n_pts), ncol(n_pts), byrow = TRUE)
  storage.mode(est_y) <- "double"
  storage.mode(est_n) <- "double"

  res <- select_mtd_estimate_cpp(
    eliminated = mtd_eliminated(n_pts, n_tox, target, cutoff_eli, extrasafe,
                                offset),
    n_pts = n_pts,
    est_y = est_y,
    est_n = est_n,
    target = as.numeric(target),
    min_mtd_sample = as.integer(min_mtd_sample)
  )

  reasons <- c("selected", "lowest_dose_eliminated", "no_admissible_dose")

  data.frame(
    trial = seq_len(nrow(n_pts)),
    mtd = res$mtd,
    reason = reasons[res$reason + 1L],
    stringsAsFactors = FALSE
  )
}
