#' Select the MTD of a gBOIN Trial
#'
#' @description
#'   Select the maximum tolerated dose from the final data of one or more
#'   trials of the gBOIN design of Mu et al. (2019), in which the toxicity
#'   of a dose is the total of the normalized scores of its patients.
#'
#' @param n_pts
#'   Integer vector or matrix. Number of patients treated at each dose. A
#'   matrix holds one trial per row.
#'
#' @param score_sum
#'   Numeric vector or matrix of the same shape. Total of the normalized
#'   toxicity scores, each between 0 and 1, of the patients treated at each
#'   dose.
#'
#' @param target
#'   Numeric scalar. Target value of the mean normalized toxicity score.
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
#'   The selection is that of \code{\link{boin_select_mtd}} with the total
#'   score in the place of the number of DLTs: doses whose posterior
#'   probability of a mean score above the target exceeds \code{cutoff_eli}
#'   are excluded together with the doses above them, the mean scores of
#'   the remaining doses that were given to patients are made monotone by
#'   isotonic regression, and the dose whose estimate is closest to the
#'   target is selected. This is the selection of the 'UnifiedDoseFinding'
#'   and 'TITEgBOIN' packages. With whole numbers of DLTs as
#'   \code{score_sum} the result is that of \code{\link{boin_select_mtd}}.
#'
#' @references
#'   Mu, R., Yuan, Y., Xu, J., Mandrekar, S. J. and Yin, J. (2019). gBOIN: a
#'   unified model-assisted phase I trial design accounting for toxicity
#'   grades, and binary or continuous end points. Journal of the Royal
#'   Statistical Society: Series C, 68, 289-308. \doi{10.1111/rssc.12263}
#'
#' @examples
#' # Equivalent toxicity scores divided by their maximum of 1.5
#' n_pts <- c(3, 6, 12, 9, 0, 0)
#' score_sum <- c(0.5, 1.5, 5.5, 6.0, 0, 0) / 1.5
#' gboin_select_mtd(n_pts, score_sum, target = 0.47 / 1.5)
#'
#' @seealso \code{\link{sim_gboin}}, \code{\link{boin_select_mtd}}
#'
#' @export
gboin_select_mtd <- function(n_pts, score_sum, target, cutoff_eli = 0.95,
                             extrasafe = FALSE, offset = 0.05,
                             min_mtd_sample = 1) {

  n_pts <- as_count_matrix(n_pts, "n_pts")
  if (is.null(dim(score_sum))) score_sum <- matrix(score_sum, nrow = 1L)
  if (!is.numeric(score_sum) || any(!is.finite(score_sum)) ||
      any(score_sum < 0)) {
    stop("'score_sum' must contain non-negative finite numbers",
         call. = FALSE)
  }
  if (!identical(dim(n_pts), dim(score_sum))) {
    stop("'n_pts' and 'score_sum' must have the same dimensions",
         call. = FALSE)
  }
  if (any(score_sum > n_pts + 1e-9)) {
    stop("'score_sum' must not exceed 'n_pts' at any dose; divide the ",
         "scores by their maximum", call. = FALSE)
  }
  check_scalar_prob(target, "target")
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_flag(extrasafe, "extrasafe")
  check_count(min_mtd_sample, "min_mtd_sample", 1L)

  score_sum <- unname(score_sum)
  storage.mode(score_sum) <- "double"
  est_n <- unname(n_pts)
  storage.mode(est_n) <- "double"

  res <- select_mtd_estimate_cpp(
    eliminated = mtd_eliminated(n_pts, score_sum, target, cutoff_eli,
                                extrasafe, offset),
    n_pts = n_pts,
    est_y = score_sum,
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
