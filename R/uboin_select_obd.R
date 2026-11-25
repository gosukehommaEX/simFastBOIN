#' Select the OBD from Completed U-BOIN Trials
#'
#' @description
#'   Choose the optimal biological dose (OBD) of the U-BOIN design for each of
#'   a set of completed trials from the numbers of patients with each outcome
#'   at each dose.
#'
#' @param n_outcome
#'   Integer matrix with one row per dose and four columns, the numbers of
#'   patients with the outcomes (no response, DLT), (no response, no DLT),
#'   (response, DLT) and (response, no DLT) in one trial, or an integer array
#'   with dimensions trial by dose by outcome for several trials, as returned
#'   by \code{\link{uboin_simulate}}.
#'
#' @param tox_upper
#'   Numeric scalar. Highest acceptable DLT probability.
#'
#' @param resp_lower
#'   Numeric scalar. Lowest acceptable response probability.
#'
#' @param eliminated
#'   Logical vector with one element per dose for one trial, a logical matrix
#'   with one row per trial for several, or \code{NULL}. Doses eliminated in
#'   stage I, which cannot be selected. Defaults to \code{NULL}, none.
#'
#' @param utility
#'   Numeric vector of length four. Utility of the four outcomes, in the order
#'   of the columns of \code{n_outcome}. Defaults to \code{c(0, 30, 50, 100)}.
#'
#' @param cutoff_eli
#'   Numeric scalar. Cutoff of the toxicity rule. Defaults to 0.95.
#'
#' @param cutoff_eli_resp
#'   Numeric scalar. Cutoff of the futility rule. Defaults to 0.90.
#'
#' @param prior_tox
#'   Numeric vector of length two. Beta prior of the toxicity rule.
#'   Defaults to \code{c(1, 1)}.
#'
#' @param prior_resp
#'   Numeric vector of length two. Beta prior of the futility rule.
#'   Defaults to \code{c(1, 1)}.
#'
#' @return
#'   A data frame with one row per trial and columns
#'   \item{trial}{Trial index.}
#'   \item{obd}{Selected dose level, or \code{NA} when no dose is admissible.}
#'   \item{reason}{\code{"selected"} or \code{"no_admissible_dose"}.}
#'   \item{utility_DL1, ...}{Posterior mean utility of each dose.}
#'   \item{admissible_DL1, ...}{Whether each dose is admissible: tried, not
#'     eliminated, neither toxic nor futile.}
#'
#' @details
#'   The OBD is the admissible dose with the largest posterior mean utility,
#'   the lowest among ties, with the rules of rule B2 described in
#'   \code{\link{uboin_simulate}}.
#'
#' @references
#'   Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian optimal
#'   interval (U-BOIN) phase I/II design to identify the optimal biological
#'   dose for targeted and immune therapies. Statistics in Medicine, 38(28),
#'   S5299-S5316. \doi{10.1002/sim.8361}
#'
#' @examples
#' # Table S3 of the supporting information of the article, with its utilities
#' n_outcome <- rbind(c(0, 3, 0, 0),
#'                    c(0, 2, 0, 1),
#'                    c(0, 3, 1, 1),
#'                    c(1, 3, 1, 4),
#'                    c(1, 2, 1, 2))
#' uboin_select_obd(n_outcome, tox_upper = 0.30, resp_lower = 0.20,
#'                  utility = c(0, 15, 25, 100))
#'
#' @seealso \code{\link{uboin_simulate}}, \code{\link{sim_uboin}}
#'
#' @export
uboin_select_obd <- function(n_outcome, tox_upper, resp_lower,
                             eliminated = NULL, utility = c(0, 30, 50, 100),
                             cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
                             prior_tox = c(1, 1), prior_resp = c(1, 1)) {

  if (!is.numeric(n_outcome) || any(is.na(n_outcome)) ||
      any(n_outcome < 0) || any(n_outcome != round(n_outcome))) {
    stop("'n_outcome' must hold non-negative whole numbers", call. = FALSE)
  }
  dims <- dim(n_outcome)
  if (length(dims) == 2L) {
    n_outcome <- array(n_outcome, c(1L, dims))
  } else if (length(dims) != 3L) {
    stop("'n_outcome' must be a matrix for one trial, or an array with ",
         "dimensions trial by dose by outcome", call. = FALSE)
  }
  dims <- dim(n_outcome)
  if (dims[3L] != 4L) {
    stop("'n_outcome' must have four outcomes", call. = FALSE)
  }
  n_trials <- dims[1L]
  n_doses <- dims[2L]
  storage.mode(n_outcome) <- "integer"

  if (is.null(eliminated)) {
    eliminated <- matrix(FALSE, n_trials, n_doses)
  }
  if (!is.logical(eliminated) || any(is.na(eliminated))) {
    stop("'eliminated' must be logical without missing values", call. = FALSE)
  }
  if (is.null(dim(eliminated))) {
    eliminated <- matrix(eliminated, nrow = 1L)
  }
  if (!identical(dim(eliminated), c(n_trials, n_doses))) {
    stop("'eliminated' must have one row per trial and one column per dose",
         call. = FALSE)
  }

  des <- uboin_check_design(tox_upper, resp_lower, NULL, NULL, NULL, utility,
                            cutoff_eli, cutoff_eli_resp, prior_tox,
                            prior_resp, stage1 = FALSE)
  max_n <- max(1L, max(apply(n_outcome, c(1L, 2L), sum)))
  tab <- uboin_tables(des, max_n)

  res <- uboin_select_cpp(
    n_outcome = n_outcome, eliminated = eliminated, utility = des$utility,
    prior_a = 0.25, max_n = as.integer(max_n), toxic = tab$toxic,
    futile = tab$futile
  )

  dose_names <- paste0("DL", seq_len(n_doses))
  utility_est <- res$utility
  admissible <- res$admissible
  colnames(utility_est) <- paste0("utility_", dose_names)
  colnames(admissible) <- paste0("admissible_", dose_names)

  out <- data.frame(
    trial = seq_len(n_trials),
    obd = res$obd,
    reason = ifelse(is.na(res$obd), "no_admissible_dose", "selected"),
    stringsAsFactors = FALSE
  )
  cbind(out, as.data.frame(utility_est), as.data.frame(admissible))
}
