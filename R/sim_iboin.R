#' Operating Characteristics of the iBOIN Design
#'
#' @description
#'   Simulate trials of the iBOIN design of Zhou et al. (2021), which brings
#'   prior information on the DLT probability of every dose into the BOIN
#'   design, and summarize how often each dose is selected as the MTD and how
#'   the patients are distributed over the doses.
#'
#' @inheritParams iboin_simulate
#'
#' @param select_with_prior
#'   Logical scalar. Add the prior information to the estimates of the DLT
#'   probabilities from which the MTD is selected, as
#'   \code{\link{iboin_select_mtd}} does. Defaults to \code{TRUE}, which
#'   reproduces the operating characteristics of the article. With
#'   \code{FALSE} the MTD is selected from the data alone, by
#'   \code{\link{boin_select_mtd}}.
#'
#' @param min_mtd_sample
#'   Integer scalar. Smallest number of patients a dose must have been given
#'   to be selected as the MTD. Defaults to 1.
#'
#' @param overdose_cutoff
#'   Numeric scalar or \code{NULL}. Doses with a true DLT probability above
#'   this value count as overdoses in the summary. Defaults to \code{NULL},
#'   the target.
#'
#' @param keep_trials
#'   Logical scalar. Keep the trial level data in the result. Defaults to
#'   \code{FALSE}.
#'
#' @param verbose
#'   Logical scalar. Report the progress in messages. Defaults to
#'   \code{FALSE}.
#'
#' @return
#'   An object of class \code{c("iboin_oc", "boin_oc")}, a list with the
#'   components of the result of \code{\link{sim_boin}} and in addition
#'   \item{design}{The name of the design, \code{"iBOIN"}.}
#'   \item{skeleton}{The skeleton.}
#'   \item{prior_ess}{The prior effective sample sizes used, after the robust
#'     prior has been applied.}
#'   \item{select_with_prior}{Whether the selection of the MTD used the
#'     prior information.}
#'   \item{boundary}{The \code{iboin_boundary} object used.}
#'   The components \code{lambda_e} and \code{lambda_d} are the boundaries of
#'   the BOIN design, which the iBOIN boundaries approach as patients are
#'   treated at a dose.
#'
#' @details
#'   The trials are simulated by \code{\link{iboin_simulate}} with the
#'   boundaries of \code{\link{iboin_boundary}}. The article gives two ways
#'   to select the MTD at the end of the trial, with and without the prior
#'   information in the estimates of the DLT probabilities, and does not say
#'   which one its simulation study uses. Its operating characteristics are
#'   reproduced with the prior information (see
#'   \code{inst/reproduce/reproduce_Zhou2021.R}), which is therefore the
#'   default. The selection from the data alone
#'   (\code{select_with_prior = FALSE}) keeps a misspecified prior from
#'   biasing the selection, at the price of a lower probability of selecting
#'   the true MTD when the prior is right.
#'
#'   With a prior effective sample size of 0 at every dose the result is,
#'   for the same seed, that of \code{\link{sim_boin}}.
#'
#'   The result can be passed to \code{\link{oc_results}} and
#'   \code{\link{oc_compare}}, where it appears under the name iBOIN.
#'
#' @references
#'   Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
#'   Incorporating historical information to improve phase I clinical trials.
#'   Pharmaceutical Statistics, 20, 1017-1034. \doi{10.1002/pst.2121}
#'
#' @examples
#' # Scenario 3 of Zhou et al. (2021): the prior MTD is the true MTD
#' oc <- sim_iboin(
#'   target = 0.30,
#'   p_true = c(0.08, 0.15, 0.31, 0.45, 0.55),
#'   skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
#'   prior_ess = 3,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_earlystop = 100,
#'   n_trials = 500,
#'   seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{iboin_boundary}}, \code{\link{iboin_simulate}},
#'   \code{\link{sim_iboin_multi}}, \code{\link{sim_boin}}
#'
#' @export
sim_iboin <- function(target, p_true, skeleton, prior_ess, n_cohort,
                      cohort_size, n_trials = 10000, start_dose = 1,
                      n_earlystop = 18, p_saf = NULL, p_tox = NULL,
                      cutoff_eli = 0.95, extrasafe = FALSE, offset = 0.05,
                      robust = FALSE, select_with_prior = TRUE,
                      min_mtd_sample = 1, overdose_cutoff = NULL,
                      n_earlystop_rule = c("with_stay", "simple"),
                      keep_trials = FALSE, verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  n_earlystop_rule <- match.arg(n_earlystop_rule)

  if (verbose) message("Simulating ", n_trials, " trials ...")

  trials <- iboin_simulate(
    target = target, p_true = p_true, skeleton = skeleton,
    prior_ess = prior_ess, n_cohort = n_cohort, cohort_size = cohort_size,
    n_trials = n_trials, start_dose = start_dose, n_earlystop = n_earlystop,
    p_saf = p_saf, p_tox = p_tox, cutoff_eli = cutoff_eli,
    extrasafe = extrasafe, offset = offset, robust = robust,
    n_earlystop_rule = n_earlystop_rule, seed = seed
  )

  if (verbose) message("Selecting the MTD ...")

  check_flag(select_with_prior, "select_with_prior")
  selection <- iboin_select_mtd(
    n_pts = trials$n_pts, n_tox = trials$n_tox, target = target,
    skeleton = skeleton, prior_ess = prior_ess,
    with_prior = select_with_prior, cutoff_eli = cutoff_eli,
    extrasafe = extrasafe, offset = offset, robust = robust,
    min_mtd_sample = min_mtd_sample
  )

  out <- mtd_oc_summary(trials, selection, p_true, target,
                        overdose_cutoff = overdose_cutoff,
                        keep_trials = keep_trials)

  if (verbose) message("Done.")

  out$lambda_e <- trials$boundary$lambda_e_boin
  out$lambda_d <- trials$boundary$lambda_d_boin
  out$design <- "iBOIN"
  out$select_with_prior <- select_with_prior
  out$skeleton <- skeleton
  out$prior_ess <- trials$boundary$prior_ess
  out$boundary <- trials$boundary
  out$call <- call_expr

  structure(out, class = c("iboin_oc", "boin_oc"))
}
