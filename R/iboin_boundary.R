#' Decision Boundaries of the iBOIN Design
#'
#' @description
#'   Tabulate the escalation and de-escalation boundaries of the iBOIN design
#'   of Zhou et al. (2021), which brings prior information on the DLT
#'   probability of every dose into the BOIN design. The prior information is
#'   a skeleton, the prior DLT probability of each dose, and a prior effective
#'   sample size (PESS), the number of patients the prior is worth. The
#'   boundaries depend on the dose and on the number of patients treated
#'   there.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param skeleton
#'   Numeric vector. Prior DLT probability at each dose, strictly between 0
#'   and 1.
#'
#' @param prior_ess
#'   Whole number or vector of whole numbers, one per dose. Prior effective
#'   sample size. A value of 0 gives the boundaries of the BOIN design at
#'   that dose.
#'
#' @param max_n
#'   Integer scalar. Largest number of patients that can be treated at a
#'   single dose, which is the number of columns of the resulting tables.
#'
#' @param p_saf
#'   Numeric scalar. Highest DLT probability deemed subtherapeutic.
#'   Defaults to \code{0.6 * target}.
#'
#' @param p_tox
#'   Numeric scalar. Lowest DLT probability deemed overly toxic.
#'   Defaults to \code{1.4 * target}.
#'
#' @param cutoff_eli
#'   Numeric scalar. Posterior probability cutoff above which a dose is
#'   eliminated for toxicity. Defaults to 0.95.
#'
#' @param extrasafe
#'   Logical scalar. When \code{TRUE}, also tabulate the safety stopping
#'   boundary applied at the lowest dose. Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is
#'   relaxed for the safety stopping rule. Defaults to 0.05.
#'
#' @param robust
#'   Logical scalar. Use the robust prior of the article: when the prior MTD,
#'   the dose whose skeleton value is closest to \code{target}, has a number
#'   of at least half the number of doses, the doses above it get a prior
#'   effective sample size of 0. Defaults to \code{FALSE}.
#'
#' @return
#'   An object of class \code{iboin_boundary}, which is a list with components
#'   \item{prior}{Matrix with one row per dose: the prior probabilities that
#'     the DLT probability equals \code{target}, \code{p_saf} and \code{p_tox}.}
#'   \item{lambda_e, lambda_d}{Matrices, doses by \code{1:max_n} patients, of
#'     the escalation and de-escalation boundaries of the DLT rate.}
#'   \item{n}{Number of patients, 1 to \code{max_n}.}
#'   \item{b_esc}{Integer matrix. Escalate when the number of DLTs is at most
#'     this value; -1 when no number of DLTs escalates.}
#'   \item{b_deesc}{Integer matrix. De-escalate when the number of DLTs
#'     reaches this value; a value above the number of patients means that no
#'     number of DLTs de-escalates.}
#'   \item{b_elim}{Eliminate the dose when the number of DLTs reaches this
#'     value, \code{NA} when no elimination is possible at that sample size.}
#'   \item{b_stop}{Safety stopping boundary at the lowest dose, all \code{NA}
#'     unless \code{extrasafe} is \code{TRUE}.}
#'   \item{searched}{Logical matrix marking the boundaries found by search
#'     rather than by the closed form, see the details.}
#'   \item{prior_ess}{The prior effective sample sizes used, after the
#'     robust prior has been applied.}
#'   \item{prior_mtd}{The dose whose skeleton value is closest to
#'     \code{target}.}
#'   together with the design parameters used.
#'
#' @details
#'   The BOIN design chooses its boundaries to minimize the probability of an
#'   incorrect decision under three point hypotheses: the DLT probability of
#'   the current dose equals \code{target} (stay), \code{p_saf} (escalate) or
#'   \code{p_tox} (de-escalate). With equal prior probabilities of the three
#'   hypotheses the boundaries are those of \code{\link{boin_lambda}}. The
#'   iBOIN design derives the prior probabilities of the hypotheses at each
#'   dose from the skeleton and the PESS (equation 5 of the article): the
#'   posterior probabilities of the hypotheses given \code{x} DLTs among
#'   \code{prior_ess} patients, averaged over the binomial distribution of
#'   \code{x} with the skeleton value as DLT probability. The boundaries then
#'   follow from equation 3 of the article, in which the prior probabilities
#'   shift the BOIN boundaries by an amount that decreases as one over the
#'   number of patients treated at the dose, so that the data override the
#'   prior as they accumulate.
#'
#'   The design escalates when the observed DLT rate is at most
#'   \code{lambda_e} and de-escalates when it exceeds \code{lambda_d}; the
#'   integer boundaries are cut as in \code{\link{boin_boundary}}. With the
#'   skeleton (0.10, 0.19, 0.30, 0.42, 0.54), a PESS of 3 at every dose and a
#'   target of 0.30 they reproduce Table 1 of the article.
#'
#'   Where the closed form gives \code{lambda_e > lambda_d}, the article
#'   determines the boundaries by a numerical search. Here the search is over
#'   the integer boundaries \code{e < d} themselves: the pair that minimizes
#'   the probability of an incorrect decision is used, and among pairs with
#'   the same probability the smallest \code{e} and then the smallest
#'   \code{d}. This happens mostly for a single patient at a dose.
#'
#'   The elimination boundary and the safety stopping boundary are those of
#'   the BOIN design: they use a uniform prior, not the skeleton, so that a
#'   misspecified prior cannot weaken the safety rules. As in the BOIN
#'   design the de-escalation boundary is capped at the elimination
#'   boundary. The escalation boundary is kept below the de-escalation
#'   boundary, so that a dose that has just been eliminated is never
#'   followed by an escalation.
#'
#'   The PESS must be a whole number, because the prior probabilities sum
#'   over the possible numbers of DLTs among \code{prior_ess} patients. The
#'   article recommends a PESS between one third and one half of the number
#'   of patients per dose, for example 2 or 3 for 30 patients and 5 doses.
#'
#' @references
#'   Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
#'   Incorporating historical information to improve phase I clinical trials.
#'   Pharmaceutical Statistics, 20, 1017-1034. \doi{10.1002/pst.2121}
#'
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I
#'   Clinical Trials. Journal of the Royal Statistical Society: Series C, 64,
#'   507-523.
#'
#' @examples
#' # Table 1 of Zhou et al. (2021)
#' bd <- iboin_boundary(
#'   target = 0.30, skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
#'   prior_ess = 3, max_n = 30
#' )
#' print(bd, cohort_size = 3)
#'
#' # Without prior information the boundaries are those of the BOIN design
#' iboin_boundary(0.30, c(0.10, 0.19, 0.30, 0.42, 0.54), 0, 9)$b_esc
#' boin_boundary(0.30, 9)$b_esc
#'
#' @seealso \code{\link{sim_iboin}}, \code{\link{boin_boundary}}
#'
#' @export
iboin_boundary <- function(target, skeleton, prior_ess, max_n, p_saf = NULL,
                           p_tox = NULL, cutoff_eli = 0.95,
                           extrasafe = FALSE, offset = 0.05, robust = FALSE) {

  base <- boin_boundary(target, max_n, p_saf = p_saf, p_tox = p_tox,
                        cutoff_eli = cutoff_eli, extrasafe = extrasafe,
                        offset = offset)
  p_saf <- base$p_saf
  p_tox <- base$p_tox
  prior <- iboin_prior(target, skeleton, prior_ess, p_saf, p_tox,
                       robust = robust)

  n_doses <- length(skeleton)
  n <- base$n
  max_n <- length(n)
  dose_names <- paste0("DL", seq_len(n_doses))
  shape <- list(dose_names, as.character(n))

  den_e <- log(target * (1 - p_saf) / (p_saf * (1 - target)))
  den_d <- log(p_tox * (1 - target) / (target * (1 - p_tox)))

  lambda_e <- matrix(NA_real_, n_doses, max_n, dimnames = shape)
  lambda_d <- matrix(NA_real_, n_doses, max_n, dimnames = shape)
  b_esc <- matrix(NA_integer_, n_doses, max_n, dimnames = shape)
  b_deesc <- matrix(NA_integer_, n_doses, max_n, dimnames = shape)
  searched <- matrix(FALSE, n_doses, max_n, dimnames = shape)

  for (j in seq_len(n_doses)) {
    pj <- unname(prior$prob[j, ])
    le <- base$lambda_e + log(pj[2L] / pj[1L]) / (n * den_e)
    ld <- base$lambda_d + log(pj[1L] / pj[3L]) / (n * den_d)

    e <- floor(le * n)
    scaled_d <- ld * n
    is_whole <- abs(round(scaled_d) - scaled_d) < 1e-12
    d <- ifelse(is_whole, round(scaled_d) + 1, ceiling(scaled_d))

    # The closed form is not valid where it puts the escalation boundary
    # above the de-escalation boundary.
    for (k in which(le > ld)) {
      found <- iboin_search(n[k], pj, target, p_saf, p_tox)
      e[k] <- found[1L]
      d[k] <- found[2L]
      searched[j, k] <- TRUE
    }

    e <- pmin(pmax(e, -1), n)
    d <- pmin(pmax(d, 0), n + 1)
    capped <- !is.na(base$b_elim) & d > base$b_elim
    d[capped] <- base$b_elim[capped]
    e <- pmin(e, d - 1)

    lambda_e[j, ] <- le
    lambda_d[j, ] <- ld
    b_esc[j, ] <- as.integer(e)
    b_deesc[j, ] <- as.integer(d)
  }

  structure(
    list(
      prior = prior$prob,
      lambda_e = lambda_e,
      lambda_d = lambda_d,
      n = n,
      b_esc = b_esc,
      b_deesc = b_deesc,
      b_elim = base$b_elim,
      b_stop = base$b_stop,
      searched = searched,
      skeleton = skeleton,
      prior_ess = prior$prior_ess,
      prior_ess_input = prior$prior_ess_input,
      prior_mtd = prior$prior_mtd,
      robust = robust,
      target = target,
      p_saf = p_saf,
      p_tox = p_tox,
      lambda_e_boin = base$lambda_e,
      lambda_d_boin = base$lambda_d,
      cutoff_eli = cutoff_eli,
      extrasafe = extrasafe,
      offset = offset
    ),
    class = "iboin_boundary"
  )
}
