#' Integer Decision Boundaries for the BOIN Design
#'
#' @description
#'   Tabulate, for every attainable number of patients treated at a dose, the
#'   number of DLTs that triggers escalation, de-escalation, elimination and, when
#'   requested, the stricter safety stopping rule at the lowest dose.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param max_n
#'   Integer scalar. Largest number of patients that can be treated at a single
#'   dose, which is the number of columns of the resulting table. For a trial with
#'   \code{n_cohort} cohorts of size \code{cohort_size} this is
#'   \code{n_cohort * cohort_size}.
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
#'   Numeric scalar. Posterior probability cutoff above which a dose is eliminated
#'   for toxicity. Defaults to 0.95.
#'
#' @param extrasafe
#'   Logical scalar. When \code{TRUE}, also tabulate the safety stopping boundary
#'   applied at the lowest dose. Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is relaxed
#'   for the safety stopping rule. Defaults to 0.05.
#'
#' @param stay_on_1_of_3
#'   Logical scalar. When \code{TRUE}, one DLT out of three patients leads to
#'   staying at the current dose rather than de-escalating. Defaults to
#'   \code{FALSE}. See the details.
#'
#' @return
#'   An object of class \code{boin_boundary}, which is a list with components
#'   \item{lambda_e}{Escalation interval boundary.}
#'   \item{lambda_d}{De-escalation interval boundary.}
#'   \item{n}{Number of patients, 1 to \code{max_n}.}
#'   \item{b_esc}{Escalate when the number of DLTs is at most this value.}
#'   \item{b_deesc}{De-escalate when the number of DLTs reaches this value.}
#'   \item{b_elim}{Eliminate the dose when the number of DLTs reaches this value,
#'     \code{NA} when no elimination is possible at that sample size.}
#'   \item{b_stop}{Safety stopping boundary at the lowest dose, all \code{NA}
#'     unless \code{extrasafe} is \code{TRUE}.}
#'   \item{stay_on_1_of_3}{Whether the modification was requested.}
#'   \item{stay_on_1_of_3_applied}{Whether it changed the boundaries.}
#'   together with the design parameters used.
#'
#' @details
#'   The escalation boundary is \code{floor(lambda_e * n)} and the de-escalation
#'   boundary is \code{ceiling(lambda_d * n)}, with the convention that when
#'   \code{lambda_d * n} is a whole number de-escalation requires one more DLT. The
#'   elimination boundary is the smallest positive number of DLTs for which
#'   \code{Pr(p > target | data)} exceeds \code{cutoff_eli} under a uniform
#'   Beta(1, 1) prior, and is not evaluated before three patients have been
#'   treated. The de-escalation boundary is capped at the elimination boundary.
#'
#'   With the default thresholds these definitions reproduce
#'   \code{BOIN::get.boundary()} exactly.
#'
#' @section Staying on one DLT out of three:
#'   For some target rates the optimal BOIN decision after one DLT in three
#'   patients is to de-escalate, because the likelihood of overdosing then exceeds
#'   the likelihood of proper dosing. Long practice with the 3+3 design has
#'   nonetheless made staying at the dose the widely accepted choice, and
#'   \code{stay_on_1_of_3 = TRUE} aligns the design with that practice by raising
#'   the de-escalation boundary at three patients from one DLT to two. Nothing
#'   else in the table changes.
#'
#'   The modification is applied only where it is meaningful, that is when one DLT
#'   out of three currently triggers de-escalation. It never overrides an
#'   escalation, and it never overrides an elimination, which is a safety rule.
#'   Whether it took effect is reported in \code{stay_on_1_of_3_applied}. With the
#'   default thresholds it takes effect for target rates from about 0.098 to
#'   0.279; below that range one DLT out of three already eliminates the dose, and
#'   above it the design already stays.
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' bd <- boin_boundary(target = 0.30, max_n = 18)
#' bd
#'
#' # At a target of 0.25 one DLT out of three de-escalates by default
#' boin_boundary(target = 0.25, max_n = 18)$b_deesc[3]
#'
#' # and stays once the modification is switched on
#' modified <- boin_boundary(target = 0.25, max_n = 18, stay_on_1_of_3 = TRUE)
#' modified$b_deesc[3]
#' modified$stay_on_1_of_3_applied
#'
#' # At a target of 0.30 the design already stays, so nothing changes
#' boin_boundary(target = 0.30, max_n = 18, stay_on_1_of_3 = TRUE)$stay_on_1_of_3_applied
#'
#' @seealso \code{\link{boin_lambda}}, \code{\link{boin_decision_table}}
#'
#' @importFrom stats pbeta
#'
#' @export
boin_boundary <- function(target, max_n, p_saf = NULL, p_tox = NULL,
                          cutoff_eli = 0.95, extrasafe = FALSE, offset = 0.05,
                          stay_on_1_of_3 = FALSE) {

  if (is.null(p_saf)) p_saf <- 0.6 * target
  if (is.null(p_tox)) p_tox <- 1.4 * target

  check_count(max_n, "max_n", 1L)
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_flag(extrasafe, "extrasafe")
  check_flag(stay_on_1_of_3, "stay_on_1_of_3")
  if (!is.numeric(offset) || length(offset) != 1L || !is.finite(offset) ||
      offset <= 0 || offset >= 0.5) {
    stop("'offset' must be a single number strictly between 0 and 0.5", call. = FALSE)
  }

  lambda <- boin_lambda(target, p_saf = p_saf, p_tox = p_tox)
  max_n <- as.integer(max_n)
  n <- seq_len(max_n)

  b_esc <- as.integer(floor(lambda$lambda_e * n))

  scaled_d <- lambda$lambda_d * n
  is_whole <- abs(round(scaled_d) - scaled_d) < 1e-12
  b_deesc <- as.integer(ifelse(is_whole, round(scaled_d) + 1, ceiling(scaled_d)))

  b_elim <- rep(NA_integer_, max_n)
  b_stop <- rep(NA_integer_, max_n)
  for (k in n) {
    if (k < 3L) next
    n_tox <- seq_len(k)
    prob_above <- 1 - pbeta(target, n_tox + 1, k - n_tox + 1)
    hit_elim <- which(prob_above > cutoff_eli)
    if (length(hit_elim) > 0L) b_elim[k] <- n_tox[hit_elim[1L]]
    if (extrasafe) {
      hit_stop <- which(prob_above > cutoff_eli - offset)
      if (length(hit_stop) > 0L) b_stop[k] <- n_tox[hit_stop[1L]]
    }
  }

  capped <- !is.na(b_elim) & b_deesc > b_elim
  b_deesc[capped] <- b_elim[capped]

  # One DLT out of three: raise the de-escalation boundary so that the decision
  # becomes stay. Only where one DLT currently de-escalates, never over an
  # escalation and never over an elimination.
  applied <- FALSE
  if (stay_on_1_of_3 && max_n >= 3L) {
    de_escalates <- b_deesc[3L] <= 1L
    not_escalation <- b_esc[3L] < 1L
    not_elimination <- is.na(b_elim[3L]) || b_elim[3L] > 1L
    if (de_escalates && not_escalation && not_elimination) {
      b_deesc[3L] <- 2L
      applied <- TRUE
    }
  }

  structure(
    list(
      lambda_e = lambda$lambda_e,
      lambda_d = lambda$lambda_d,
      n = n,
      b_esc = b_esc,
      b_deesc = b_deesc,
      b_elim = b_elim,
      b_stop = b_stop,
      target = target,
      p_saf = p_saf,
      p_tox = p_tox,
      cutoff_eli = cutoff_eli,
      extrasafe = extrasafe,
      offset = offset,
      stay_on_1_of_3 = stay_on_1_of_3,
      stay_on_1_of_3_applied = applied
    ),
    class = "boin_boundary"
  )
}
