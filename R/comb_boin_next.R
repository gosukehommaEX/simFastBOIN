#' Next Dose Combination in a BOIN Drug-Combination Trial
#'
#' @description
#'   Decide the dose combination for the next cohort of a trial run with the
#'   BOIN drug-combination design of Lin and Yin (2017), from the numbers of
#'   patients and DLTs observed so far at every combination and the combination
#'   that treated the last cohort. The rules are those of
#'   \code{\link{comb_boin_simulate}}, so the decisions taken during the trial
#'   are the ones whose operating characteristics \code{\link{sim_comb_boin}}
#'   reports.
#'
#' @param target
#'   Numeric scalar. Target DLT probability, for example 0.30.
#'
#' @param n_pts
#'   Integer matrix. Number of patients treated so far at each dose combination,
#'   with one row per dose level of drug A and one column per dose level of
#'   drug B.
#'
#' @param n_tox
#'   Integer matrix of the same shape as \code{n_pts}. Number of patients with a
#'   DLT at each combination.
#'
#' @param current
#'   Integer vector of length two. Dose levels of drug A and drug B of the
#'   combination that treated the last cohort.
#'
#' @param n_earlystop
#'   Integer scalar. The trial stops once this many patients have been treated at
#'   the current combination and the design would neither escalate nor
#'   de-escalate, or cannot. Defaults to 18. Use the value of the simulations
#'   that support the protocol.
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
#'   Numeric scalar. Posterior probability cutoff for eliminating a combination.
#'   Defaults to 0.95.
#'
#' @param extrasafe
#'   Logical scalar. Apply the stricter safety stopping rule at the lowest
#'   combination. Defaults to \code{FALSE}.
#'
#' @param offset
#'   Numeric scalar between 0 and 0.5. Amount by which \code{cutoff_eli} is
#'   relaxed for the safety stopping rule. Defaults to 0.05.
#'
#' @return
#'   An object of class \code{comb_boin_next}, which is a list with components
#'   \item{decision}{One of \code{"escalate"}, \code{"de-escalate"}, \code{"stay"} and \code{"stop"}.}
#'   \item{next_dose}{Integer matrix with columns \code{dose_a} and \code{dose_b}: one row for the next combination, two rows when two candidates tie, and no row when the trial stops.}
#'   \item{tie}{Whether two candidates tie, in which case the design picks one of them at random.}
#'   \item{stop_reason}{\code{NA}, or why the trial stops: \code{"lowest_dose_eliminated"}, \code{"lowest_dose_too_toxic"} or \code{"n_earlystop"}.}
#'   \item{note}{A sentence explaining a decision to stay when the rule called for a move that is not available, otherwise \code{NA}.}
#'   \item{candidates}{Data frame of the two combinations considered for an escalation or a de-escalation, one level of drug A and one level of drug B away, with their data, the posterior probability that their DLT rate lies between \code{lambda_e} and \code{lambda_d}, the score used to choose between them and whether they are available; \code{NULL} when no move was considered.}
#'   \item{eliminated}{Logical matrix of the combinations eliminated by the data.}
#'   \item{boundary}{The escalation, de-escalation and elimination boundaries at the number of patients treated at the current combination.}
#'   together with the data, the current combination and the design parameters.
#'
#' @details
#'   The decision follows the same steps as a simulated trial after each
#'   cohort. A combination whose DLT count reaches its elimination boundary is
#'   eliminated together with every combination with at least the same level of
#'   both drugs, and the trial stops if the lowest combination is eliminated or,
#'   with \code{extrasafe}, judged too toxic. The trial also stops when
#'   \code{n_earlystop} patients have been treated at the current combination
#'   and the design would neither escalate nor de-escalate. Otherwise the
#'   observed DLT rate at the current combination is compared with the
#'   boundaries of \code{\link{boin_boundary}}, and an escalation or a
#'   de-escalation moves one level of drug A or of drug B to the candidate with
#'   the higher score: the posterior probability, under a Beta(0.5, 0.5) prior,
#'   that its DLT rate lies between \code{lambda_e} and \code{lambda_d}, plus
#'   0.0005 per patient already treated there. An eliminated combination is not
#'   a candidate for escalation. When the two scores are equal the design
#'   chooses at random, and the function returns both.
#'
#'   During a simulated trial a combination is eliminated as soon as its data
#'   reach the boundary and is never treated again, so the combinations
#'   eliminated by the data are exactly those eliminated in the course of the
#'   trial.
#'
#'   \code{BOIN::next.comb()} applies somewhat different rules from the
#'   simulations of \code{BOIN::get.oc.comb()}: it stops as soon as
#'   \code{n.earlystop} patients have been treated at the current combination,
#'   it does not escalate drug A when a combination with the next level of drug
#'   A and at most the current level of drug B has an observed DLT rate at or
#'   above \code{lambda_d} (and likewise for drug B), and it draws a random
#'   number to break ties. This function follows the simulations instead, so
#'   that a trial is run with the rules whose operating characteristics were
#'   evaluated. Where the rules coincide, the two functions agree.
#'
#' @references
#'   Lin, R. and Yin, G. (2017). Bayesian Optimal Interval Design for Dose
#'   Finding in Drug-Combination Trials. Statistical Methods in Medical
#'   Research, 26(5), 2155-2167.
#'
#' @examples
#' # Three levels of drug A and four of drug B; the last cohort was treated at
#' # A2B2, where two of three patients had a DLT
#' n_pts <- matrix(c(6, 6, 0, 0,
#'                   3, 3, 0, 0,
#'                   0, 0, 0, 0), nrow = 3, byrow = TRUE)
#' n_tox <- matrix(c(0, 1, 0, 0,
#'                   1, 2, 0, 0,
#'                   0, 0, 0, 0), nrow = 3, byrow = TRUE)
#'
#' nxt <- comb_boin_next(target = 0.25, n_pts = n_pts, n_tox = n_tox,
#'                       current = c(2, 2))
#' nxt
#' nxt$next_dose
#'
#' @seealso \code{\link{plot.comb_boin_next}}, \code{\link{boin_decision_table}},
#'   \code{\link{sim_comb_boin}}
#'
#' @export
comb_boin_next <- function(target, n_pts, n_tox, current, n_earlystop = 18,
                           p_saf = NULL, p_tox = NULL, cutoff_eli = 0.95,
                           extrasafe = FALSE, offset = 0.05) {

  if (!is.matrix(n_pts) || !is.matrix(n_tox)) {
    stop("'n_pts' and 'n_tox' must be matrices, with one row per level of ",
         "drug A and one column per level of drug B", call. = FALSE)
  }
  if (!identical(dim(n_pts), dim(n_tox))) {
    stop("'n_pts' and 'n_tox' must have the same dimensions", call. = FALSE)
  }
  n_rows <- nrow(n_pts)
  n_cols <- ncol(n_pts)
  if (n_rows * n_cols < 2L) {
    stop("'n_pts' must have at least two combinations", call. = FALSE)
  }
  n_pts <- matrix(as_comb_count_array(n_pts, "n_pts"), n_rows, n_cols)
  n_tox <- matrix(as_comb_count_array(n_tox, "n_tox"), n_rows, n_cols)
  if (any(n_tox > n_pts)) {
    stop("'n_tox' must not exceed 'n_pts' at any combination", call. = FALSE)
  }
  if (!is.numeric(current) || length(current) != 2L ||
      any(!is.finite(current)) || any(current != round(current)) ||
      any(current < 1) || current[1L] > n_rows || current[2L] > n_cols) {
    stop("'current' must give the dose levels of drug A and drug B of a ",
         "combination in 'n_pts'", call. = FALSE)
  }
  a <- as.integer(current[1L])
  b <- as.integer(current[2L])
  if (n_pts[a, b] == 0L) {
    stop("no patient has been treated at the 'current' combination",
         call. = FALSE)
  }
  check_count(n_earlystop, "n_earlystop", 1L)
  check_flag(extrasafe, "extrasafe")

  bound <- boin_boundary(target, max(n_pts), p_saf = p_saf, p_tox = p_tox,
                         cutoff_eli = cutoff_eli, extrasafe = extrasafe,
                         offset = offset)
  b_elim <- bound$b_elim
  b_elim[is.na(b_elim)] <- 0L

  res <- comb_boin_next_cpp(
    n_pts = n_pts, n_tox = n_tox, current_row = a, current_col = b,
    n_earlystop = as.integer(n_earlystop), extrasafe = extrasafe,
    target = as.numeric(target), cutoff_eli = as.numeric(cutoff_eli),
    offset = as.numeric(offset), lambda_e = bound$lambda_e,
    lambda_d = bound$lambda_d, b_esc = bound$b_esc, b_deesc = bound$b_deesc,
    b_elim = as.integer(b_elim)
  )

  stop_reasons <- c("lowest_dose_eliminated", "lowest_dose_too_toxic",
                    "n_earlystop")
  decision <- "stay"
  next_dose <- matrix(c(a, b), nrow = 1L)
  tie <- FALSE
  stop_reason <- NA_character_
  note <- NA_character_
  candidates <- NULL

  if (res$stop_code >= 0L) {
    decision <- "stop"
    stop_reason <- stop_reasons[res$stop_code + 1L]
    next_dose <- matrix(integer(0), ncol = 2L)
  } else if (res$action != 0L) {
    # 1 is an escalation and 2 a de-escalation: the candidates are one level
    # of drug A and one level of drug B away, in that order.
    step <- if (res$action == 1L) 1L else -1L
    cand_a <- c(a + step, a)
    cand_b <- c(b, b + step)
    inside <- cand_a >= 1L & cand_a <= n_rows & cand_b >= 1L & cand_b <= n_cols
    available <- res$score > 0
    cand_n <- rep(NA_integer_, 2L)
    cand_y <- rep(NA_integer_, 2L)
    prob <- rep(NA_real_, 2L)
    for (k in which(inside)) {
      cand_n[k] <- n_pts[cand_a[k], cand_b[k]]
      cand_y[k] <- n_tox[cand_a[k], cand_b[k]]
      shape1 <- cand_y[k] + 0.5
      shape2 <- (cand_n[k] - cand_y[k]) + 0.5
      prob[k] <- stats::pbeta(bound$lambda_d, shape1, shape2) -
        stats::pbeta(bound$lambda_e, shape1, shape2)
    }
    candidates <- data.frame(
      direction = c("drug A", "drug B"),
      dose_a = ifelse(inside, cand_a, NA_integer_),
      dose_b = ifelse(inside, cand_b, NA_integer_),
      n_pts = cand_n,
      n_tox = cand_y,
      prob_in_interval = prob,
      score = ifelse(available, res$score, NA_real_),
      available = available,
      stringsAsFactors = FALSE
    )
    if (!any(available)) {
      note <- if (res$action == 1L) {
        paste("The rule calls for escalation, but no combination one level",
              "up of either drug is available, so the dose stays.")
      } else {
        paste("The rule calls for de-escalation, but the current combination",
              "is the lowest, so the dose stays.")
      }
    } else {
      decision <- if (res$action == 1L) "escalate" else "de-escalate"
      best <- which(available & res$score == max(res$score))
      tie <- length(best) > 1L
      next_dose <- cbind(cand_a[best], cand_b[best])
    }
  }
  storage.mode(next_dose) <- "integer"
  dimnames(next_dose) <- list(NULL, c("dose_a", "dose_b"))

  labels <- list(paste0("A", seq_len(n_rows)), paste0("B", seq_len(n_cols)))
  eliminated <- res$eliminated
  dimnames(eliminated) <- labels
  dimnames(n_pts) <- labels
  dimnames(n_tox) <- labels

  n_current <- n_pts[a, b]
  structure(
    list(
      decision = decision,
      next_dose = next_dose,
      tie = tie,
      stop_reason = stop_reason,
      note = note,
      candidates = candidates,
      eliminated = eliminated,
      boundary = list(
        n = n_current,
        escalate_at_most = bound$b_esc[n_current],
        deescalate_at_least = bound$b_deesc[n_current],
        eliminate_at_least = bound$b_elim[n_current]
      ),
      current = c(dose_a = a, dose_b = b),
      n_pts = n_pts,
      n_tox = n_tox,
      target = target,
      lambda_e = bound$lambda_e,
      lambda_d = bound$lambda_d,
      settings = list(
        n_earlystop = n_earlystop,
        p_saf = bound$p_saf,
        p_tox = bound$p_tox,
        cutoff_eli = cutoff_eli,
        extrasafe = extrasafe,
        offset = offset
      )
    ),
    class = "comb_boin_next"
  )
}
