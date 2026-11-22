#' Prior Probabilities of the Three Hypotheses of the iBOIN Design
#'
#' @description
#'   Internal helper. Check the skeleton and the prior effective sample
#'   sizes, apply the robust prior when requested, and turn them into the
#'   prior probabilities of the three point hypotheses at every dose, as in
#'   equation (5) of Zhou et al. (2021).
#'
#' @param target Target DLT probability.
#' @param skeleton Prior DLT probability at each dose.
#' @param prior_ess Prior effective sample size, one value or one per dose.
#' @param p_saf Highest DLT probability deemed subtherapeutic.
#' @param p_tox Lowest DLT probability deemed overly toxic.
#' @param robust Whether to drop the prior information above the prior MTD
#'   when the prior MTD lies in the upper half of the doses.
#'
#' @return A list with the matrix \code{prob} of prior probabilities (one row
#'   per dose, columns \code{target}, \code{p_saf} and \code{p_tox}), the
#'   effective sample sizes used (\code{prior_ess}) and given
#'   (\code{prior_ess_input}), and the prior MTD (\code{prior_mtd}).
#'
#' @importFrom stats dbinom
#'
#' @noRd
iboin_prior <- function(target, skeleton, prior_ess, p_saf, p_tox,
                        robust = FALSE) {

  if (!is.numeric(skeleton) || length(skeleton) < 1L ||
      any(!is.finite(skeleton)) || any(skeleton <= 0) || any(skeleton >= 1)) {
    stop("'skeleton' must be a numeric vector of prior DLT probabilities ",
         "strictly between 0 and 1", call. = FALSE)
  }
  n_doses <- length(skeleton)
  if (!is.numeric(prior_ess) || !(length(prior_ess) %in% c(1L, n_doses)) ||
      any(!is.finite(prior_ess)) || any(prior_ess < 0) ||
      any(prior_ess != round(prior_ess))) {
    stop("'prior_ess' must be one non-negative whole number or one for ",
         "each dose of 'skeleton'", call. = FALSE)
  }
  check_flag(robust, "robust")

  ess_input <- as.integer(rep_len(prior_ess, n_doses))
  ess <- ess_input
  prior_mtd <- which.min(abs(skeleton - target))
  if (robust && prior_mtd >= n_doses / 2 && prior_mtd < n_doses) {
    ess[(prior_mtd + 1L):n_doses] <- 0L
  }

  phi <- c(target, p_saf, p_tox)
  prob <- matrix(NA_real_, n_doses, 3L,
                 dimnames = list(paste0("DL", seq_len(n_doses)),
                                 c("target", "p_saf", "p_tox")))
  for (j in seq_len(n_doses)) {
    x <- 0:ess[j]
    lik <- outer(x, phi, function(x, p) p^x * (1 - p)^(ess[j] - x))
    post <- lik / rowSums(lik)
    prob[j, ] <- colSums(post * dbinom(x, ess[j], skeleton[j]))
  }

  list(prob = prob, prior_ess = ess, prior_ess_input = ess_input,
       prior_mtd = prior_mtd)
}
