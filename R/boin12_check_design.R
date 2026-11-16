#' Check and Complete the Design Parameters of BOIN12
#'
#' @param target,target_resp,p_saf,p_tox
#'   See \code{\link{boin12_simulate}}.
#'
#' @param psi00,psi11,benchmark_weight,utility_prior
#'   See \code{\link{boin12_simulate}}.
#'
#' @return A list with the checked parameters, the boundaries \code{lambda_e}
#'   and \code{lambda_d}, and the lowest acceptable utility \code{u_low} and
#'   the benchmark \code{u_benchmark}, both on the scale 0 to 100.
#'
#' @noRd
boin12_check_design <- function(target, target_resp, p_saf = NULL,
                                p_tox = NULL, psi00 = 40, psi11 = 60,
                                benchmark_weight = 0.5,
                                utility_prior = c(1, 1)) {
  check_scalar_prob(target, "target")
  check_scalar_prob(target_resp, "target_resp")
  if (is.null(p_saf)) p_saf <- 0.6 * target
  if (is.null(p_tox)) p_tox <- 1.4 * target
  lam <- boin_lambda(target, p_saf = p_saf, p_tox = p_tox)
  for (nm in c("psi00", "psi11")) {
    value <- get(nm)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value > 100) {
      stop("'", nm, "' must be a single number from 0 to 100", call. = FALSE)
    }
  }
  if (!is.numeric(benchmark_weight) || length(benchmark_weight) != 1L ||
      !is.finite(benchmark_weight) || benchmark_weight <= 0 ||
      benchmark_weight > 1) {
    stop("'benchmark_weight' must be a single number above 0 and at most 1",
         call. = FALSE)
  }
  if (!is.numeric(utility_prior) || length(utility_prior) != 2L ||
      any(!is.finite(utility_prior)) || any(utility_prior <= 0)) {
    stop("'utility_prior' must hold two positive numbers", call. = FALSE)
  }
  u_low <- 100 * target_resp * (1 - target) +
    psi00 * (1 - target) * (1 - target_resp) + psi11 * target * target_resp
  u_benchmark <- benchmark_weight * u_low + (1 - benchmark_weight) * 100
  list(target = target, target_resp = target_resp, p_saf = p_saf,
       p_tox = p_tox, lambda_e = lam$lambda_e, lambda_d = lam$lambda_d,
       psi00 = psi00, psi11 = psi11, benchmark_weight = benchmark_weight,
       utility_prior = utility_prior, u_low = u_low,
       u_benchmark = u_benchmark)
}
