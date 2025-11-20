#' Generate Time-to-Event Data for TITE-BOIN Trials
#'
#' @description
#' Simulates time-to-event data for dose-limiting toxicity (DLT) in phase I trials
#' using time-to-event Bayesian optimal interval (TITE-BOIN) designs.
#' Supports multiple distributions for modeling the time to toxicity.
#'
#' @param n Integer. Number of patients to simulate.
#' @param prob Numeric. True DLT probability (scalar value between 0 and 1).
#' @param alpha1 Numeric. Late-onset parameter (0 < alpha1 < 1). Represents the
#'   proportion of toxicity that occurs in the last fraction (alpha2) of the
#'   assessment window. Default is 0.5.
#' @param alpha2 Numeric. Weibull shape parameter (alpha2 > 0). Defines the
#'   last fraction of the assessment window where alpha1 proportion of toxicity
#'   occurs. Default is 0.5.
#' @param distribution Character. Distribution for time-to-toxicity. One of:
#'   \itemize{
#'     \item "weibull": Weibull distribution (default, most commonly used)
#'     \item "log-logistic": Log-logistic distribution
#'     \item "uniform": Uniform distribution (simplest assumption)
#'   }
#' @param maxt Numeric. Length of the DLT assessment window (e.g., 28 days).
#'   Toxicities occurring after this time are not considered as DLTs.
#'   Default is 1 (standardized time unit).
#'
#' @return A list with two components:
#'   \describe{
#'     \item{tox}{Integer vector. DLT indicator (1 = DLT observed within assessment window,
#'                0 = no DLT or DLT occurred after assessment window)}
#'     \item{t_tox}{Numeric vector. Time to toxicity or censoring. For patients without
#'                  DLT within the assessment window, this is set to maxt.}
#'   }
#'
#' @details
#' The function generates time-to-event data based on the specified distribution:
#'
#' \strong{Uniform Distribution:}
#' DLT times are uniformly distributed over \code{[0, maxt]}.
#'
#' \strong{Weibull Distribution:}
#' Parameters are calculated such that:
#' \itemize{
#'   \item The cumulative DLT probability at time maxt equals prob
#'   \item The cumulative DLT probability at time (1 - alpha2) * maxt equals (1 - alpha1) * prob
#' }
#'
#' \strong{Log-logistic Distribution:}
#' Similar parameterization to Weibull, with shape and scale parameters derived
#' from alpha1, alpha2, and prob.
#'
#' @references
#' Yuan, Y., Lin, R., Li, D., Nie, L., & Warren, K. E. (2018).
#' Time-to-Event Bayesian Optimal Interval Design to Accelerate Phase I Trials.
#' \emph{Clinical Cancer Research}, 24(20), 4921-4930.
#'
#' Lin, R., & Yuan, Y. (2020).
#' Time-to-event model-assisted designs for dose-finding trials with delayed toxicity.
#' \emph{Biostatistics}, 21(4), 807-824.
#'
#' @importFrom stats rbinom runif
#'
#' @examples
#' # Example 1: Generate data for 10 patients with 30% DLT rate
#' # using Weibull distribution (default) and 28-day assessment window
#' set.seed(123)
#' data1 <- rtite(n = 10, prob = 0.30, maxt = 28)
#' data1$tox  # DLT indicators
#' data1$t_tox  # Time to DLT or censoring
#'
#' # Example 2: Use uniform distribution
#' set.seed(123)
#' data2 <- rtite(n = 10, prob = 0.30, distribution = "uniform", maxt = 28)
#'
#' # Example 3: Custom alpha1 and alpha2 parameters
#' # Assume 80% of toxicity occurs in the last 30% of assessment window
#' set.seed(123)
#' data3 <- rtite(
#'   n = 20,
#'   prob = 0.25,
#'   alpha1 = 0.8,
#'   alpha2 = 0.3,
#'   distribution = "weibull",
#'   maxt = 21
#' )
#'
#' # Example 4: Compare distributions
#' set.seed(123)
#' n_sim <- 1000
#' data_weibull <- rtite(n_sim, prob = 0.30, distribution = "weibull", maxt = 28)
#' data_uniform <- rtite(n_sim, prob = 0.30, distribution = "uniform", maxt = 28)
#' data_loglogistic <- rtite(n_sim, prob = 0.30, distribution = "log-logistic", maxt = 28)
#'
#' # Compare DLT rates (should all be close to 0.30)
#' mean(data_weibull$tox)
#' mean(data_uniform$tox)
#' mean(data_loglogistic$tox)
#'
#' @export
rtite <- function(
    n,
    prob,
    alpha1 = 0.5,
    alpha2 = 0.5,
    distribution = "weibull",
    maxt = 1
) {

  # Input validation
  if (!distribution %in% c("weibull", "log-logistic", "uniform")) {
    stop("distribution must be one of: 'weibull', 'log-logistic', 'uniform'")
  }
  if (alpha1 <= 0 || alpha1 >= 1) {
    stop("alpha1 must be between 0 and 1 (exclusive)")
  }
  if (alpha2 <= 0 || alpha2 >= 1) {
    stop("alpha2 must be between 0 and 1 (exclusive)")
  }
  if (maxt <= 0) {
    stop("maxt must be positive")
  }
  if (n <= 0) {
    stop("n must be positive")
  }

  # Extract DLT probability
  pi <- prob

  if (pi <= 0 || pi >= 1) {
    stop("DLT probability (prob) must be between 0 and 1 (exclusive)")
  }

  # Generate time-to-event data based on distribution
  if (distribution == "uniform") {
    # Uniform distribution: DLT times uniformly distributed over [0, maxt]
    tox_indicator <- rbinom(n, 1, pi)
    t_temp <- runif(n, 0, maxt)
    t_tox <- ifelse(tox_indicator == 1, t_temp, maxt)
    tox_obs <- tox_indicator

  } else if (distribution == "weibull") {
    # Weibull distribution with parameters derived from alpha1 and alpha2
    # pihalft: cumulative DLT probability at time (1 - alpha2) * maxt
    pihalft <- (1 - alpha1) * pi

    # Calculate Weibull shape parameter (alpha)
    # Based on the relationship between pi and pihalft
    alpha <- log(log(1 - pi) / log(1 - pihalft)) / log(1 / (1 - alpha2))

    # Calculate Weibull scale parameter (lambda)
    lambda <- -log(1 - pi) / (maxt ^ alpha)

    # Generate Weibull-distributed time-to-toxicity
    u <- runif(n)
    t_tox <- (-log(1 - u) / lambda) ^ (1 / alpha)

    # Determine if DLT occurred within assessment window
    tox_obs <- as.integer(t_tox <= maxt)

    # Censor times beyond maxt
    t_tox[tox_obs == 0] <- maxt

  } else if (distribution == "log-logistic") {
    # Log-logistic distribution with parameters derived from alpha1 and alpha2
    pihalft <- (1 - alpha1) * pi

    # Calculate log-logistic shape parameter (alpha)
    alpha <- log(log(1 - pi) / log(1 - pihalft)) / log(1 / (1 - alpha2))

    # Calculate log-logistic scale parameter
    # Avoid division by zero when alpha = 1
    if (abs(alpha - 1) < 1e-10) {
      # When alpha ≈ 1, use limiting case
      scale_param <- maxt / log((1 - pihalft) / (1 - pi))
    } else {
      denom <- alpha ^ (1 / alpha) - alpha ^ (-1 / alpha)
      if (abs(denom) < 1e-10) {
        # Fallback to avoid numerical issues
        scale_param <- maxt
      } else {
        scale_param <- maxt / denom
      }
    }

    # Generate log-logistic-distributed time-to-toxicity
    u <- runif(n)
    t_tox <- scale_param * ((u / (1 - u)) ^ (1 / alpha))

    # Determine if DLT occurred within assessment window
    tox_obs <- as.integer(t_tox <= maxt)

    # Censor times beyond maxt
    t_tox[tox_obs == 0] <- maxt
  }

  # Return results
  return(list(
    tox = as.integer(tox_obs),
    t_tox = t_tox
  ))
}
