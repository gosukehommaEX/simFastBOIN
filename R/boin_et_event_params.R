#' Parameters of the Event Times of the Time-to-Event BOIN-ET Designs
#'
#' @description
#'   Shape and scale of the Weibull distribution of the time to an event, set
#'   so that an event occurs within the window with probability \code{p} and
#'   a share \code{late_fraction} of these events occurs in the second half of
#'   the window, as in \code{boinet::tite.boinet()} and
#'   \code{bfboinet::get.oc.backboinet()}. Unused entries, for a probability
#'   of zero or uniform event times, are set to one.
#'
#' @param p Numeric vector of event probabilities.
#' @param window Numeric scalar. Length of the assessment window.
#' @param late_fraction Numeric scalar strictly between 0 and 1.
#' @param event_time Character scalar, \code{"weibull"} or \code{"uniform"}.
#' @param name Character scalar used in error messages.
#'
#' @return A list with numeric vectors \code{shape} and \code{scale}.
#'
#' @noRd
boin_et_event_params <- function(p, window, late_fraction, event_time, name) {
  shape <- rep(1, length(p))
  scale <- rep(1, length(p))
  if (event_time == "weibull") {
    if (any(p >= 1)) {
      stop("'", name, "' must be below 1 with event_time = \"weibull\"",
           call. = FALSE)
    }
    use <- p > 0
    q <- p[use]
    shape[use] <- log(log(1 - q) / log(1 - q + late_fraction * q)) /
      log(1 / (1 - 0.5))
    scale[use] <- window / (-log(1 - q))^(1 / shape[use])
  }
  list(shape = shape, scale = scale)
}
