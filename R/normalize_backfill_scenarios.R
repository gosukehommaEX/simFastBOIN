#' Scenarios of DLT and Response Probabilities
#'
#' @description
#'   Internal helper of \code{\link{sim_bf_boin_multi}} and
#'   \code{\link{sim_be_boin_multi}}. Check a list of scenarios, each a list
#'   with \code{p_true} and \code{p_resp} and optionally a \code{name}, and give
#'   every scenario a name.
#'
#' @param scenarios
#'   A list of scenarios.
#'
#' @return
#'   A list of lists with components \code{name}, \code{p_true} and
#'   \code{p_resp}.
#'
#' @noRd
normalize_backfill_scenarios <- function(scenarios) {
  if (!is.list(scenarios) || length(scenarios) == 0L) {
    stop("'scenarios' must be a non-empty list", call. = FALSE)
  }
  has_both <- vapply(scenarios, function(element) {
    is.list(element) && !is.null(element[["p_true"]]) &&
      !is.null(element[["p_resp"]])
  }, logical(1))
  if (!all(has_both)) {
    stop("each scenario must be a list with 'p_true' and 'p_resp' elements",
         call. = FALSE)
  }

  out <- normalize_scenarios(scenarios)
  for (i in seq_along(out)) {
    p_resp <- scenarios[[i]][["p_resp"]]
    if (!is.numeric(p_resp) || length(p_resp) != length(out[[i]]$p_true) ||
        any(!is.finite(p_resp)) || any(p_resp < 0) || any(p_resp > 1)) {
      stop("'p_resp' of scenario '", out[[i]]$name, "' must hold one ",
           "response probability between 0 and 1 per dose", call. = FALSE)
    }
    out[[i]]$p_resp <- p_resp
  }
  out
}
