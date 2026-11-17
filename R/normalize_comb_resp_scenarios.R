#' Normalize Drug-Combination Scenarios with Response Probabilities
#'
#' @param scenarios List of scenarios, each a list with the matrices
#'   \code{p_true} and \code{p_resp} and optionally \code{name}.
#'
#' @return A list of \code{list(name, p_true, p_resp)}, the names and the
#'   DLT probabilities as from \code{normalize_comb_scenarios()}.
#'
#' @noRd
normalize_comb_resp_scenarios <- function(scenarios) {
  out <- normalize_comb_scenarios(scenarios)
  for (i in seq_along(out)) {
    element <- scenarios[[i]]
    if (!is.list(element) || is.null(element$p_resp)) {
      stop("each scenario must be a list with the matrices 'p_true' and ",
           "'p_resp'", call. = FALSE)
    }
    check_p_resp_matrix(element$p_resp, out[[i]]$p_true)
    out[[i]]$p_resp <- element$p_resp
  }
  out
}
