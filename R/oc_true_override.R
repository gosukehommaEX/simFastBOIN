#' True MTDs or OBDs Given by the User
#'
#' @description
#'   Internal helper of \code{\link{oc_scenarios}}. Check the argument
#'   \code{true_mtd} or \code{true_obd} and return, for every scenario, the
#'   doses given by the user, or \code{NULL} where the doses are to be
#'   computed.
#'
#' @param value The argument: \code{NULL}, or a list with one element per
#'   scenario, in the order of the scenarios or named by them.
#' @param name Character scalar, the name of the argument.
#' @param scenario_names Character vector of the scenario names.
#' @param dims Integer vector: the number of doses, or the numbers of levels
#'   of drug A and drug B.
#'
#' @return A list with one element per scenario: \code{NULL}, or a logical
#'   vector (single agent) or matrix (combination) marking the true doses.
#'
#' @noRd
oc_true_override <- function(value, name, scenario_names, dims) {
  out <- vector("list", length(scenario_names))
  if (is.null(value)) return(out)
  if (!is.list(value)) {
    stop("'", name, "' must be NULL or a list with one element per scenario",
         call. = FALSE)
  }
  if (!is.null(names(value)) && all(nzchar(names(value)))) {
    unknown <- setdiff(names(value), scenario_names)
    if (length(unknown)) {
      stop("'", name, "' names scenarios that do not exist: ",
           paste0("'", unknown, "'", collapse = ", "), call. = FALSE)
    }
    idx <- match(names(value), scenario_names)
  } else {
    if (length(value) != length(scenario_names)) {
      stop("'", name, "' must have one element per scenario", call. = FALSE)
    }
    idx <- seq_along(value)
  }
  for (k in seq_along(value)) {
    v <- value[[k]]
    if (is.null(v)) next
    if (length(dims) == 1L) {
      if (!is.numeric(v) || any(!is.finite(v)) || any(v != round(v)) ||
          any(v < 1) || any(v > dims)) {
        stop("each element of '", name, "' must hold dose levels between ",
             "1 and ", dims, call. = FALSE)
      }
      mark <- logical(dims)
      mark[v] <- TRUE
    } else {
      if (length(v) == 0L) v <- matrix(numeric(0), 0L, 2L)
      if (!is.numeric(v) || !is.matrix(v) || ncol(v) != 2L ||
          any(!is.finite(v)) || any(v != round(v)) || any(v < 1) ||
          any(v[, 1L] > dims[1L]) || any(v[, 2L] > dims[2L])) {
        stop("each element of '", name, "' must be a two-column matrix of ",
             "the levels of drug A and drug B", call. = FALSE)
      }
      mark <- matrix(FALSE, dims[1L], dims[2L])
      mark[v] <- TRUE
    }
    out[[idx[k]]] <- mark
  }
  out
}
