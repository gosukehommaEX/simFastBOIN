#' Normalize Drug-Combination Scenarios
#'
#' @description
#'   Accept either a list of \code{list(name, p_true)} elements or a named list
#'   of matrices, check every matrix and require them all to have the same
#'   dimensions. Internal.
#'
#' @param scenarios
#'   A non-empty list of scenarios.
#'
#' @return
#'   A list of \code{list(name, p_true)} elements.
#'
#' @keywords internal
#' @noRd
normalize_comb_scenarios <- function(scenarios) {
  if (!is.list(scenarios) || length(scenarios) == 0L) {
    stop("'scenarios' must be a non-empty list", call. = FALSE)
  }
  nms <- names(scenarios)
  out <- vector("list", length(scenarios))
  for (i in seq_along(scenarios)) {
    element <- scenarios[[i]]
    if (is.matrix(element)) {
      name <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste("Scenario", i)
      out[[i]] <- list(name = name, p_true = element)
    } else if (is.list(element) && !is.null(element$p_true)) {
      name <- if (!is.null(element$name)) {
        as.character(element$name)
      } else if (!is.null(nms) && nzchar(nms[i])) {
        nms[i]
      } else {
        paste("Scenario", i)
      }
      out[[i]] <- list(name = name, p_true = element$p_true)
    } else {
      stop("each scenario must be a matrix or a list with a 'p_true' element",
           call. = FALSE)
    }
    check_p_true_matrix(out[[i]]$p_true)
  }
  first_dim <- dim(out[[1L]]$p_true)
  same_dim <- vapply(out, function(z) identical(dim(z$p_true), first_dim),
                     logical(1))
  if (!all(same_dim)) {
    stop("all scenarios must have the same numbers of levels of both drugs",
         call. = FALSE)
  }
  if (anyDuplicated(vapply(out, function(z) z$name, character(1)))) {
    stop("scenario names must be unique", call. = FALSE)
  }
  out
}
