#' Safety Stopping Boundary as a Two-Row Table
#'
#' @description
#'   Present the safety stopping boundary at the lowest dose as a compact table
#'   with one column per sample size, which is the shape usually wanted in a
#'   protocol or a report.
#'
#' @param x
#'   An object of class \code{boin_boundary}, built with \code{extrasafe = TRUE}.
#'
#' @param cohort_size
#'   Integer scalar or \code{NULL}. When supplied, only the sample sizes that are
#'   multiples of \code{cohort_size} are kept.
#'
#' @return
#'   A data frame with two rows, labelled by the quantity they hold, and one
#'   column per sample size. Sample sizes at which no number of DLTs triggers the
#'   rule are dropped.
#'
#' @details
#'   The same numbers are available one row per sample size from
#'   \code{as.data.frame()} on the boundary object. This function transposes them
#'   and labels the rows, which suits a table pasted into a document or rendered
#'   with \code{knitr::kable()} or \code{DT::datatable()}.
#'
#' @examples
#' bd <- boin_boundary(target = 0.30, max_n = 18, extrasafe = TRUE)
#'
#' boin_stopping_table(bd)
#'
#' # Only the sample sizes reached at the end of a cohort of three
#' boin_stopping_table(bd, cohort_size = 3)
#'
#' @seealso \code{\link{boin_boundary}}
#'
#' @export
boin_stopping_table <- function(x, cohort_size = NULL) {

  if (!inherits(x, "boin_boundary")) {
    stop("'x' must be a 'boin_boundary' object, as returned by boin_boundary()",
         call. = FALSE)
  }
  if (!isTRUE(x$extrasafe)) {
    stop("'x' was built without 'extrasafe = TRUE', so it holds no safety ",
         "stopping boundary", call. = FALSE)
  }

  keep <- which(!is.na(x$b_stop))
  if (!is.null(cohort_size)) {
    check_count(cohort_size, "cohort_size", 1L)
    keep <- keep[x$n[keep] %% as.integer(cohort_size) == 0L]
  }
  if (length(keep) == 0L) {
    stop("no sample size in the table has a safety stopping boundary",
         call. = FALSE)
  }

  out <- as.data.frame(rbind(x$n[keep], x$b_stop[keep]))
  colnames(out) <- as.character(x$n[keep])
  rownames(out) <- c(
    "Number of evaluable patients treated at the lowest dose level",
    "Stop the trial if # of DLT >="
  )
  out
}
