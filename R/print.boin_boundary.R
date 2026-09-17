#' Print BOIN Decision Boundaries
#'
#' @description
#'   Display the interval boundaries and the integer decision boundaries produced
#'   by \code{\link{boin_boundary}}.
#'
#' @param x
#'   An object of class \code{boin_boundary}.
#'
#' @param cohort_size
#'   Integer scalar or \code{NULL}. When supplied, only the sample sizes that are
#'   multiples of \code{cohort_size} are shown, which is how the boundaries are
#'   consulted during a trial with equally sized cohorts.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' bd <- boin_boundary(target = 0.30, max_n = 18, extrasafe = TRUE)
#' print(bd, cohort_size = 3)
#'
#' @export
print.boin_boundary <- function(x, cohort_size = NULL, ...) {

  # Validate before printing anything, so that a rejected argument does not
  # leave a half-written table behind.
  keep <- seq_along(x$n)
  if (!is.null(cohort_size)) {
    check_count(cohort_size, "cohort_size", 1L)
    keep <- keep[x$n %% as.integer(cohort_size) == 0L]
  }

  cat("BOIN decision boundaries\n")
  cat("  target DLT rate : ", format(x$target), "\n", sep = "")
  cat("  p_saf / p_tox   : ", format(x$p_saf), " / ", format(x$p_tox), "\n", sep = "")
  cat("  lambda_e        : ", format(round(x$lambda_e, 4)), "\n", sep = "")
  cat("  lambda_d        : ", format(round(x$lambda_d, 4)), "\n", sep = "")
  cat("  cutoff_eli      : ", format(x$cutoff_eli), "\n", sep = "")
  if (x$extrasafe) {
    cat("  safety cutoff   : ", format(x$cutoff_eli - x$offset), "\n", sep = "")
  }
  cat("\n")

  if (length(keep) == 0L) {
    cat("No sample size in the table is a multiple of the cohort size.\n")
    return(invisible(x))
  }

  rows <- list(
    "Number of patients treated" = x$n[keep],
    "Escalate if # of DLT <=" = x$b_esc[keep],
    "Deescalate if # of DLT >=" = x$b_deesc[keep],
    "Eliminate if # of DLT >=" = x$b_elim[keep]
  )
  if (x$extrasafe) {
    rows[["Stop at lowest dose if # of DLT >="]] <- x$b_stop[keep]
  }

  tab <- do.call(rbind, rows)
  colnames(tab) <- rep("", ncol(tab))
  print(tab, na.print = "NA", quote = FALSE)

  invisible(x)
}
