#' Print Operating Characteristics Across Scenarios
#'
#' @description
#'   Display the combined summary table produced by \code{\link{sim_boin_multi}}.
#'
#' @param x
#'   An object of class \code{boin_oc_multi}.
#'
#' @param kable
#'   Logical scalar. Return the table through \code{knitr::kable()} rather than
#'   printing it directly. Requires the \pkg{knitr} package. Defaults to
#'   \code{FALSE}.
#'
#' @param kable_format
#'   Character scalar passed to \code{knitr::kable()}, for example \code{"pipe"},
#'   \code{"html"} or \code{"latex"}. Defaults to \code{"pipe"}.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @examples
#' oc <- sim_boin_multi(
#'   target = 0.30,
#'   scenarios = list(
#'     list(name = "MTD at dose 3", p_true = c(0.05, 0.15, 0.30, 0.45, 0.60)),
#'     list(name = "All toxic",     p_true = c(0.35, 0.45, 0.55, 0.65, 0.75))
#'   ),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_trials = 200,
#'   seed = 123
#' )
#'
#' print(oc)
#'
#' @export
print.boin_oc_multi <- function(x, kable = FALSE, kable_format = "pipe", ...) {

  tab <- x$summary_table
  tab[] <- lapply(tab, function(column) {
    if (is.numeric(column)) {
      out <- format(column, nsmall = 1, trim = TRUE)
      out[is.na(column)] <- ""
      out
    } else {
      column
    }
  })

  if (kable) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package 'knitr' is needed when 'kable' is TRUE", call. = FALSE)
    }
    print(knitr::kable(tab, format = kable_format, row.names = FALSE))
    return(invisible(x))
  }

  first <- x$results[[1L]]
  cat("BOIN operating characteristics across ", length(x$scenario_names),
      " scenarios\n", sep = "")
  cat("  target DLT rate : ", format(first$target), "\n", sep = "")
  cat("  trials each     : ", first$n_trials, "\n", sep = "")
  cat("  max sample size : ", first$settings$max_total_pts, "\n", sep = "")
  cat("\n")

  print(tab, row.names = FALSE, right = TRUE)

  invisible(x)
}
