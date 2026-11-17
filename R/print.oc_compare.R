#' Print a Comparison of Designs
#'
#' @description
#'   Print the comparison built by \code{\link{oc_compare}}: one block per
#'   scenario with one row per design.
#'
#' @param x
#'   An object of class \code{oc_compare}.
#'
#' @param digits
#'   Integer scalar. Number of decimals. Defaults to 1.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   \code{x}, invisibly.
#'
#' @examples
#' truth <- oc_scenarios(list(c(0.10, 0.25, 0.40, 0.55)), target = 0.30)
#' boin <- sim_boin(target = 0.30, p_true = c(0.10, 0.25, 0.40, 0.55),
#'                  n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 1)
#' print(oc_compare(BOIN = boin, "3+3" = oc_3p3(c(0.10, 0.25, 0.40, 0.55)),
#'                  truth = truth))
#'
#' @seealso \code{\link{oc_compare}}, \code{\link{plot.oc_compare}}
#'
#' @export
print.oc_compare <- function(x, digits = 1, ...) {
  if (!inherits(x, "oc_compare") ||
      !all(c("design", "target_type", "scenario", "pct_correct") %in%
             names(x))) {
    stop("'x' must be the result of oc_compare()", call. = FALSE)
  }
  x_in <- x
  x$design <- droplevels(x$design)
  x$scenario <- droplevels(x$scenario)
  check_count(digits, "digits", min_value = 0L)
  cols <- c(pct_correct = "Correct (%)",
            pct_patients_correct = "Pts at correct (%)",
            pct_overdose_selected = "Overdose sel. (%)",
            n_pts_overdose = "Pts at overdose",
            pct_no_selection = "No selection (%)",
            total_n_pts = "Patients", duration = "Duration")
  if (all(is.na(x$duration))) cols <- cols[names(cols) != "duration"]
  fmt <- function(v) ifelse(is.na(v), "-",
                            formatC(v, format = "f", digits = digits))
  cat("Comparison of", nlevels(x$design), "design(s) with the true doses\n")
  cat("Correct: the true MTD for designs selecting an MTD, the true OBD for",
      "designs selecting an OBD\n")
  if (!is.null(attr(x, "overdose_cutoff"))) {
    cat("Overdose: true DLT probability above", attr(x, "overdose_cutoff"),
        "\n")
  }
  for (s in levels(x$scenario)) {
    d <- x[x$scenario == s, , drop = FALSE]
    tab <- vapply(names(cols), function(m) fmt(d[[m]]),
                  character(nrow(d)))
    tab <- matrix(tab, nrow = nrow(d),
                  dimnames = list(paste0(as.character(d$design), " (",
                                         d$target_type, ")"),
                                  unname(cols)))
    cat("\n", s, "\n", sep = "")
    print(tab, quote = FALSE, right = TRUE)
  }
  invisible(x_in)
}
