#' Print True Scenarios
#'
#' @description
#'   Print the scenarios tabulated by \code{\link{oc_scenarios}}: for a single
#'   agent, one block per scenario with the DLT and response probabilities and
#'   the utility of every dose and marks under the true MTDs and OBDs; for a
#'   combination, the same quantities as matrices (rows for drug A, columns
#'   for drug B) with the true MTDs and OBDs marked by an asterisk.
#'
#' @param x
#'   An object of class \code{oc_scenarios}.
#'
#' @param digits
#'   Integer scalar. Number of decimals of the probabilities; the utility is
#'   printed with one decimal. Defaults to 2.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   \code{x}, invisibly.
#'
#' @examples
#' truth <- oc_scenarios(
#'   list(list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'             p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))),
#'   target = 0.35, target_resp = 0.25
#' )
#' print(truth)
#'
#' @seealso \code{\link{oc_scenarios}}, \code{\link{plot.oc_scenarios}}
#'
#' @export
print.oc_scenarios <- function(x, digits = 2, ...) {
  if (!inherits(x, "oc_scenarios") ||
      !all(c("scenario", "p_true", "true_mtd") %in% names(x))) {
    stop("'x' must be the result of oc_scenarios()", call. = FALSE)
  }
  x_in <- x
  check_count(digits, "digits", min_value = 0L)
  # The type and size come from the columns, so that a subset prints too.
  type <- if ("dose_a" %in% names(x)) "combination" else "single"
  dims <- if (type == "single") length(unique(x$dose)) else
    c(max(x$dose_a), max(x$dose_b))
  has_resp <- all(c("p_resp", "utility", "true_obd") %in% names(x))
  x$scenario <- droplevels(x$scenario)
  fmt <- function(v, d) formatC(v, format = "f", digits = d)

  cat("True scenarios:", nlevels(x$scenario), "scenario(s),",
      if (type == "single") paste(dims, "doses") else
        paste(dims[1L], "x", dims[2L], "combinations"), "\n")
  if (!is.null(attr(x, "target"))) {
    cat("True MTD: DLT probability closest to", attr(x, "target"), "\n")
  }
  if (has_resp && !is.null(attr(x, "tox_limit"))) {
    cat("True OBD: highest utility among the doses with DLT probability <=",
        attr(x, "tox_limit"),
        if (!is.null(attr(x, "target_resp"))) {
          paste("and response probability >=", attr(x, "target_resp"))
        }, "\n")
    cat("Utilities: 100 (response, no DLT),", attr(x, "psi00"),
        "(neither),", attr(x, "psi11"), "(both), 0 (DLT, no response)\n")
  }

  for (s in levels(x$scenario)) {
    d <- x[x$scenario == s, , drop = FALSE]
    cat("\n", s, "\n", sep = "")
    if (type == "single") {
      tab <- rbind(p_true = fmt(d$p_true, digits))
      if (has_resp) {
        tab <- rbind(tab, p_resp = fmt(d$p_resp, digits),
                     utility = fmt(d$utility, 1L))
      }
      tab <- rbind(tab, "true MTD" = ifelse(d$true_mtd, "*", ""))
      if (has_resp) tab <- rbind(tab, "true OBD" = ifelse(d$true_obd, "*", ""))
      colnames(tab) <- d$dose
      print(tab, quote = FALSE, right = TRUE)
    } else {
      mat <- function(v, mark) {
        matrix(paste0(v, ifelse(mark, "*", " ")), dims[1L], dims[2L],
               dimnames = list(paste0("A", seq_len(dims[1L])),
                               paste0("B", seq_len(dims[2L]))))
      }
      cat("DLT probability (* true MTD)\n")
      print(mat(fmt(d$p_true, digits), d$true_mtd), quote = FALSE,
            right = TRUE)
      if (has_resp) {
        cat("Response probability\n")
        print(mat(fmt(d$p_resp, digits), rep(FALSE, nrow(d))), quote = FALSE,
              right = TRUE)
        cat("Utility (* true OBD)\n")
        print(mat(fmt(d$utility, 1L), d$true_obd), quote = FALSE,
              right = TRUE)
      }
    }
  }
  invisible(x_in)
}
