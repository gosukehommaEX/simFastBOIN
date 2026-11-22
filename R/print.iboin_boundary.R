#' Print the Decision Boundaries of the iBOIN Design
#'
#' @description
#'   Print the prior information and, for every dose, the numbers of DLTs
#'   that lead to escalation and de-escalation, in the layout of Table 1 of
#'   Zhou et al. (2021).
#'
#' @param x
#'   An object of class \code{iboin_boundary}, from
#'   \code{\link{iboin_boundary}}.
#'
#' @param cohort_size
#'   Integer scalar or \code{NULL}. When given, only the numbers of patients
#'   that are multiples of the cohort size are shown. Defaults to \code{NULL},
#'   every number of patients.
#'
#' @param ...
#'   Unused.
#'
#' @return
#'   The object \code{x}, invisibly.
#'
#' @details
#'   An escalation boundary printed as \code{NA} means that no number of
#'   DLTs leads to escalation at that dose and sample size, and a
#'   de-escalation boundary printed as \code{NA} that none leads to
#'   de-escalation.
#'
#' @examples
#' bd <- iboin_boundary(
#'   target = 0.30, skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
#'   prior_ess = 3, max_n = 30
#' )
#' print(bd, cohort_size = 3)
#'
#' @seealso \code{\link{iboin_boundary}}
#'
#' @export
print.iboin_boundary <- function(x, cohort_size = NULL, ...) {

  # Validate before printing anything, so that a rejected argument does not
  # leave a half-written table behind.
  keep <- seq_along(x$n)
  if (!is.null(cohort_size)) {
    check_count(cohort_size, "cohort_size", 1L)
    keep <- keep[x$n %% as.integer(cohort_size) == 0L]
  }

  n_doses <- length(x$skeleton)
  dose_names <- paste0("DL", seq_len(n_doses))

  cat("iBOIN decision boundaries\n")
  cat("  target DLT rate : ", format(x$target), "\n", sep = "")
  cat("  p_saf / p_tox   : ", format(x$p_saf), " / ", format(x$p_tox), "\n",
      sep = "")
  cat("  cutoff_eli      : ", format(x$cutoff_eli), "\n", sep = "")
  if (x$extrasafe) {
    cat("  safety cutoff   : ", format(x$cutoff_eli - x$offset), "\n",
        sep = "")
  }
  cat("  robust prior    : ", if (x$robust) "yes" else "no", "\n", sep = "")
  cat("\n")

  prior_tab <- cbind(
    "Skeleton" = x$skeleton,
    "Prior ESS" = x$prior_ess,
    "Pr(target)" = round(x$prior[, "target"], 3),
    "Pr(p_saf)" = round(x$prior[, "p_saf"], 3),
    "Pr(p_tox)" = round(x$prior[, "p_tox"], 3)
  )
  rownames(prior_tab) <- dose_names
  cat("Prior information:\n")
  print(prior_tab)
  cat("\n")

  if (length(keep) == 0L) {
    cat("No sample size in the table is a multiple of the cohort size.\n")
    return(invisible(x))
  }

  n_keep <- x$n[keep]
  esc <- x$b_esc[, keep, drop = FALSE]
  deesc <- x$b_deesc[, keep, drop = FALSE]
  esc[esc < 0L] <- NA_integer_
  deesc[deesc > matrix(n_keep, n_doses, length(keep), byrow = TRUE)] <-
    NA_integer_

  tab <- matrix(NA_integer_, 2L * n_doses + 1L, length(keep))
  tab[2L * seq_len(n_doses) - 1L, ] <- esc
  tab[2L * seq_len(n_doses), ] <- deesc
  tab[2L * n_doses + 1L, ] <- x$b_elim[keep]
  rownames(tab) <- c(
    as.vector(rbind(paste(dose_names, "escalate if # of DLT <="),
                    paste(dose_names, "de-escalate if # of DLT >="))),
    "Eliminate if # of DLT >="
  )
  colnames(tab) <- n_keep
  cat("Number of patients treated at the current dose (columns):\n")
  print(tab, na.print = "NA")

  invisible(x)
}
