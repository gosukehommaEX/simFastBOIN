#' Print the Operating Characteristics of Several Designs
#'
#' @description
#'   Print the table built by \code{\link{oc_results}}: for every scenario,
#'   the chosen per-dose quantities with one row per design and one column
#'   per dose (for a combination, one matrix per design with rows for drug A
#'   and columns for drug B), followed by the quantities of the whole trial
#'   with one row per design.
#'
#' @param x
#'   An object of class \code{oc_results}.
#'
#' @param metric
#'   Character vector. Per-dose quantities to print, among
#'   \code{"sel_percent"}, \code{"n_pts"}, \code{"n_tox"} and
#'   \code{"n_resp"}. Defaults to \code{c("sel_percent", "n_pts")}.
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
#' boin <- sim_boin(target = 0.30, p_true = c(0.10, 0.25, 0.40, 0.55),
#'                  n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 1)
#' print(oc_results(BOIN = boin, "3+3" = oc_3p3(c(0.10, 0.25, 0.40, 0.55))))
#'
#' @seealso \code{\link{oc_results}}, \code{\link{plot.oc_results}}
#'
#' @export
print.oc_results <- function(x, metric = c("sel_percent", "n_pts"),
                             digits = 1, ...) {
  if (!inherits(x, "oc_results") ||
      !all(c("design", "scenario", "metric", "value") %in% names(x))) {
    stop("'x' must be the result of oc_results()", call. = FALSE)
  }
  dose_metrics <- c("sel_percent", "n_pts", "n_tox", "n_resp")
  if (!is.character(metric) || length(metric) == 0L ||
      !all(metric %in% dose_metrics)) {
    stop("'metric' must name per-dose quantities among ",
         paste0("'", dose_metrics, "'", collapse = ", "), call. = FALSE)
  }
  x_in <- x
  check_count(digits, "digits", min_value = 0L)
  labels <- c(sel_percent = "Selection (%)", n_pts = "Patients",
              n_tox = "DLTs", n_resp = "Responses",
              percent_no_selection = "No selection (%)",
              total_n_pts = "Patients", total_n_tox = "DLTs",
              total_n_resp = "Responses", duration = "Duration")
  type <- if ("dose_a" %in% names(x)) "combination" else "single"
  x$design <- droplevels(x$design)
  x$scenario <- droplevels(x$scenario)
  fmt <- function(v) ifelse(is.na(v), "-",
                            formatC(v, format = "f", digits = digits))

  cat("Operating characteristics of", nlevels(x$design), "design(s) in",
      nlevels(x$scenario), "scenario(s)\n")
  for (s in levels(x$scenario)) {
    d <- x[x$scenario == s, , drop = FALSE]
    cat("\n", s, "\n", sep = "")
    for (m in metric) {
      dm <- d[d$metric == m, , drop = FALSE]
      if (nrow(dm) == 0L) next
      if (type == "single") {
        doses <- sort(unique(dm$dose))
        tab <- matrix("-", nlevels(x$design), length(doses),
                      dimnames = list(levels(x$design), doses))
        tab[cbind(as.integer(dm$design), match(dm$dose, doses))] <-
          fmt(dm$value)
        tab <- tab[levels(x$design) %in% dm$design, , drop = FALSE]
        cat(labels[[m]], "by dose\n")
        print(tab, quote = FALSE, right = TRUE)
      } else {
        for (g in levels(x$design)) {
          dg <- dm[dm$design == g, , drop = FALSE]
          if (nrow(dg) == 0L) next
          tab <- matrix("-", max(dg$dose_a), max(dg$dose_b),
                        dimnames = list(paste0("A", seq_len(max(dg$dose_a))),
                                        paste0("B", seq_len(max(dg$dose_b)))))
          tab[cbind(dg$dose_a, dg$dose_b)] <- fmt(dg$value)
          cat(labels[[m]], "by combination,", g, "\n")
          print(tab, quote = FALSE, right = TRUE)
        }
      }
    }
    trial <- d[d$metric %in% names(labels)[5:9], , drop = FALSE]
    if (nrow(trial)) {
      ms <- intersect(names(labels)[5:9], as.character(trial$metric))
      tab <- matrix("-", nlevels(x$design), length(ms),
                    dimnames = list(levels(x$design), unname(labels[ms])))
      tab[cbind(as.integer(trial$design),
                match(as.character(trial$metric), ms))] <- fmt(trial$value)
      tab <- tab[levels(x$design) %in% trial$design, , drop = FALSE]
      cat("Per trial (averages)\n")
      print(tab, quote = FALSE, right = TRUE)
    }
  }
  invisible(x_in)
}
