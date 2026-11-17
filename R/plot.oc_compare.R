#' Plot a Comparison of Designs
#'
#' @description
#'   Draw the comparison built by \code{\link{oc_compare}}: one panel per
#'   quantity, with a point for every design and scenario, the points of a
#'   design joined by a line.
#'
#' @param x
#'   An object of class \code{oc_compare}.
#'
#' @param metrics
#'   Character vector. Columns of \code{x} to draw. Defaults to
#'   \code{c("pct_correct", "pct_patients_correct", "pct_overdose_selected",
#'   "n_pts_overdose")}.
#'
#' @param text_size
#'   Numeric scalar. Multiplier applied to every text size in the figure.
#'   Defaults to 1.
#'
#' @param ...
#'   Further arguments, currently ignored.
#'
#' @return
#'   A \pkg{ggplot2} object.
#'
#' @details
#'   Requires the \pkg{ggplot2} package, which is only suggested by
#'   \pkg{simFastBOIN}. Missing values (for example the correct selection in
#'   a scenario without a true OBD) are left out.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' truth <- oc_scenarios(list(c(0.10, 0.25, 0.40, 0.55)), target = 0.30)
#' boin <- sim_boin(target = 0.30, p_true = c(0.10, 0.25, 0.40, 0.55),
#'                  n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 1)
#' plot(oc_compare(BOIN = boin, "3+3" = oc_3p3(c(0.10, 0.25, 0.40, 0.55)),
#'                 truth = truth))
#'
#' @seealso \code{\link{oc_compare}}, \code{\link{print.oc_compare}}
#'
#' @export
plot.oc_compare <- function(x,
                            metrics = c("pct_correct", "pct_patients_correct",
                                        "pct_overdose_selected",
                                        "n_pts_overdose"),
                            text_size = 1, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed to plot the comparison", call. = FALSE)
  }
  if (!inherits(x, "oc_compare") ||
      !all(c("design", "scenario", "pct_correct") %in% names(x))) {
    stop("'x' must be the result of oc_compare()", call. = FALSE)
  }
  labels <- c(pct_correct = "Correct selection (%)",
              n_pts_correct = "Patients at the correct dose",
              pct_patients_correct = "Patients at the correct dose (%)",
              pct_overdose_selected = "Overdose selected (%)",
              n_pts_overdose = "Patients at overdoses",
              pct_patients_overdose = "Patients at overdoses (%)",
              pct_no_selection = "No selection (%)",
              total_n_pts = "Patients", duration = "Duration")
  if (!is.character(metrics) || length(metrics) == 0L ||
      !all(metrics %in% names(labels))) {
    stop("'metrics' must name columns among ",
         paste0("'", names(labels), "'", collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(text_size) || length(text_size) != 1L ||
      !is.finite(text_size) || text_size <= 0) {
    stop("'text_size' must be a single positive number", call. = FALSE)
  }
  long <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(design = x$design, scenario = x$scenario,
               metric = labels[[m]], value = x[[m]],
               stringsAsFactors = FALSE)
  }))
  long$metric <- factor(long$metric, levels = unname(labels[metrics]))
  long <- long[!is.na(long$value), , drop = FALSE]
  p <- ggplot2::ggplot(long, ggplot2::aes(x = scenario, y = value,
                                          colour = design, group = design))
  # Lines join the scenarios of a design only when there are several.
  if (nlevels(droplevels(long$scenario)) > 1L) {
    p <- p + ggplot2::geom_line(linewidth = 0.6)
  }
  p +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, colour = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10 * text_size),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(size = 11 * text_size,
                                         face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10 * text_size)
    )
}
