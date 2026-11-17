#' Plot the Operating Characteristics of Several Designs
#'
#' @description
#'   Draw one per-dose quantity of the table built by
#'   \code{\link{oc_results}}. For a single agent, one panel per scenario
#'   shows bars for every dose, colored by design and labeled with their
#'   values. For a combination, one panel per design and scenario shows the
#'   value of every combination as a tile.
#'
#' @param x
#'   An object of class \code{oc_results}.
#'
#' @param metric
#'   Character scalar. The per-dose quantity: \code{"sel_percent"} (the
#'   default), \code{"n_pts"}, \code{"n_tox"} or \code{"n_resp"}.
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
#'   \pkg{simFastBOIN}.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' boin <- sim_boin(target = 0.30, p_true = c(0.10, 0.25, 0.40, 0.55),
#'                  n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 1)
#' res <- oc_results(BOIN = boin, "3+3" = oc_3p3(c(0.10, 0.25, 0.40, 0.55)))
#' plot(res)
#' plot(res, metric = "n_pts")
#'
#' @seealso \code{\link{oc_results}}, \code{\link{print.oc_results}}
#'
#' @export
plot.oc_results <- function(x, metric = "sel_percent", text_size = 1, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed to plot the results", call. = FALSE)
  }
  if (!inherits(x, "oc_results") ||
      !all(c("design", "scenario", "metric", "value") %in% names(x))) {
    stop("'x' must be the result of oc_results()", call. = FALSE)
  }
  dose_metrics <- c("sel_percent", "n_pts", "n_tox", "n_resp")
  if (!is.character(metric) || length(metric) != 1L ||
      !metric %in% dose_metrics) {
    stop("'metric' must be one of ",
         paste0("'", dose_metrics, "'", collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(text_size) || length(text_size) != 1L ||
      !is.finite(text_size) || text_size <= 0) {
    stop("'text_size' must be a single positive number", call. = FALSE)
  }
  d <- x[x$metric == metric, , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("no design records '", metric, "'", call. = FALSE)
  }
  d <- as.data.frame(d)
  d$label <- formatC(d$value, format = "f", digits = 1)
  y_label <- c(sel_percent = "Selection (%)", n_pts = "Patients",
               n_tox = "DLTs", n_resp = "Responses")[[metric]]
  base <- ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10 * text_size),
      axis.title = ggplot2::element_text(size = 11 * text_size),
      strip.text = ggplot2::element_text(size = 11 * text_size,
                                         face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10 * text_size)
    )

  if (!"dose_a" %in% names(d)) {
    d$dose <- factor(d$dose, levels = sort(unique(d$dose)))
    dodge <- ggplot2::position_dodge(width = 0.9)
    ggplot2::ggplot(d, ggplot2::aes(x = dose, y = value, fill = design)) +
      ggplot2::geom_col(position = dodge, width = 0.85) +
      ggplot2::geom_text(ggplot2::aes(label = label), position = dodge,
                         vjust = -0.3, size = 2.5 * text_size) +
      ggplot2::facet_wrap(~ scenario) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.1))
      ) +
      ggplot2::labs(x = "Dose level", y = y_label, fill = NULL) +
      base
  } else {
    d$dose_a <- factor(d$dose_a)
    d$dose_b <- factor(d$dose_b)
    ggplot2::ggplot(d, ggplot2::aes(x = dose_b, y = dose_a, fill = value)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 3 * text_size) +
      ggplot2::scale_fill_gradient(low = "#F2F6FC", high = "#2A78D6",
                                   name = y_label) +
      ggplot2::facet_grid(design ~ scenario) +
      ggplot2::labs(x = "Drug B level", y = "Drug A level") +
      base
  }
}
