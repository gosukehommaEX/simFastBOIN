#' Plot True Scenarios
#'
#' @description
#'   Draw the scenarios tabulated by \code{\link{oc_scenarios}}. For a single
#'   agent, one panel per scenario shows the DLT probability, the response
#'   probability and the utility divided by 100 against the dose, with the
#'   true MTDs labeled "MTD" on the DLT curve and the true OBDs labeled "OBD"
#'   on the utility curve. For a combination, one panel per scenario shows the
#'   DLT probability of every combination as a tile carrying its value, the
#'   true MTDs outlined with a solid line and the true OBDs with a dashed
#'   line.
#'
#' @param x
#'   An object of class \code{oc_scenarios}.
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
#' truth <- oc_scenarios(
#'   list(list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
#'             p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))),
#'   target = 0.35, target_resp = 0.25
#' )
#' plot(truth)
#'
#' @seealso \code{\link{oc_scenarios}}, \code{\link{print.oc_scenarios}}
#'
#' @export
plot.oc_scenarios <- function(x, text_size = 1, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed to plot the scenarios", call. = FALSE)
  }
  if (!inherits(x, "oc_scenarios") ||
      !all(c("scenario", "p_true", "true_mtd") %in% names(x))) {
    stop("'x' must be the result of oc_scenarios()", call. = FALSE)
  }
  if (!is.numeric(text_size) || length(text_size) != 1L ||
      !is.finite(text_size) || text_size <= 0) {
    stop("'text_size' must be a single positive number", call. = FALSE)
  }
  has_resp <- all(c("p_resp", "utility", "true_obd") %in% names(x))
  x$scenario <- droplevels(x$scenario)
  base <- ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10 * text_size),
      axis.title = ggplot2::element_text(size = 11 * text_size),
      strip.text = ggplot2::element_text(size = 11 * text_size,
                                         face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10 * text_size)
    )

  if (!"dose_a" %in% names(x)) {
    series <- c("DLT probability", "Response probability", "Utility / 100")
    long <- data.frame(scenario = x$scenario, dose = x$dose,
                       value = x$p_true,
                       series = factor(series[1L], levels = series))
    marks <- data.frame(scenario = x$scenario[x$true_mtd],
                        dose = x$dose[x$true_mtd],
                        value = x$p_true[x$true_mtd],
                        label = rep("MTD", sum(x$true_mtd)),
                        stringsAsFactors = FALSE)
    if (has_resp) {
      long <- rbind(
        long,
        data.frame(scenario = x$scenario, dose = x$dose, value = x$p_resp,
                   series = factor(series[2L], levels = series)),
        data.frame(scenario = x$scenario, dose = x$dose,
                   value = x$utility / 100,
                   series = factor(series[3L], levels = series))
      )
      marks <- rbind(marks,
                     data.frame(scenario = x$scenario[x$true_obd],
                                dose = x$dose[x$true_obd],
                                value = x$utility[x$true_obd] / 100,
                                label = rep("OBD", sum(x$true_obd)),
                                stringsAsFactors = FALSE))
    }
    long$series <- droplevels(long$series)
    p <- ggplot2::ggplot(long,
                         ggplot2::aes(x = dose, y = value, colour = series)) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::geom_text(data = marks,
                         ggplot2::aes(x = dose, y = value, label = label),
                         inherit.aes = FALSE, vjust = -0.9,
                         size = 3 * text_size, fontface = "bold")
    if (!is.null(attr(x, "target"))) {
      p <- p + ggplot2::geom_hline(yintercept = attr(x, "target"),
                                   linetype = "dotted", colour = "grey40")
    }
    p +
      ggplot2::scale_x_continuous(breaks = sort(unique(x$dose))) +
      ggplot2::scale_y_continuous(limits = c(0, 1.08)) +
      ggplot2::scale_colour_manual(
        values = c("DLT probability" = "#D03B3B",
                   "Response probability" = "#2A78D6",
                   "Utility / 100" = "#0CA30C"),
        name = NULL
      ) +
      ggplot2::facet_wrap(~ scenario) +
      ggplot2::labs(x = "Dose level", y = "Probability or utility / 100") +
      base
  } else {
    tiles <- data.frame(scenario = x$scenario,
                        dose_a = factor(x$dose_a), dose_b = factor(x$dose_b),
                        p_true = x$p_true,
                        label = formatC(x$p_true, format = "f", digits = 2),
                        stringsAsFactors = FALSE)
    if (has_resp) {
      tiles$label <- paste0(tiles$label, "\n",
                            formatC(x$p_resp, format = "f", digits = 2))
    }
    p <- ggplot2::ggplot(tiles,
                         ggplot2::aes(x = dose_b, y = dose_a, fill = p_true)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 3 * text_size) +
      ggplot2::geom_tile(data = tiles[x$true_mtd, , drop = FALSE],
                         fill = NA, colour = "black", linewidth = 1)
    if (has_resp) {
      p <- p + ggplot2::geom_tile(data = tiles[x$true_obd, , drop = FALSE],
                                  fill = NA, colour = "black",
                                  linewidth = 1, linetype = "dashed")
    }
    p +
      ggplot2::scale_fill_gradient(low = "#FFF5E6", high = "#D03B3B",
                                   limits = c(0, 1),
                                   name = "DLT probability") +
      ggplot2::facet_wrap(~ scenario) +
      ggplot2::labs(x = "Drug B level", y = "Drug A level",
                    caption = if (has_resp) {
                      paste("Tiles: DLT and response probabilities;",
                            "solid outline: true MTD; dashed: true OBD")
                    } else {
                      "Solid outline: true MTD"
                    }) +
      base
  }
}
