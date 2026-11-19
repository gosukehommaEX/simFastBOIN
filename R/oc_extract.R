#' Summaries of One Design Result
#'
#' @description
#'   Internal helper of \code{\link{oc_results}} and \code{\link{oc_compare}}.
#'   Turn one result of a simulation function (a single scenario or a
#'   \code{_multi} result), of \code{\link{oc_3p3}}, or a named list of such
#'   single-scenario results, into the quantities shared by all the designs.
#'
#' @param x The result.
#' @param label Character scalar or \code{NULL}, the name of the design.
#'
#' @return A list with \code{design} (the name), \code{target_type}
#'   (\code{"MTD"} or \code{"OBD"}), \code{type} (\code{"single"} or
#'   \code{"combination"}) and \code{scenarios}, a named list with one element
#'   per scenario holding \code{p_true}, \code{sel}, \code{n_pts},
#'   \code{n_tox}, \code{n_resp} (\code{NULL} when not recorded),
#'   \code{no_sel}, \code{total_n_pts}, \code{total_n_tox},
#'   \code{total_n_resp} and \code{duration} (\code{NA} when not recorded).
#'
#' @noRd
oc_extract <- function(x, label = NULL) {
  single_classes <- c("boin_oc", "oc_3p3", "comb_boin_oc", "boin_et_oc",
                      "bf_boin_et_oc", "boin12_oc", "comb_boin12_oc",
                      "boin_etc_oc")
  is_single <- function(z) inherits(z, single_classes)

  if (is_single(x)) {
    results <- list("Scenario 1" = x)
  } else if (inherits(x, "list") || !is.object(x)) {
    if (!is.list(x) || length(x) == 0L || !all(vapply(x, is_single,
                                                       logical(1)))) {
      stop("each design must be the result of a simulation function, of ",
           "oc_3p3(), or a list of such results for single scenarios",
           call. = FALSE)
    }
    results <- x
    nms <- names(results)
    if (is.null(nms)) nms <- rep("", length(results))
    nms[!nzchar(nms)] <- paste("Scenario", which(!nzchar(nms)))
    names(results) <- nms
  } else if (is.list(x) && is.list(x[["results"]]) &&
             length(x[["results"]]) > 0L &&
             all(vapply(x[["results"]], is_single, logical(1)))) {
    results <- x[["results"]]
  } else {
    stop("each design must be the result of a simulation function, of ",
         "oc_3p3(), or a list of such results for single scenarios",
         call. = FALSE)
  }
  if (anyDuplicated(names(results))) {
    stop("the scenario names of a design must be unique", call. = FALSE)
  }

  first <- results[[1L]]
  design <- if (!is.null(label) && nzchar(label)) {
    label
  } else if (!is.null(first[["design"]])) {
    first[["design"]]
  } else if (inherits(first, "bf_boin_oc")) {
    "BF-BOIN"
  } else if (inherits(first, "be_boin_oc")) {
    "BE-BOIN"
  } else if (inherits(first, "tite_boin_oc")) {
    "TITE-BOIN"
  } else if (inherits(first, "comb_boin_oc")) {
    "BOIN combination"
  } else if (inherits(first, "oc_3p3")) {
    "3+3"
  } else {
    "BOIN"
  }
  classes <- vapply(results, function(z) class(z)[1L], character(1))
  if (length(unique(classes)) > 1L) {
    stop("the results of design '", design, "' come from different ",
         "functions", call. = FALSE)
  }
  target_type <- if (inherits(first, c("boin_oc", "oc_3p3", "comb_boin_oc"))) {
    "MTD"
  } else {
    "OBD"
  }

  scen <- lapply(results, function(z) {
    pick <- function(a, b = NULL) {
      if (!is.null(z[[a]])) z[[a]] else if (!is.null(b) && !is.null(z[[b]]))
        z[[b]] else NULL
    }
    no_sel <- as.numeric(pick("percent_no_mtd", "percent_no_obd"))
    duration <- as.numeric(pick("duration_mean", "duration"))
    list(
      p_true = z[["p_true"]],
      sel = z[["sel_percent"]],
      n_pts = z[["n_pts_dose"]],
      n_tox = z[["n_tox_dose"]],
      n_resp = z[["n_resp_dose"]],
      no_sel = no_sel,
      total_n_pts = as.numeric(z[["total_n_pts"]]),
      total_n_tox = as.numeric(z[["total_n_tox"]]),
      total_n_resp = if (is.null(z[["total_n_resp"]])) NA_real_ else
        as.numeric(z[["total_n_resp"]]),
      duration = if (length(duration) == 0L) NA_real_ else duration
    )
  })
  dims <- vapply(scen, function(z) length(z$sel), integer(1))
  if (length(unique(dims)) > 1L) {
    stop("the scenarios of design '", design, "' have different numbers ",
         "of doses", call. = FALSE)
  }

  list(design = design, target_type = target_type,
       type = if (is.matrix(first[["sel_percent"]])) "combination" else
         "single",
       scenarios = scen)
}
