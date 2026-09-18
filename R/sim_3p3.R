#' Simulate the 3+3 Design
#'
#' @description
#'   Simulate the traditional 3+3 dose-finding design and summarize the same
#'   quantities that \code{\link{oc_3p3}} computes exactly. The two are meant to
#'   be used together: \code{oc_3p3()} gives the answer without Monte Carlo
#'   error, and this function confirms it from independent trials.
#'
#' @param p_true
#'   Numeric vector. True DLT probability at each dose level, in increasing dose
#'   order.
#'
#' @param n_trials
#'   Integer scalar. Number of trials to simulate. Defaults to 10000.
#'
#' @param mtd_rule
#'   Character scalar, either \code{"previous"} or \code{"expand"}. See
#'   \code{\link{oc_3p3}}. Defaults to \code{"previous"}.
#'
#' @param start_dose
#'   Integer scalar. Dose level for the first cohort. Defaults to 1.
#'
#' @param overdose_cutoff
#'   Numeric scalar or \code{NULL}. Doses whose true DLT probability exceeds this
#'   value count as overdoses in the \code{overdose} component of the result.
#'   Defaults to \code{NULL}, which uses one third.
#'
#' @param seed
#'   Integer scalar or \code{NULL}. Random seed. The state of the calling session
#'   is restored on exit. Defaults to 123.
#'
#' @return
#'   An object of class \code{oc_3p3}, with the components described in
#'   \code{\link{oc_3p3}}. The \code{method} component is the number of trials
#'   rather than \code{"exact"}.
#'
#' @details
#'   The design is described under \code{\link{oc_3p3}}. Since the operating
#'   characteristics of the 3+3 design can be written in closed form, this
#'   function is a check on that closed form rather than the way to obtain the
#'   numbers. Simulation is the only option when a quantity is wanted that the
#'   closed form does not provide, such as the distribution of the sample size
#'   rather than its mean.
#'
#' @examples
#' exact <- oc_3p3(p_true = c(0.30, 0.48, 0.67))
#' simulated <- sim_3p3(p_true = c(0.30, 0.48, 0.67), n_trials = 2000, seed = 1)
#'
#' round(exact$sel_percent, 2)
#' round(simulated$sel_percent, 2)
#'
#' @seealso \code{\link{oc_3p3}}, \code{\link{sim_boin}}
#'
#' @importFrom stats rbinom
#'
#' @export
sim_3p3 <- function(p_true, n_trials = 10000, mtd_rule = c("previous", "expand"),
                    start_dose = 1, overdose_cutoff = NULL, seed = 123) {

  call_expr <- match.call()
  check_p_true(p_true)
  mtd_rule <- match.arg(mtd_rule)
  check_count(n_trials, "n_trials")
  check_count(start_dose, "start_dose")
  n_trials <- as.integer(n_trials)
  start_dose <- as.integer(start_dose)
  n_doses <- length(p_true)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of dose levels", call. = FALSE)
  }
  if (is.null(overdose_cutoff)) overdose_cutoff <- 1 / 3
  check_scalar_prob(overdose_cutoff, "overdose_cutoff")

  old_seed <- get_random_seed()
  on.exit(restore_random_seed(old_seed), add = TRUE)
  if (!is.null(seed)) set.seed(seed)

  above <- p_true > overdose_cutoff

  n_total <- numeric(n_doses)
  y_total <- numeric(n_doses)
  selected <- integer(n_trials)
  n_above_trial <- numeric(n_trials)
  n_trial <- numeric(n_trials)

  for (trial in seq_len(n_trials)) {

    n <- integer(n_doses)
    y <- integer(n_doses)
    dose <- start_dose
    toxic <- 0L
    ran_out <- FALSE

    repeat {
      first <- rbinom(1L, 3L, p_true[dose])
      n[dose] <- n[dose] + 3L
      y[dose] <- y[dose] + first
      moves_up <- first == 0L
      if (first == 1L) {
        second <- rbinom(1L, 3L, p_true[dose])
        n[dose] <- n[dose] + 3L
        y[dose] <- y[dose] + second
        moves_up <- second == 0L
      }
      if (moves_up) {
        if (dose == n_doses) {
          ran_out <- TRUE
          break
        }
        dose <- dose + 1L
      } else {
        toxic <- dose
        break
      }
    }

    if (mtd_rule == "previous") {
      mtd <- if (ran_out) {
        n_doses
      } else if (toxic > start_dose) {
        toxic - 1L
      } else {
        0L
      }
    } else {
      candidate <- if (ran_out) n_doses else toxic - 1L
      mtd <- 0L
      while (candidate >= start_dose) {
        if (n[candidate] < 6L) {
          extra <- rbinom(1L, 6L - n[candidate], p_true[candidate])
          y[candidate] <- y[candidate] + extra
          n[candidate] <- 6L
        }
        if (y[candidate] <= 1L) {
          mtd <- candidate
          break
        }
        candidate <- candidate - 1L
      }
    }

    n_total <- n_total + n
    y_total <- y_total + y
    selected[trial] <- mtd
    n_trial[trial] <- sum(n)
    n_above_trial[trial] <- sum(n[above])
  }

  n_pts_dose <- n_total / n_trials
  n_tox_dose <- y_total / n_trials
  sel_percent <- vapply(seq_len(n_doses),
                        function(d) mean(selected == d) * 100, numeric(1))
  percent_no_mtd <- mean(selected == 0L) * 100

  any_selected <- selected > 0L
  mtd_above <- logical(n_trials)
  mtd_above[any_selected] <- above[selected[any_selected]]

  overdose <- list(
    cutoff = overdose_cutoff,
    doses = which(above),
    pct_patients = sum(n_above_trial) / sum(n_trial) * 100,
    avg_n_patients = mean(n_above_trial),
    pct_trials_any = mean(n_above_trial > 0) * 100,
    pct_trials_mtd_above = mean(mtd_above) * 100,
    pct_mtd_above_when_selected = if (any(any_selected)) {
      mean(mtd_above[any_selected]) * 100
    } else {
      NA_real_
    }
  )

  structure(
    list(
      sel_percent = sel_percent,
      percent_no_mtd = percent_no_mtd,
      n_pts_dose = n_pts_dose,
      n_tox_dose = n_tox_dose,
      total_n_pts = sum(n_pts_dose),
      total_n_tox = sum(n_tox_dose),
      overdose = overdose,
      p_true = p_true,
      n_doses = n_doses,
      mtd_rule = mtd_rule,
      start_dose = start_dose,
      method = "simulated",
      n_trials = n_trials,
      call = call_expr
    ),
    class = "oc_3p3"
  )
}
