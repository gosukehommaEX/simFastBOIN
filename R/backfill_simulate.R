#' Simulate Trials of a BOIN Design with Backfilling
#'
#' @description
#'   Internal workhorse of \code{\link{bf_boin_simulate}} and
#'   \code{\link{be_boin_simulate}}. Validate the arguments, build the decision
#'   boundaries and run the C++ engine.
#'
#' @param design
#'   Character scalar, \code{"bf_boin"} or \code{"be_boin"}.
#'
#' @param ...
#'   The arguments of \code{\link{bf_boin_simulate}} or
#'   \code{\link{be_boin_simulate}}, already matched.
#'
#' @return
#'   An object of class \code{c("bf_boin_trials", "backfill_trials")} or
#'   \code{c("be_boin_trials", "backfill_trials")}, as described in
#'   \code{\link{bf_boin_simulate}}.
#'
#' @noRd
backfill_simulate <- function(design, target, p_true, p_resp, n_cohort,
                              cohort_size, window, accrual_rate, n_cap,
                              backfill_dose, conflict_dose, no_slot, accrual,
                              dlt_time, late_fraction, resp_window,
                              resp_late_fraction, resp_cor, max_pending_ratio,
                              min_completed, min_follow_up, n_trials,
                              start_dose, n_earlystop, p_saf, p_tox,
                              cutoff_eli, extrasafe, offset, stay_on_1_of_3,
                              n_earlystop_rule, seed) {

  imputation <- identical(design, "be_boin")
  if (imputation) {
    rules <- tite_rules("imputation", max_pending_ratio, min_completed,
                        min_follow_up)
  }

  check_p_true(p_true)
  if (!is.numeric(p_resp) || length(p_resp) != length(p_true) ||
      any(!is.finite(p_resp)) || any(p_resp < 0) || any(p_resp > 1)) {
    stop("'p_resp' must be a vector of response probabilities between 0 and 1, ",
         "one per dose in 'p_true'", call. = FALSE)
  }
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(n_earlystop, "n_earlystop", 1L)
  check_count(start_dose, "start_dose", 1L)
  check_count(n_cap, "n_cap", 1L)
  check_flag(extrasafe, "extrasafe")
  check_flag(stay_on_1_of_3, "stay_on_1_of_3")
  for (arg in c("window", "accrual_rate", "resp_window")) {
    value <- get(arg)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) {
      stop("'", arg, "' must be a single positive number", call. = FALSE)
    }
  }
  check_scalar_prob(late_fraction, "late_fraction")
  check_scalar_prob(resp_late_fraction, "resp_late_fraction")
  check_scalar_prob(resp_cor, "resp_cor", lower = -1, upper = 1)
  if (dlt_time == "weibull" && any(p_true >= 1)) {
    stop("'p_true' must be below 1 when 'dlt_time' is \"weibull\"", call. = FALSE)
  }
  if (dlt_time == "weibull" && any(p_resp >= 1)) {
    stop("'p_resp' must be below 1 when 'dlt_time' is \"weibull\"", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         !is.finite(seed) || seed != round(seed))) {
    stop("'seed' must be a single whole number or NULL", call. = FALSE)
  }
  if (n_earlystop <= 6) {
    warning("'n_earlystop' is low; values between 9 and 18 are recommended",
            call. = FALSE)
  }

  n_doses <- length(p_true)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of doses in 'p_true'", call. = FALSE)
  }

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)

  # Pooled estimates can involve every patient of the trial, so the boundaries
  # cover the escalation patients plus a full backfill cap at every dose.
  bound <- boin_boundary(target, max_total_pts + n_cap * n_doses,
                         p_saf = p_saf, p_tox = p_tox, cutoff_eli = cutoff_eli,
                         extrasafe = extrasafe, offset = offset,
                         stay_on_1_of_3 = stay_on_1_of_3)

  # The engine encodes a missing elimination boundary as zero.
  b_elim <- bound$b_elim
  b_elim[is.na(b_elim)] <- 0L

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
    stream_seed <- seed
  } else {
    stream_seed <- floor(stats::runif(1L) * 2147483647)
  }

  res <- backfill_simulate_cpp(
    n_trials = as.integer(n_trials),
    estimate = if (imputation) 1L else 0L,
    p_true = as.numeric(p_true),
    p_resp = as.numeric(p_resp),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    n_earlystop = as.integer(n_earlystop),
    early_stop_simple = identical(n_earlystop_rule, "simple"),
    extrasafe = extrasafe,
    target = as.numeric(target),
    cutoff_eli = as.numeric(cutoff_eli),
    offset = as.numeric(offset),
    b_esc = bound$b_esc,
    b_deesc = bound$b_deesc,
    b_elim = as.integer(b_elim),
    max_total_pts = as.integer(max_total_pts),
    lambda_e = bound$lambda_e,
    lambda_d = bound$lambda_d,
    max_pending_ratio = if (imputation) rules$max_pending_ratio else 1,
    min_completed = if (imputation) rules$min_completed else 0L,
    min_follow_up = if (imputation) rules$min_follow_up else 0,
    window = as.numeric(window),
    accrual = match(accrual, c("exponential", "uniform", "fixed")) - 1L,
    accrual_rate = as.numeric(accrual_rate),
    dlt_time = match(dlt_time, c("weibull", "uniform")) - 1L,
    late_fraction = as.numeric(late_fraction),
    resp_window = as.numeric(resp_window),
    resp_late_fraction = as.numeric(resp_late_fraction),
    resp_cor = as.numeric(resp_cor),
    n_cap = as.integer(n_cap),
    backfill_dose = match(backfill_dose, c("highest", "lowest")) - 1L,
    conflict_dose = match(conflict_dose, c("highest", "lowest")) - 1L,
    no_slot = match(no_slot, c("wait", "leave")) - 1L,
    stream_seed = as.numeric(stream_seed)
  )

  reasons <- c("lowest_dose_eliminated", "lowest_dose_too_toxic",
               "n_earlystop", "max_sample_size", "max_cohorts")

  dose_names <- paste0("DL", seq_len(n_doses))
  for (component in c("n_pts", "n_tox", "n_bf", "n_resp", "eliminated")) {
    colnames(res[[component]]) <- dose_names
  }

  settings <- list(
    design = design,
    target = target,
    p_true = p_true,
    p_resp = p_resp,
    n_cohort = n_cohort,
    cohort_size = cohort_size_vec,
    window = window,
    accrual_rate = accrual_rate,
    n_cap = n_cap,
    backfill_dose = backfill_dose,
    conflict_dose = conflict_dose,
    no_slot = no_slot,
    accrual = accrual,
    dlt_time = dlt_time,
    late_fraction = late_fraction,
    resp_window = resp_window,
    resp_late_fraction = resp_late_fraction,
    resp_cor = resp_cor
  )
  if (imputation) {
    settings$max_pending_ratio <- rules$max_pending_ratio
    settings$min_completed <- rules$min_completed
    settings$min_follow_up <- rules$min_follow_up
  } else {
    settings$stay_on_1_of_3 <- stay_on_1_of_3
  }
  settings <- c(settings, list(
    n_trials = n_trials,
    start_dose = start_dose,
    n_earlystop = n_earlystop,
    n_earlystop_rule = n_earlystop_rule,
    p_saf = bound$p_saf,
    p_tox = bound$p_tox,
    cutoff_eli = cutoff_eli,
    extrasafe = extrasafe,
    offset = offset,
    max_total_pts = max_total_pts,
    seed = seed
  ))

  structure(
    list(
      n_pts = res$n_pts,
      n_tox = res$n_tox,
      n_bf = res$n_bf,
      n_resp = res$n_resp,
      eliminated = res$eliminated,
      cohorts_used = res$cohorts_used,
      stop_reason = reasons[res$stop_code + 1L],
      duration = res$duration,
      n_suspensions = res$n_suspensions,
      time_suspended = res$time_suspended,
      n_turned_away = res$n_turned_away,
      boundary = bound,
      settings = settings
    ),
    class = c(paste0(design, "_trials"), "backfill_trials")
  )
}
