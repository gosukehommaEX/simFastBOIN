#' Operating Characteristics of the BF-BOIN-ET Design
#'
#' @description
#'   Simulate trials of the backfill BOIN-ET (BF-BOIN-ET) design of Takeda,
#'   Zhu and Hirakawa (2025) under one scenario of DLT and response
#'   probabilities, select the optimal biological dose (OBD) of each trial,
#'   and summarize how often each dose is selected, how many patients are
#'   treated and backfilled at each dose, and how long the trials last.
#'
#' @inheritParams tite_boin_et_simulate
#' @inheritParams sim_boin_et
#'
#' @param target
#'   Numeric scalar. Target DLT probability.
#'
#' @param target_resp
#'   Numeric scalar. Lowest response probability of interest. A dose is
#'   excluded when its response rate is likely below it, and the backfill
#'   doses start from the lowest dose whose response rate reaches the
#'   boundary \code{eta1} derived from it.
#'
#' @param n_earlystop
#'   Integer scalar. Number of patients at a dose, dose escalation and
#'   backfill together, that closes it to backfilling, and that ends the trial
#'   when it is reached at the current dose and the decision would keep the
#'   trial there. Defaults to 12, the default of the 'bfboinet' package.
#'
#' @param n_elimination
#'   Integer scalar. Number of patients a dose needs before it can be
#'   excluded. Defaults to 6.
#'
#' @param suspend
#'   Numeric scalar from 0 to 1. Share of the patients at the current dose
#'   whose DLT assessment may still be running when the next cohort starts.
#'   Defaults to 0, which waits until the DLT outcomes of the cohort are
#'   known.
#'
#' @param p_saf
#'   Numeric scalar or \code{NULL}. Highest DLT probability deemed
#'   subtherapeutic, used for the escalation boundary. Defaults to
#'   \code{0.6 * target}.
#'
#' @param p_tox
#'   Numeric scalar or \code{NULL}. Lowest DLT probability deemed overly
#'   toxic, used for the de-escalation boundary. Defaults to
#'   \code{1.4 * target}.
#'
#' @param p_saf_resp
#'   Numeric scalar or \code{NULL}. Response probability deemed
#'   subtherapeutic, used for the response boundary. Defaults to
#'   \code{0.6 * target_resp}.
#'
#' @param cutoff_eli
#'   Numeric scalar. A dose with \code{n_elimination} patients or more is
#'   excluded, together with every higher dose, when the posterior probability
#'   that its DLT rate exceeds \code{target} is above this cutoff. Defaults
#'   to 0.95.
#'
#' @param cutoff_eli_resp
#'   Numeric scalar. Such a dose is excluded alone when the posterior
#'   probability that its response rate is below \code{target_resp} is above
#'   this cutoff. Defaults to 0.90.
#'
#' @param psi00,psi11
#'   Numeric scalars. Scores of a patient without DLT or response and of a
#'   patient with both, which define, with \code{psi_resp_no_dlt} and
#'   \code{psi_dlt_no_resp}, the utility used to choose the backfill dose and
#'   the OBD. Default to 40 and 60.
#'
#' @param psi_resp_no_dlt,psi_dlt_no_resp
#'   Numeric scalars. Scores of a patient with a response and no DLT and of a
#'   patient with a DLT and no response. Default to 100 and 0. The 'bfboinet'
#'   package calls them \code{ppsi10} and \code{ppsi01}.
#'
#' @param backfill_extension
#'   Numeric scalar. Time, in the unit of \code{window}, during which
#'   backfilling continues after the last cohort has been enrolled. Defaults
#'   to 100, the value fixed in the 'bfboinet' package.
#'
#' @return
#'   An object of class \code{bf_boin_et_oc}, which is a list with components
#'   \item{design}{\code{"BF-BOIN-ET"}.}
#'   \item{target, target_resp, p_true, p_resp}{The targets and the scenario.}
#'   \item{n_doses, n_trials}{Numbers of doses and of simulated trials.}
#'   \item{lambda1, lambda2, eta1}{The escalation and de-escalation boundaries
#'     of the DLT rate and the boundary of the response rate.}
#'   \item{sel_percent}{Percentage of trials selecting each dose as the OBD.}
#'   \item{percent_no_obd}{Percentage of trials selecting no OBD.}
#'   \item{percent_early_stop}{Percentage of trials stopped because every dose
#'     was excluded.}
#'   \item{n_pts_dose, n_bf_dose}{Average numbers of patients treated and of
#'     patients backfilled at each dose.}
#'   \item{n_tox_dose, n_resp_dose}{Average numbers of DLTs and responses at
#'     each dose observed at the last decision.}
#'   \item{n_tox_bf_dose, n_resp_bf_dose}{The same among the backfilled
#'     patients.}
#'   \item{total_n_pts, total_n_bf, total_n_tox, total_n_resp}{Their averages
#'     per trial.}
#'   \item{duration}{Average time until the last patient finishes both
#'     assessments.}
#'   \item{stop_reason_percent}{Percentage of trials ending for each reason, see
#'     Details.}
#'   \item{settings}{List of the design parameters.}
#'   \item{call}{The matched call.}
#'   \item{trials}{When \code{keep_trials} is \code{TRUE}, a list with the
#'     matrices \code{n_pts}, \code{n_bf}, \code{n_tox}, \code{n_resp},
#'     \code{n_tox_bf} and \code{n_resp_bf} (one row per trial), and the
#'     vectors \code{cohorts_used}, \code{stop_reason}, \code{early_stop},
#'     \code{obd} (\code{NA} when none) and \code{duration}.}
#'
#' @details
#'   The rules are those of \code{bfboinet::get.oc.backboinet()} (version
#'   1.1.0). A cohort is treated at the current dose. The next cohort starts
#'   with the first patient who arrives once the DLT outcomes of the cohort
#'   are known (or, with \code{suspend}, once enough of them are), and after
#'   the last cohort backfilling goes on for \code{backfill_extension}. The
#'   patients who arrive in between are backfilled one by one, when the
#'   current dose is not the lowest, to the backfill dose with the highest
#'   utility (ties to the higher dose), the DLT rates being made monotone by
#'   isotonic regression; doses with \code{n_earlystop} patients take no more.
#'
#'   The decision uses the DLT rate at the current dose, which counts a
#'   patient still under assessment by the fraction of the window followed,
#'   with the boundaries \code{lambda1} and \code{lambda2} of BOIN; after 1 DLT
#'   in 3 patients the trial stays. The backfill doses for the next cohort are
#'   the admissible doses from the lowest dose whose response rate is at
#'   least \code{eta1} up to the current dose (escalation), the dose below it
#'   (stay) or two doses below it (de-escalation). Once excluded, a dose stays
#'   excluded, and an escalation moves past a dose excluded by its response
#'   data only. The OBD is the admissible dose up to the estimated MTD with
#'   the highest utility, from the observed rates at the last decision; ties go
#'   to the lower dose.
#'
#'   Possible values of \code{stop_reason} are \code{"no_admissible_dose"} (the
#'   only early stop, which selects no OBD when fewer patients than the
#'   maximum sample size have been treated), \code{"lowest_dose_excluded"},
#'   \code{"top_dose_excluded"}, \code{"no_lower_admissible_dose"},
#'   \code{"n_earlystop"} and \code{"max_cohorts"}.
#'
#'   The 'bfboinet' package reports the position of the OBD among the doses
#'   tried, which is the dose itself when the trial starts at the lowest dose;
#'   the dose is reported here in every case.
#'
#' @references
#'   Takeda, K., Zhu, J. and Hirakawa, A. (2025). BF-BOIN-ET: A backfill
#'   Bayesian optimal interval design using efficacy and toxicity outcomes for
#'   dose optimization. Pharmaceutical Statistics, 24, e2470.
#'   \doi{10.1002/pst.2470}
#'
#' @examples
#' # The example of bfboinet::get.oc.backboinet(): windows of one month and
#' # three patients a month
#' oc <- sim_bf_boin_et(
#'   target = 0.30, target_resp = 0.25,
#'   p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
#'   p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70),
#'   n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
#'   n_trials = 300, seed = 123
#' )
#' oc
#'
#' @seealso \code{\link{sim_bf_boin_et_multi}}, \code{\link{sim_tite_boin_et}}
#' @export
sim_bf_boin_et <- function(target, target_resp, p_true, p_resp, n_cohort,
                           cohort_size, window, accrual_rate,
                           resp_window = window,
                           accrual = c("uniform", "exponential", "fixed"),
                           event_time = c("weibull", "uniform"),
                           late_fraction = 0.5, resp_late_fraction = 0.5,
                           resp_cor = 0.2, n_trials = 10000, start_dose = 1,
                           n_earlystop = 12, n_elimination = 6, suspend = 0,
                           p_saf = NULL, p_tox = NULL, p_saf_resp = NULL,
                           cutoff_eli = 0.95, cutoff_eli_resp = 0.90,
                           psi00 = 40, psi11 = 60, psi_resp_no_dlt = 100,
                           psi_dlt_no_resp = 0,
                           backfill_extension = 100, keep_trials = FALSE,
                           verbose = FALSE, seed = 123) {

  call_expr <- match.call()
  accrual <- match.arg(accrual)
  event_time <- match.arg(event_time)
  check_flag(keep_trials, "keep_trials")
  check_flag(verbose, "verbose")

  check_p_true(p_true)
  n_doses <- length(p_true)
  if (!is.numeric(p_resp) || length(p_resp) != n_doses ||
      any(!is.finite(p_resp)) || any(p_resp < 0) || any(p_resp > 1)) {
    stop("'p_resp' must hold one response probability between 0 and 1 ",
         "per dose in 'p_true'", call. = FALSE)
  }
  check_scalar_prob(target, "target")
  check_scalar_prob(target_resp, "target_resp")
  if (is.null(p_saf)) p_saf <- 0.6 * target
  if (is.null(p_tox)) p_tox <- 1.4 * target
  if (is.null(p_saf_resp)) p_saf_resp <- 0.6 * target_resp
  check_scalar_prob(p_saf, "p_saf")
  check_scalar_prob(p_tox, "p_tox")
  check_scalar_prob(p_saf_resp, "p_saf_resp")
  if (!(p_saf < target && target < p_tox)) {
    stop("'p_saf' < 'target' < 'p_tox' must hold", call. = FALSE)
  }
  if (!(p_saf_resp < target_resp)) {
    stop("'p_saf_resp' must be below 'target_resp'", call. = FALSE)
  }
  check_count(n_cohort, "n_cohort", 1L)
  check_count(n_trials, "n_trials", 1L)
  check_count(start_dose, "start_dose", 1L)
  if (start_dose > n_doses) {
    stop("'start_dose' must not exceed the number of doses in 'p_true'",
         call. = FALSE)
  }
  check_count(n_earlystop, "n_earlystop", 1L)
  check_count(n_elimination, "n_elimination", 1L)
  if (!is.numeric(suspend) || length(suspend) != 1L || !is.finite(suspend) ||
      suspend < 0 || suspend > 1) {
    stop("'suspend' must be a single number from 0 to 1", call. = FALSE)
  }
  if (!is.numeric(backfill_extension) || length(backfill_extension) != 1L ||
      !is.finite(backfill_extension) || backfill_extension < 0) {
    stop("'backfill_extension' must be a single non-negative number",
         call. = FALSE)
  }
  check_scalar_prob(cutoff_eli, "cutoff_eli")
  check_scalar_prob(cutoff_eli_resp, "cutoff_eli_resp")
  check_scalar_prob(resp_cor, "resp_cor", lower = -1, upper = 1)
  for (nm in c("psi00", "psi11", "psi_resp_no_dlt", "psi_dlt_no_resp")) {
    value <- get(nm)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop("'", nm, "' must be a single finite number", call. = FALSE)
    }
  }
  boin_et_check_timing(window, resp_window, accrual_rate, late_fraction,
                       resp_late_fraction)
  par_t <- boin_et_event_params(p_true, window, late_fraction, event_time,
                                "p_true")
  par_e <- boin_et_event_params(p_resp, resp_window, resp_late_fraction,
                                event_time, "p_resp")

  cohort_size_vec <- expand_cohort_size(cohort_size, n_cohort)
  max_total_pts <- sum(cohort_size_vec)

  # Boundaries of BOIN for the DLT rate and the same formula for the
  # response rate, written as in bfboinet::get.oc.backboinet().
  lambda1 <- log((1 - p_saf) / (1 - target)) /
    log(target * (1 - p_saf) / (p_saf * (1 - target)))
  lambda2 <- log((1 - target) / (1 - p_tox)) /
    log(p_tox * (1 - target) / (target * (1 - p_tox)))
  eta1 <- log((1 - p_saf_resp) / (1 - target_resp)) /
    log(target_resp * (1 - p_saf_resp) / (p_saf_resp * (1 - target_resp)))

  max_n <- max_total_pts + n_earlystop
  adm <- boin_et_admissible_table(target, target_resp, cutoff_eli,
                                  cutoff_eli_resp, max_n)

  if (!is.null(seed)) {
    old_seed <- get_random_seed()
    on.exit(restore_random_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  if (verbose) message("Simulating ", n_trials, " trials ...")

  res <- bf_boin_et_simulate_cpp(
    n_trials = as.integer(n_trials),
    p_true = as.numeric(p_true),
    p_resp = as.numeric(p_resp),
    resp_cor = as.numeric(resp_cor),
    cohort_size = cohort_size_vec,
    start_dose = as.integer(start_dose),
    lambda1 = lambda1,
    lambda2 = lambda2,
    eta1 = eta1,
    n_elimination = as.integer(n_elimination),
    n_earlystop = as.integer(n_earlystop),
    suspend = as.numeric(suspend),
    extension = as.numeric(backfill_extension),
    tox_ok = as.integer(adm$tox_ok),
    resp_ok = as.integer(adm$resp_ok),
    max_n = as.integer(max_n),
    target = as.numeric(target),
    psi00 = as.numeric(psi00),
    psi01 = as.numeric(psi_dlt_no_resp),
    psi10 = as.numeric(psi_resp_no_dlt),
    psi11 = as.numeric(psi11),
    window = as.numeric(window),
    resp_window = as.numeric(resp_window),
    accrual_rate = as.numeric(accrual_rate),
    arrival = match(accrual, c("uniform", "exponential", "fixed")) - 1L,
    event_time = match(event_time, c("weibull", "uniform")) - 1L,
    shape_t = par_t$shape,
    scale_t = par_t$scale,
    shape_e = par_e$shape,
    scale_e = par_e$scale
  )

  dose_names <- paste0("DL", seq_len(n_doses))
  for (nm in c("n_pts", "n_bf", "n_tox", "n_resp", "n_tox_bf", "n_resp_bf")) {
    colnames(res[[nm]]) <- dose_names
  }
  reasons <- c("no_admissible_dose", "lowest_dose_excluded",
               "top_dose_excluded", "no_lower_admissible_dose",
               "n_earlystop", "max_cohorts")
  stop_reason <- reasons[res$stop_code + 1L]
  obd <- res$obd
  obd[obd == 0L] <- NA_integer_

  means <- function(m) {
    out <- colMeans(m)
    names(out) <- dose_names
    out
  }
  selected <- factor(obd, levels = seq_len(n_doses))
  sel_percent <- as.numeric(table(selected)) / n_trials * 100
  names(sel_percent) <- dose_names
  stop_reasons <- as.numeric(table(factor(stop_reason, levels = reasons))) /
    n_trials * 100
  names(stop_reasons) <- reasons

  n_pts_dose <- means(res$n_pts)
  n_bf_dose <- means(res$n_bf)
  n_tox_dose <- means(res$n_tox)
  n_resp_dose <- means(res$n_resp)

  out <- list(
    design = "BF-BOIN-ET",
    target = target,
    target_resp = target_resp,
    p_true = p_true,
    p_resp = p_resp,
    n_doses = n_doses,
    n_trials = n_trials,
    lambda1 = lambda1,
    lambda2 = lambda2,
    eta1 = eta1,
    sel_percent = sel_percent,
    percent_no_obd = mean(is.na(obd)) * 100,
    percent_early_stop = mean(res$early_stop) * 100,
    n_pts_dose = n_pts_dose,
    n_bf_dose = n_bf_dose,
    n_tox_dose = n_tox_dose,
    n_resp_dose = n_resp_dose,
    n_tox_bf_dose = means(res$n_tox_bf),
    n_resp_bf_dose = means(res$n_resp_bf),
    total_n_pts = sum(n_pts_dose),
    total_n_bf = sum(n_bf_dose),
    total_n_tox = sum(n_tox_dose),
    total_n_resp = sum(n_resp_dose),
    duration = mean(res$duration),
    stop_reason_percent = stop_reasons,
    settings = list(
      target = target, target_resp = target_resp, p_true = p_true,
      p_resp = p_resp, n_cohort = n_cohort, cohort_size = cohort_size_vec,
      window = window, resp_window = resp_window,
      accrual_rate = accrual_rate, accrual = accrual,
      event_time = event_time, late_fraction = late_fraction,
      resp_late_fraction = resp_late_fraction, resp_cor = resp_cor,
      n_trials = n_trials, start_dose = start_dose,
      n_earlystop = n_earlystop, n_elimination = n_elimination,
      suspend = suspend, p_saf = p_saf, p_tox = p_tox,
      p_saf_resp = p_saf_resp, cutoff_eli = cutoff_eli,
      cutoff_eli_resp = cutoff_eli_resp, psi00 = psi00, psi11 = psi11,
      psi_resp_no_dlt = psi_resp_no_dlt, psi_dlt_no_resp = psi_dlt_no_resp,
      backfill_extension = backfill_extension,
      max_total_pts = max_total_pts, seed = seed
    ),
    call = call_expr
  )
  if (keep_trials) {
    out$trials <- list(
      n_pts = res$n_pts, n_bf = res$n_bf, n_tox = res$n_tox,
      n_resp = res$n_resp, n_tox_bf = res$n_tox_bf,
      n_resp_bf = res$n_resp_bf, cohorts_used = res$cohorts_used,
      stop_reason = stop_reason, early_stop = as.logical(res$early_stop),
      obd = obd, duration = res$duration
    )
  }

  structure(out, class = "bf_boin_et_oc")
}
