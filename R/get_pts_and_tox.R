#' Generate Patient and Toxicity Data for BOIN Simulations
#'
#' @description
#'   Conduct vectorized BOIN trial simulations with optional titration phase
#'   to generate patient enrollment and toxicity (DLT) data across multiple trials.
#'   This function focuses on Step 1 of the operating characteristic evaluation workflow:
#'   generating n_pts, n_tox, and eliminated_mat matrices. All computations are fully
#'   vectorized using matrix operations. Titration uses parallel Bernoulli trials across
#'   all doses simultaneously, following the exact BOIN algorithm.
#'
#' @param n_trials
#'   Numeric. Number of trials to simulate. Default is 10000.
#'   Typically 1000-10000 for detailed operating characteristics.
#'
#' @param target
#'   Numeric. The target toxicity probability (e.g., 0.30 for 30%).
#'
#' @param p_true
#'   Numeric vector. True toxicity probabilities for each dose.
#'   Length determines number of doses evaluated. Must be increasing vector.
#'
#' @param n_cohort
#'   Numeric. Maximum number of cohorts per trial.
#'
#' @param cohort_size
#'   Numeric vector or scalar specifying patients per cohort.
#'   If vector (e.g., c(3, 3, 3)), each element specifies size for corresponding cohort.
#'   If scalar, all cohorts use the same size.
#'
#' @param n_earlystop
#'   Numeric. Sample size at current dose triggering trial termination.
#'   Default is 18. Note that the actual maximum sample size per dose may exceed
#'   this value when using n_earlystop_rule = "with_stay".
#'
#' @param cutoff_eli
#'   Numeric. Cutoff probability for dose elimination. Default is 0.95.
#'   If Pr(toxicity > target | data) > cutoff_eli, eliminate the dose and higher doses.
#'
#' @param extrasafe
#'   Logical. If TRUE, apply additional safety stopping rule at the lowest dose.
#'   Default is FALSE. When TRUE, the trial stops early if the lowest dose is
#'   overly toxic based on cutoff_eli - offset.
#'
#' @param offset
#'   Numeric. Offset for safety stopping boundary when extrasafe = TRUE.
#'   Default is 0.05. The safety cutoff becomes cutoff_eli - offset.
#'
#' @param n_earlystop_rule
#'   Character. Rule for early stopping at n_earlystop. Options are "simple"
#'   (stop when n >= n_earlystop) or "with_stay" (stop when n >= n_earlystop
#'   AND decision = "Stay"). Default is "with_stay". "with_stay" follows the
#'   BOIN standard implementation.
#'
#' @param titration
#'   Logical. If TRUE, perform accelerated dose escalation with parallel Bernoulli
#'   trials at all doses simultaneously, then add (cohort_size - 1) patients to the
#'   dose with first toxicity. Default is FALSE.
#'
#' @param seed
#'   Numeric. Random seed for reproducibility. Default is 123.
#'
#' @return
#'   A list containing:
#'   \item{n_pts_all}{Matrix of size n_trials x n_doses: patient counts at each dose}
#'   \item{n_tox_all}{Matrix of size n_trials x n_doses: DLT counts at each dose}
#'   \item{eliminated_mat}{Matrix of size n_trials x n_doses: logical, TRUE if dose eliminated}
#'   \item{cohorts_completed}{Integer vector of length n_trials: cohorts completed per trial}
#'   \item{stop_reason}{Character vector of length n_trials: reason for stopping}
#'
#' @details
#'   **Decision Table Size Calculation:**
#'
#'   When n_earlystop_rule = "with_stay", a single dose can exceed n_earlystop
#'   because the trial only stops when both n >= n_earlystop AND decision = "Stay".
#'   The theoretical maximum sample size per dose is n_cohort * max(cohort_size).
#'   Decision tables and stopping boundaries are generated using this value to
#'   ensure all possible patient counts are covered.
#'
#'   **Titration Phase Logic (exactly matching BOIN):**
#'   When titration = TRUE:
#'   \enumerate{
#'     \item Generate Bernoulli trials for all doses simultaneously
#'     \item Find first dose d with DLT
#'     \item Allocate 1 patient to each dose from 1 to d
#'     \item Add (cohort_size - 1) patients to dose d
#'     \item Execute dose decision (escalate/deescalate/eliminate)
#'     \item Exit titration phase and proceed with normal cohort_size
#'   }
#'
#'   **Key Performance Features:**
#'   \itemize{
#'     \item Automatic generation of decision tables and boundaries
#'     \item Optimized random number generation using runif() for DLT generation
#'     \item Matrix-based decision table lookups via vectorized indexing
#'     \item Early termination of inactive trials to reduce overhead
#'     \item Full vectorization: all trials processed simultaneously per cohort
#'   }
#'
#'   **DLT Generation Optimization:**
#'   DLT outcomes are generated using runif() for significant performance gains.
#'   For single-patient cohorts: dlt <- as.integer(runif(n) < p_true).
#'   For multi-patient cohorts: generate matrix and compare row-wise.
#'   This approach is faster than rbinom() because runif() is computationally
#'   simpler, enables fully vectorized threshold comparison, and reduces
#'   function call overhead.
#'
#'   **Stopping Reasons:**
#'   Trials stop for one of these reasons:
#'   - "n_earlystop_simple": Sample size reached (simple rule)
#'   - "n_earlystop_with_stay": Sample size reached with Stay decision (BOIN standard)
#'   - "lowest_dose_too_toxic": Safety stopping at lowest dose
#'   - "lowest_dose_eliminated": Lowest dose eliminated (cannot continue)
#'   - "max_cohorts_reached": Maximum number of cohorts completed
#'   - "max_sample_size_reached": Maximum total sample size exceeded
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#'   Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020).BOIN: An R Package
#'   for Designing Single-Agent and Drug-Combination Dose-Finding Trials Using Bayesian
#'   Optimal Interval Designs. Journal of Statistical Software, 94(13), 1-32.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic BOIN simulation without titration
#' target <- 0.30
#' p_true <- c(0.10, 0.25, 0.40, 0.55, 0.70)
#'
#' result <- get_pts_and_tox(
#'   n_trials = 1000,
#'   target = target,
#'   p_true = p_true,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   n_earlystop = 18,
#'   seed = 123
#' )
#'
#' # Extract results
#' n_pts <- result$n_pts_all
#' n_tox <- result$n_tox_all
#' eliminated <- result$eliminated_mat
#'
#' # Check operating characteristics
#' colMeans(n_pts)  # Average patients at each dose
#' colMeans(n_tox)  # Average DLTs at each dose
#' table(result$stop_reason)  # Distribution of stopping reasons
#'
#' # Example 2: BOIN with titration
#' result_tit <- get_pts_and_tox(
#'   n_trials = 1000,
#'   target = target,
#'   p_true = p_true,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   titration = TRUE,
#'   seed = 123
#' )
#'
#' # Compare titration vs non-titration
#' colMeans(result_tit$n_pts_all)  # Fewer patients at lower doses
#'
#' # Example 3: BOIN with extra safety
#' result_safe <- get_pts_and_tox(
#'   n_trials = 1000,
#'   target = target,
#'   p_true = p_true,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   extrasafe = TRUE,
#'   offset = 0.05,
#'   seed = 123
#' )
#'
#' table(result_safe$stop_reason)  # More safety stops
#' }
#'
#' @importFrom stats runif rbinom
#' @export
get_pts_and_tox <- function(
    n_trials = 10000,
    target,
    p_true,
    n_cohort,
    cohort_size,
    n_earlystop = 18,
    cutoff_eli = 0.95,
    extrasafe = FALSE,
    offset = 0.05,
    n_earlystop_rule = c("with_stay", "simple"),
    titration = FALSE,
    seed = 123
) {

  n_earlystop_rule <- match.arg(n_earlystop_rule)
  set.seed(seed)

  n_doses <- length(p_true)

  # Calculate max_sample_size for decision tables
  if (length(cohort_size) == 1) {
    max_cohort_size <- cohort_size
  } else {
    max_cohort_size <- max(cohort_size)
  }

  max_sample_size_for_tables <- n_cohort * max_cohort_size

  # Generate decision table
  boin_bound <- get_boin_boundary(target)
  decision_table <- get_boin_decision(
    target = target,
    lambda_e = boin_bound$lambda_e,
    lambda_d = boin_bound$lambda_d,
    max_sample_size = max_sample_size_for_tables,
    cutoff_eli = cutoff_eli
  )

  if (extrasafe) {
    cutoff_stop <- cutoff_eli - offset
    stopping_boundaries <- get_boin_stopping_boundaries(
      target = target,
      max_sample_size = max_sample_size_for_tables,
      cutoff_stop = cutoff_stop
    )
  } else {
    stopping_boundaries <- NULL
  }

  # Initialize matrices
  n_pts_mat <- matrix(0L, nrow = n_trials, ncol = n_doses)
  n_tox_mat <- matrix(0L, nrow = n_trials, ncol = n_doses)
  current_dose_vec <- rep(1L, n_trials)
  eliminated_mat <- matrix(FALSE, nrow = n_trials, ncol = n_doses)
  active_trials <- rep(TRUE, n_trials)
  cohorts_completed <- rep(0L, n_trials)
  stop_reason <- rep(NA_character_, n_trials)

  # Titration tracking
  in_titration <- rep(titration, n_trials)
  ft_flag <- rep(titration, n_trials)

  # Pre-compute cohort size vector
  if (length(cohort_size) == 1) {
    cohort_size_vec <- rep(cohort_size, n_cohort)
  } else if (length(cohort_size) < n_cohort) {
    cohort_size_vec <- c(
      cohort_size,
      rep(cohort_size[length(cohort_size)], n_cohort - length(cohort_size))
    )
  } else {
    cohort_size_vec <- cohort_size[1:n_cohort]
  }

  max_total_pts <- n_cohort * cohort_size_vec[1]
  max_col_decision <- ncol(decision_table)
  if (extrasafe) {
    max_col_stopping <- ncol(stopping_boundaries)
  }

  # TITRATION INITIALIZATION (Parallel Bernoulli)
  if (titration) {
    z_mat <- matrix(runif(n_trials * n_doses), nrow = n_trials, ncol = n_doses) <
      matrix(rep(p_true, n_trials), nrow = n_trials, byrow = TRUE)

    first_dlt_dose <- apply(z_mat, 1, function(row) {
      idx <- which(row)[1]
      if (is.na(idx)) n_doses + 1 else idx
    })

    for (trial in 1:n_trials) {
      if (first_dlt_dose[trial] > n_doses) {
        n_pts_mat[trial, 1:n_doses] <- 1L
        current_dose_vec[trial] <- n_doses
      } else {
        d <- first_dlt_dose[trial]
        n_pts_mat[trial, 1:d] <- 1L
        n_tox_mat[trial, d] <- 1L
        current_dose_vec[trial] <- d
      }
    }

    cohorts_completed[1:n_trials] <- 1L

    if (cohort_size > 1) {
      fillup_count <- cohort_size - 1L

      for (trial in 1:n_trials) {
        if (ft_flag[trial]) {
          d <- current_dose_vec[trial]

          if (fillup_count == 1) {
            u <- runif(1)
            dlt_fillup <- as.integer(u < p_true[d])
          } else {
            u_fillup <- runif(fillup_count)
            dlt_fillup <- sum(u_fillup < p_true[d])
          }

          n_pts_mat[trial, d] <- n_pts_mat[trial, d] + fillup_count
          n_tox_mat[trial, d] <- n_tox_mat[trial, d] + dlt_fillup

          ft_flag[trial] <- FALSE
          in_titration[trial] <- FALSE
        }
      }
    } else {
      in_titration[1:n_trials] <- FALSE
    }
  }

  # FIRST DOSE DECISION (after titration initialization)
  if (titration) {
    for (trial in 1:n_trials) {
      if (!active_trials[trial]) next

      d <- current_dose_vec[trial]
      n_pts_d <- n_pts_mat[trial, d]
      n_tox_d <- n_tox_mat[trial, d]

      if (!is.na(n_pts_d)) {
        n_val_eli <- pmin(n_pts_d, max_col_decision)
        y_val_eli <- pmin(n_tox_d, n_val_eli)

        if (!is.na(decision_table[y_val_eli + 1L, n_val_eli])) {
          elim_decision <- decision_table[y_val_eli + 1L, n_val_eli]

          if (elim_decision == "DE") {
            eliminated_mat[trial, d:n_doses] <- TRUE
            if (d > 1) {
              new_d <- d - 1L
              while (new_d > 1 && eliminated_mat[trial, new_d]) {
                new_d <- new_d - 1L
              }
              current_dose_vec[trial] <- new_d
            } else {
              active_trials[trial] <- FALSE
              stop_reason[trial] <- "lowest_dose_eliminated"
            }
          } else if (elim_decision == "E" && d < n_doses && !eliminated_mat[trial, d + 1L]) {
            current_dose_vec[trial] <- d + 1L
          } else if (elim_decision == "D" && d > 1) {
            new_d <- d - 1L
            while (new_d > 1 && eliminated_mat[trial, new_d]) {
              new_d <- new_d - 1L
            }
            current_dose_vec[trial] <- new_d
          }
        }
      }

      if (extrasafe && active_trials[trial]) {
        n_pts_dose1 <- n_pts_mat[trial, 1]
        n_tox_dose1 <- n_tox_mat[trial, 1]

        if (n_pts_dose1 >= 3) {
          n_pts_for_stop <- pmin(n_pts_dose1, max_col_stopping)
          stop_decision <- stopping_boundaries[n_tox_dose1 + 1L, n_pts_for_stop]

          if (!is.na(stop_decision) && stop_decision == "STOP") {
            active_trials[trial] <- FALSE
            stop_reason[trial] <- "lowest_dose_too_toxic"
          }
        }
      }
    }
  }

  # MAIN COHORT LOOP (Normal phase)
  start_cohort <- if (titration) 2L else 1L
  for (cohort in start_cohort:n_cohort) {

    normal_active_idx <- which(active_trials & !in_titration)
    n_normal <- length(normal_active_idx)

    if (n_normal == 0) break

    current_cohort_size <- cohort_size_vec[cohort]

    # Early stopping check (simple rule)
    if (n_earlystop_rule == "simple") {
      early_stop_check <- n_pts_mat[cbind(normal_active_idx, current_dose_vec[normal_active_idx])] >= n_earlystop
      if (any(early_stop_check)) {
        early_stop_trials <- normal_active_idx[early_stop_check]
        active_trials[early_stop_trials] <- FALSE
        cohorts_completed[early_stop_trials] <- cohort
        stop_reason[early_stop_trials] <- "n_earlystop_simple"

        normal_active_idx <- normal_active_idx[!early_stop_check]
        n_normal <- length(normal_active_idx)
        if (n_normal == 0) next
      }
    }

    # Check max total patients
    current_total_pts <- rowSums(n_pts_mat[normal_active_idx, , drop = FALSE])
    will_exceed <- (current_total_pts + current_cohort_size) > max_total_pts

    if (any(will_exceed)) {
      exceed_idx <- which(will_exceed)
      exceed_trials <- normal_active_idx[exceed_idx]
      remaining_pts <- max_total_pts - current_total_pts[exceed_idx]

      for (i in seq_along(exceed_idx)) {
        if (remaining_pts[i] > 0) {
          trial <- exceed_trials[i]
          dose <- current_dose_vec[trial]
          dlt_count <- rbinom(1, remaining_pts[i], p_true[dose])
          n_pts_mat[trial, dose] <- n_pts_mat[trial, dose] + remaining_pts[i]
          n_tox_mat[trial, dose] <- n_tox_mat[trial, dose] + dlt_count
        }
      }

      active_trials[exceed_trials] <- FALSE
      stop_reason[exceed_trials] <- "max_sample_size_reached"
      cohorts_completed[exceed_trials] <- cohort

      normal_active_idx <- normal_active_idx[-exceed_idx]
      n_normal <- length(normal_active_idx)
      if (n_normal == 0) next
    }

    # Generate DLT
    if (n_normal > 0) {
      current_doses <- current_dose_vec[normal_active_idx]
      p_true_current <- p_true[current_doses]

      if (current_cohort_size == 1) {
        u <- runif(n_normal)
        dlt_counts <- as.integer(u < p_true_current)
      } else {
        u_mat <- matrix(runif(n_normal * current_cohort_size),
                        nrow = n_normal, ncol = current_cohort_size)
        dlt_counts <- rowSums(u_mat < p_true_current)
      }

      update_idx <- cbind(normal_active_idx, current_doses)
      n_pts_mat[update_idx] <- n_pts_mat[update_idx] + current_cohort_size
      n_tox_mat[update_idx] <- n_tox_mat[update_idx] + dlt_counts
    } else {
      next
    }

    cohorts_completed[normal_active_idx] <- cohort

    # Safety stopping at lowest dose
    if (extrasafe) {
      dose1_trials <- normal_active_idx[current_doses == 1]
      if (length(dose1_trials) > 0) {
        n_pts_dose1 <- n_pts_mat[cbind(dose1_trials, rep(1, length(dose1_trials)))]
        n_tox_dose1 <- n_tox_mat[cbind(dose1_trials, rep(1, length(dose1_trials)))]

        check_stop <- n_pts_dose1 >= 3
        if (any(check_stop)) {
          check_idx <- dose1_trials[check_stop]
          n_pts_for_stop <- pmin(n_pts_dose1[check_stop], max_col_stopping)
          stop_decisions <- stopping_boundaries[cbind(n_tox_dose1[check_stop] + 1L, n_pts_for_stop)]

          stopped <- !is.na(stop_decisions) & stop_decisions == "STOP"
          if (any(stopped)) {
            stopped_trials <- check_idx[stopped]
            active_trials[stopped_trials] <- FALSE
            cohorts_completed[stopped_trials] <- cohort
            stop_reason[stopped_trials] <- "lowest_dose_too_toxic"

            normal_active_idx <- setdiff(normal_active_idx, stopped_trials)
            n_normal <- length(normal_active_idx)
            if (n_normal == 0) next
          }
        }
      }
    }

    # Dose decision
    current_doses <- current_dose_vec[normal_active_idx]
    n_pts_current <- n_pts_mat[cbind(normal_active_idx, current_doses)]
    n_tox_current <- n_tox_mat[cbind(normal_active_idx, current_doses)]

    n_vals <- pmin(n_pts_current, max_col_decision)
    y_vals <- pmin(n_tox_current, n_vals)

    decisions <- decision_table[cbind(y_vals + 1L, n_vals)]

    # Process Escalate
    esc_idx <- which(decisions == "E")
    if (length(esc_idx) > 0) {
      esc_trials <- normal_active_idx[esc_idx]
      esc_doses <- current_doses[esc_idx]
      not_at_max <- esc_doses < n_doses

      can_escalate <- rep(FALSE, length(esc_doses))
      if (any(not_at_max)) {
        not_at_max_idx <- which(not_at_max)
        next_dose_not_elim <- !eliminated_mat[cbind(
          esc_trials[not_at_max_idx],
          esc_doses[not_at_max_idx] + 1L
        )]
        can_escalate[not_at_max_idx] <- next_dose_not_elim
      }

      if (any(can_escalate)) {
        current_dose_vec[esc_trials[can_escalate]] <- esc_doses[can_escalate] + 1L
      }
    }

    # Process De-escalate
    deesc_idx <- which(decisions %in% c("D", "DE"))
    if (length(deesc_idx) > 0) {
      deesc_trials <- normal_active_idx[deesc_idx]
      deesc_doses <- current_doses[deesc_idx]

      elim_idx <- which(decisions[deesc_idx] == "DE")
      if (length(elim_idx) > 0) {
        elim_trials <- deesc_trials[elim_idx]
        elim_doses <- deesc_doses[elim_idx]

        for (i in seq_along(elim_idx)) {
          trial_idx <- elim_trials[i]
          dose <- elim_doses[i]

          eliminated_mat[trial_idx, dose:n_doses] <- TRUE

          if (dose > 1) {
            new_dose <- dose - 1L
            while (new_dose > 1 && eliminated_mat[trial_idx, new_dose]) {
              new_dose <- new_dose - 1L
            }
            current_dose_vec[trial_idx] <- new_dose
          } else {
            active_trials[trial_idx] <- FALSE
            cohorts_completed[trial_idx] <- cohort
            stop_reason[trial_idx] <- "lowest_dose_eliminated"
          }
        }

        normal_active_idx <- which(active_trials & !in_titration)
        n_normal <- length(normal_active_idx)
        if (n_normal == 0) next
      }

      no_elim_idx <- which(decisions[deesc_idx] == "D")
      if (length(no_elim_idx) > 0) {
        no_elim_trials <- deesc_trials[no_elim_idx]
        no_elim_doses <- deesc_doses[no_elim_idx]

        can_deescalate <- no_elim_doses > 1
        if (any(can_deescalate)) {
          de_trials <- no_elim_trials[can_deescalate]
          de_doses <- no_elim_doses[can_deescalate]

          for (i in seq_along(de_trials)) {
            trial <- de_trials[i]
            new_dose <- de_doses[i] - 1L
            while (new_dose > 1 && eliminated_mat[trial, new_dose]) {
              new_dose <- new_dose - 1L
            }
            current_dose_vec[trial] <- new_dose
          }
        }
      }
    }

    # Check with_stay early stopping
    if (n_earlystop_rule == "with_stay") {
      current_doses <- current_dose_vec[normal_active_idx]
      current_n_pts <- n_pts_mat[cbind(normal_active_idx, current_doses)]
      current_n_tox <- n_tox_mat[cbind(normal_active_idx, current_doses)]

      n_check <- current_n_pts >= n_earlystop

      if (any(n_check)) {
        check_idx <- which(n_check)
        n_vals <- pmin(current_n_pts[check_idx], max_col_decision)
        y_vals <- pmin(current_n_tox[check_idx], n_vals)

        decisions_check <- decision_table[cbind(y_vals + 1L, n_vals)]
        dose_check <- current_doses[check_idx]

        is_stay <- (decisions_check == "S") |
          (dose_check == 1 & decisions_check %in% c("D", "DE")) |
          (dose_check == n_doses & decisions_check == "E")

        if (any(!is_stay & dose_check < n_doses & decisions_check == "E")) {
          elim_check_idx <- which(!is_stay & dose_check < n_doses & decisions_check == "E")
          for (i in elim_check_idx) {
            trial <- normal_active_idx[check_idx[i]]
            next_dose <- dose_check[i] + 1
            if (eliminated_mat[trial, next_dose]) {
              is_stay[i] <- TRUE
            }
          }
        }

        if (any(is_stay)) {
          stop_trials <- normal_active_idx[check_idx[is_stay]]
          active_trials[stop_trials] <- FALSE
          cohorts_completed[stop_trials] <- cohort
          stop_reason[stop_trials] <- "n_earlystop_with_stay"

          normal_active_idx <- setdiff(normal_active_idx, stop_trials)
          n_normal <- length(normal_active_idx)
          if (n_normal == 0) next
        }
      }
    }
  }

  # Mark remaining active trials
  still_active <- which(active_trials)
  if (length(still_active) > 0) {
    active_trials[still_active] <- FALSE
    stop_reason[still_active] <- "max_cohorts_reached"
  }

  return(list(
    n_pts_all = n_pts_mat,
    n_tox_all = n_tox_mat,
    eliminated_mat = eliminated_mat,
    cohorts_completed = cohorts_completed,
    stop_reason = stop_reason
  ))
}
