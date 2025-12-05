#' Generate Patient and Toxicity Data for BOIN Simulations
#'
#' @description
#'   Conduct vectorized BOIN trial simulations with optional titration phase
#'   to generate patient enrollment and toxicity (DLT) data across multiple trials.
#'
#' @param n_trials
#'   Numeric. Number of trials to simulate. Default is 10000.
#'
#' @param target
#'   Numeric. The target toxicity probability (e.g., 0.30 for 30%).
#'
#' @param p_true
#'   Numeric vector. True toxicity probabilities for each dose.
#'
#' @param n_cohort
#'   Numeric. Maximum number of cohorts per trial.
#'
#' @param cohort_size
#'   Numeric vector or scalar specifying patients per cohort.
#'
#' @param n_earlystop
#'   Numeric. Sample size at current dose triggering trial termination. Default is 18.
#'
#' @param cutoff_eli
#'   Numeric. Cutoff probability for dose elimination. Default is 0.95.
#'
#' @param extrasafe
#'   Logical. If TRUE, apply additional safety stopping rule at the lowest dose. Default is FALSE.
#'
#' @param offset
#'   Numeric. Offset for safety stopping boundary when extrasafe = TRUE. Default is 0.05.
#'
#' @param n_earlystop_rule
#'   Character. Rule for early stopping at n_earlystop. Options are "with_stay" or "simple".
#'   Default is "with_stay".
#'
#' @param titration
#'   Logical. If TRUE, perform accelerated dose escalation in titration phase. Default is FALSE.
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
#'   The function simulates multiple BOIN trials in parallel using vectorized operations.
#'   When titration = TRUE, an accelerated dose escalation phase is performed using
#'   single-patient cohorts to rapidly locate the dose-toxicity region. The main trial
#'   then proceeds with standard cohort-based dosing.
#'
#'   Early stopping rules include:
#'   \itemize{
#'     \item Simple rule: Stop when n_earlystop patients are treated at current dose
#'     \item With-stay rule: Stop when n_earlystop patients are treated and dose stays
#'     \item Extra-safety rule (if extrasafe = TRUE): Stop if lowest dose is excessively toxic
#'   }
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' # Basic BOIN simulation
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
#' # With titration phase
#' result_titration <- get_pts_and_tox(
#'   n_trials = 1000,
#'   target = target,
#'   p_true = p_true,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   titration = TRUE,
#'   seed = 123
#' )
#'
#' # With extra safety stopping
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
#' }
#'
#' @importFrom stats runif rbinom
#'
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
  # Validate and set the n_earlystop_rule parameter
  n_earlystop_rule <- match.arg(n_earlystop_rule)
  set.seed(seed)

  # Get the number of doses from the length of p_true vector
  n_doses <- length(p_true)

  # Determine the maximum cohort size (handle both scalar and vector inputs)
  if (length(cohort_size) == 1) {
    max_cohort_size <- cohort_size
  } else {
    max_cohort_size <- max(cohort_size)
  }

  # Calculate maximum sample size for generating decision tables
  max_sample_size_for_tables <- n_cohort * max_cohort_size

  # Generate BOIN decision boundaries (lambda_e for escalation, lambda_d for de-escalation)
  boin_bound <- get_boin_boundary(target)

  # Generate decision table based on BOIN boundaries and sample size
  decision_table <- get_boin_decision(
    target = target,
    lambda_e = boin_bound$lambda_e,
    lambda_d = boin_bound$lambda_d,
    max_sample_size = max_sample_size_for_tables,
    cutoff_eli = cutoff_eli
  )

  # Generate additional stopping boundaries if extrasafe option is enabled
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

  # Initialize matrices to track patient counts and DLT counts across all trials
  n_pts_mat <- matrix(0L, nrow = n_trials, ncol = n_doses)
  n_tox_mat <- matrix(0L, nrow = n_trials, ncol = n_doses)

  # Initialize vectors to track the current dose level for each trial
  current_dose_vec <- rep(1L, n_trials)

  # Initialize matrix to track which doses have been eliminated
  eliminated_mat <- matrix(FALSE, nrow = n_trials, ncol = n_doses)

  # Initialize flags to track active trials and cohorts completed
  active_trials <- rep(TRUE, n_trials)
  cohorts_completed <- rep(0L, n_trials)
  stop_reason <- rep(NA_character_, n_trials)

  # Initialize titration phase flags for accelerated dose escalation
  in_titration <- rep(titration, n_trials)
  ft_flag <- rep(titration, n_trials)

  # Expand cohort_size to match n_cohort length (handle various input formats)
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

  # Set maximum total patients allowed in a single trial
  max_total_pts <- n_cohort * cohort_size_vec[1]

  # Get dimensions of decision table for indexing
  max_col_decision <- ncol(decision_table)
  if (extrasafe) {
    max_col_stopping <- ncol(stopping_boundaries)
  }

  # ===== TITRATION PHASE (if enabled) =====
  # Accelerated dose escalation to quickly find approximate dose range
  if (titration) {
    # Generate random uniform matrix and compare to true toxicity probabilities
    z_mat <- matrix(runif(n_trials * n_doses), nrow = n_trials, ncol = n_doses) <
      matrix(rep(p_true, n_trials), nrow = n_trials, byrow = TRUE)

    # Find the first dose where DLT occurs for accelerated escalation
    first_dlt_dose <- apply(z_mat, 1, function(row) {
      idx <- which(row)[1]
      if (is.na(idx)) n_doses + 1 else idx
    })

    # Treat one patient at each dose up to (and including) the first DLT dose
    for (trial in 1:n_trials) {
      if (first_dlt_dose[trial] > n_doses) {
        # No DLT observed: treat one patient at each dose level
        n_pts_mat[trial, 1:n_doses] <- 1L
        current_dose_vec[trial] <- n_doses
      } else {
        # DLT observed: treat patients at doses 1 to d, with DLT at dose d
        d <- first_dlt_dose[trial]
        n_pts_mat[trial, 1:d] <- 1L
        n_tox_mat[trial, d] <- 1L
        current_dose_vec[trial] <- d
      }
    }

    # First titration cohort is complete
    cohorts_completed[1:n_trials] <- 1L

    # Fill up remaining patients in the first cohort at the current dose
    if (cohort_size > 1) {
      fillup_count <- cohort_size - 1L

      for (trial in 1:n_trials) {
        if (ft_flag[trial]) {
          d <- current_dose_vec[trial]

          # Generate DLT events for remaining patients in the cohort
          if (fillup_count == 1) {
            u <- runif(1)
            dlt_fillup <- as.integer(u < p_true[d])
          } else {
            u_fillup <- runif(fillup_count)
            dlt_fillup <- sum(u_fillup < p_true[d])
          }

          # Update patient and DLT counts
          n_pts_mat[trial, d] <- n_pts_mat[trial, d] + fillup_count
          n_tox_mat[trial, d] <- n_tox_mat[trial, d] + dlt_fillup

          # Exit titration phase
          ft_flag[trial] <- FALSE
          in_titration[trial] <- FALSE
        }
      }
    } else {
      in_titration[1:n_trials] <- FALSE
    }
  }

  # ===== EVALUATE DECISIONS AFTER TITRATION PHASE =====
  # Apply BOIN decision rules to determine escalation/de-escalation
  if (titration) {
    for (trial in 1:n_trials) {
      if (!active_trials[trial]) next

      d <- current_dose_vec[trial]
      n_pts_d <- n_pts_mat[trial, d]
      n_tox_d <- n_tox_mat[trial, d]

      if (!is.na(n_pts_d)) {
        # Constrain indices to decision table bounds
        n_val_eli <- pmin(n_pts_d, max_col_decision)
        y_val_eli <- pmin(n_tox_d, n_val_eli)

        if (!is.na(decision_table[y_val_eli + 1L, n_val_eli])) {
          elim_decision <- decision_table[y_val_eli + 1L, n_val_eli]

          # Handle dose elimination decision (DE)
          if (elim_decision == "DE") {
            eliminated_mat[trial, d:n_doses] <- TRUE
            if (d > 1) {
              valid_doses <- which(!eliminated_mat[trial, 1:(d - 1)])
              current_dose_vec[trial] <- if (length(valid_doses) > 0) max(valid_doses) else 1L
            } else {
              active_trials[trial] <- FALSE
              stop_reason[trial] <- "lowest_dose_eliminated"
            }
          } else if (elim_decision == "E" && d < n_doses && !eliminated_mat[trial, d + 1L]) {
            # Handle escalation decision (E)
            current_dose_vec[trial] <- d + 1L
          } else if (elim_decision == "D" && d > 1) {
            # Handle de-escalation decision (D)
            valid_doses <- which(!eliminated_mat[trial, 1:(d - 1)])
            current_dose_vec[trial] <- if (length(valid_doses) > 0) max(valid_doses) else 1L
          }
        }
      }

      # Apply extra safety stopping rule for dose 1 after titration if enabled
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

  # ===== MAIN TRIAL LOOP: Cohort-by-cohort simulation =====
  start_cohort <- if (titration) 2L else 1L
  for (cohort in start_cohort:n_cohort) {
    # Identify trials still actively running (excluding those in titration)
    normal_active_idx <- which(active_trials & !in_titration)
    n_normal <- length(normal_active_idx)

    if (n_normal == 0) break

    # Get the cohort size for the current cohort
    current_cohort_size <- cohort_size_vec[cohort]

    # ===== EARLY STOPPING RULE (simple) =====
    # Stop if sample size at current dose reaches n_earlystop threshold
    if (n_earlystop_rule == "simple") {
      early_stop_check <- n_pts_mat[cbind(normal_active_idx, current_dose_vec[normal_active_idx])] >= n_earlystop
      if (any(early_stop_check)) {
        early_stop_trials <- normal_active_idx[early_stop_check]
        active_trials[early_stop_trials] <- FALSE
        cohorts_completed[early_stop_trials] <- cohort
        stop_reason[early_stop_trials] <- "n_earlystop_simple"

        normal_active_idx <- which(active_trials & !in_titration)
        n_normal <- length(normal_active_idx)
        if (n_normal == 0) next
      }
    }

    # ===== CHECK MAXIMUM SAMPLE SIZE CONSTRAINT =====
    # Check if adding the next cohort would exceed maximum sample size
    current_total_pts <- rowSums(n_pts_mat[normal_active_idx, , drop = FALSE])
    will_exceed <- (current_total_pts + current_cohort_size) > max_total_pts

    if (any(will_exceed)) {
      exceed_idx <- which(will_exceed)
      exceed_trials <- normal_active_idx[exceed_idx]
      remaining_pts <- max_total_pts - current_total_pts[exceed_idx]

      # Treat remaining patients at current dose to reach maximum
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

      normal_active_idx <- which(active_trials & !in_titration)
      n_normal <- length(normal_active_idx)
      if (n_normal == 0) next
    }

    # ===== GENERATE DLT OUTCOMES FOR THE COHORT =====
    # Treat patients in the current cohort and generate DLT outcomes
    if (n_normal > 0) {
      current_doses <- current_dose_vec[normal_active_idx]
      p_true_current <- p_true[current_doses]

      if (current_cohort_size == 1) {
        # Single patient cohort
        u <- runif(n_normal)
        dlt_counts <- as.integer(u < p_true_current)
      } else {
        # Multi-patient cohort: generate DLT for each patient
        u_mat <- matrix(runif(n_normal * current_cohort_size),
                        nrow = n_normal, ncol = current_cohort_size)
        dlt_counts <- rowSums(u_mat < p_true_current)
      }

      # Update patient and DLT counts in matrices
      update_idx <- cbind(normal_active_idx, current_doses)
      n_pts_mat[update_idx] <- n_pts_mat[update_idx] + current_cohort_size
      n_tox_mat[update_idx] <- n_tox_mat[update_idx] + dlt_counts
    } else {
      next
    }

    # Record cohort number for each trial
    cohorts_completed[normal_active_idx] <- cohort

    # ===== EXTRA SAFETY STOPPING RULE (if enabled) =====
    # Apply safety stopping rule for dose 1 when extrasafe = TRUE
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

            normal_active_idx <- which(active_trials & !in_titration)
            n_normal <- length(normal_active_idx)
            if (n_normal == 0) next
          }
        }
      }
    }

    # ===== APPLY BOIN DECISION RULES =====
    # Determine escalation/de-escalation based on updated data
    current_doses <- current_dose_vec[normal_active_idx]
    n_pts_current <- n_pts_mat[cbind(normal_active_idx, current_doses)]
    n_tox_current <- n_tox_mat[cbind(normal_active_idx, current_doses)]

    # Constrain indices to decision table bounds
    n_vals <- pmin(n_pts_current, max_col_decision)
    y_vals <- pmin(n_tox_current, n_vals)

    # Look up decisions from decision table
    decisions <- decision_table[cbind(y_vals + 1L, n_vals)]

    # ===== HANDLE ESCALATION DECISIONS =====
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

      # Escalate to next dose for trials where escalation is possible
      if (any(can_escalate)) {
        current_dose_vec[esc_trials[can_escalate]] <- esc_doses[can_escalate] + 1L
      }
    }

    # ===== HANDLE DE-ESCALATION AND ELIMINATION DECISIONS =====
    deesc_idx <- which(decisions %in% c("D", "DE"))
    if (length(deesc_idx) > 0) {
      deesc_trials <- normal_active_idx[deesc_idx]
      deesc_doses <- current_doses[deesc_idx]

      # Handle dose elimination (DE decision)
      elim_idx <- which(decisions[deesc_idx] == "DE")
      if (length(elim_idx) > 0) {
        elim_trials <- deesc_trials[elim_idx]
        elim_doses <- deesc_doses[elim_idx]

        # Mark current dose and all higher doses as eliminated
        for (i in seq_along(elim_trials)) {
          trial <- elim_trials[i]
          dose <- elim_doses[i]
          eliminated_mat[trial, dose:n_doses] <- TRUE
        }

        # Stop trial if lowest dose is eliminated
        dose1_elim <- elim_doses == 1
        if (any(dose1_elim)) {
          stopped <- elim_trials[dose1_elim]
          active_trials[stopped] <- FALSE
          cohorts_completed[stopped] <- cohort
          stop_reason[stopped] <- "lowest_dose_eliminated"
        }

        # De-escalate to highest valid dose below eliminated dose
        need_deesc <- elim_doses > 1
        if (any(need_deesc)) {
          deesc_t <- elim_trials[need_deesc]
          deesc_d <- elim_doses[need_deesc]

          for (j in seq_along(deesc_t)) {
            trial <- deesc_t[j]
            dose <- deesc_d[j]
            valid_doses <- which(!eliminated_mat[trial, 1:(dose - 1)])
            current_dose_vec[trial] <- if (length(valid_doses) > 0) max(valid_doses) else 1L
          }
        }

        normal_active_idx <- which(active_trials & !in_titration)
        n_normal <- length(normal_active_idx)
        if (n_normal == 0) next
      }

      # Handle de-escalation without elimination (D decision)
      no_elim_idx <- which(decisions[deesc_idx] == "D")
      if (length(no_elim_idx) > 0) {
        no_elim_trials <- deesc_trials[no_elim_idx]
        no_elim_doses <- deesc_doses[no_elim_idx]

        can_deescalate <- no_elim_doses > 1
        if (any(can_deescalate)) {
          de_trials <- no_elim_trials[can_deescalate]
          de_doses <- no_elim_doses[can_deescalate]

          # De-escalate to highest valid dose below current dose
          for (i in seq_along(de_trials)) {
            trial <- de_trials[i]
            dose <- de_doses[i]
            valid_doses <- which(!eliminated_mat[trial, 1:(dose - 1)])
            current_dose_vec[trial] <- if (length(valid_doses) > 0) max(valid_doses) else 1L
          }
        }
      }
    }

    # ===== EARLY STOPPING RULE (with_stay) =====
    # Stop if sample size at current dose reaches n_earlystop threshold
    # and the decision indicates staying at current dose
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

        # Determine if trial should stop: dose stays (S decision or boundary constraints)
        is_stay <- (decisions_check == "S") |
          (dose_check == 1 & decisions_check %in% c("D", "DE")) |
          (dose_check == n_doses & decisions_check == "E")

        # Check if escalation is blocked by elimination
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

        # Stop trials where dose stays at current level
        if (any(is_stay)) {
          stop_trials <- normal_active_idx[check_idx[is_stay]]
          active_trials[stop_trials] <- FALSE
          cohorts_completed[stop_trials] <- cohort
          stop_reason[stop_trials] <- "n_earlystop_with_stay"

          normal_active_idx <- which(active_trials & !in_titration)
          n_normal <- length(normal_active_idx)
          if (n_normal == 0) next
        }
      }
    }
  }

  # ===== FINALIZE RESULTS =====
  # Mark any remaining active trials as stopped due to max cohorts reached
  still_active <- which(active_trials)
  if (length(still_active) > 0) {
    active_trials[still_active] <- FALSE
    stop_reason[still_active] <- "max_cohorts_reached"
  }

  # Return simulation results as a list
  return(list(
    n_pts_all = n_pts_mat,
    n_tox_all = n_tox_mat,
    eliminated_mat = eliminated_mat,
    cohorts_completed = cohorts_completed,
    stop_reason = stop_reason
  ))
}
