#' Run BOIN Simulation
#'
#' @description
#'   Execute multiple BOIN trial simulations using optimized batch processing for
#'   computational efficiency. The function automatically generates decision tables
#'   and stopping boundaries based on the target toxicity rate and design parameters,
#'   then processes all trials simultaneously at each cohort for maximum performance.
#'
#'   This implementation is particularly beneficial for large-scale simulations
#'   (e.g., 10,000+ trials) where computational efficiency is critical for practical
#'   protocol development and regulatory submissions.
#'
#' @param n_trials Numeric. Number of trials to simulate. Default is 10000.
#'   Typically 1000-10000 for detailed operating characteristics.
#' @param target Numeric. The target toxicity probability (e.g., 0.30 for 30%).
#' @param p_true Numeric vector. True toxicity probabilities for each dose.
#'   Length determines number of doses evaluated.
#' @param n_doses Numeric. Number of doses evaluated.
#' @param n_cohort Numeric. Maximum number of cohorts per trial.
#' @param cohort_size Numeric vector or scalar specifying patients per cohort.
#'   If vector (e.g., c(4, 3, 3)), each element specifies size for corresponding cohort.
#'   If scalar, all cohorts use the same size.
#' @param n_earlystop Numeric. Sample size at current dose triggering trial termination.
#'   Default is 18. This is also used as the maximum sample size for decision table
#'   generation.
#' @param cutoff_eli Numeric. Cutoff probability for dose elimination. Default is 0.95.
#'   If Pr(toxicity > target | data) > cutoff_eli, eliminate the dose and higher doses.
#' @param extrasafe Logical. If TRUE, apply additional safety stopping rule at the
#'   lowest dose. Default is FALSE. When TRUE, the trial stops early if the lowest
#'   dose is overly toxic based on \code{cutoff_eli - offset}.
#' @param offset Numeric. Offset for safety stopping boundary when \code{extrasafe = TRUE}.
#'   Default is 0.05. The safety cutoff becomes \code{cutoff_eli - offset}.
#' @param min_mtd_sample Numeric. Minimum sample size for MTD consideration. Default is 6.
#' @param boundMTD Logical. If TRUE, impose the condition that the isotonic estimate of
#'   toxicity probability for the selected MTD must be less than the de-escalation boundary.
#'   Default is FALSE. This provides a more conservative MTD selection.
#' @param n_earlystop_rule Character. Rule for early stopping at n_earlystop.
#'   Options are "simple" (stop when n >= n_earlystop) or "with_stay" (stop when
#'   n >= n_earlystop AND decision = "Stay"). Default is "simple" for backward compatibility.
#'   "with_stay" follows the BOIN standard implementation.
#' @param seed Numeric. Random seed for reproducibility. Default is 123.
#'
#' @return A list containing:
#'   \item{detailed_results}{List of results from each individual trial}
#'   \item{summary}{Aggregated summary statistics from \code{summarize_simulation_boin()}}
#'
#' @details
#'   This function implements an optimized approach to BOIN simulation:
#'   \enumerate{
#'     \item Generate BOIN boundaries using \code{get_boin_boundary(target)}
#'     \item Generate decision table using \code{get_boin_decision()}
#'     \item If \code{extrasafe = TRUE}, generate safety stopping boundaries
#'     \item Initialize matrices to store state for all trials simultaneously
#'     \item Loop over cohorts (not trials), processing all active trials at each cohort
#'     \item Generate DLT data for all active trials in a single \code{rbinom()} call
#'     \item Apply dose decisions using vectorized matrix operations
#'     \item Perform MTD selection using optimized batch processing
#'   }
#'
#'   The optimized approach reduces the number of iterations from
#'   \code{n_trials * n_cohort} to just \code{n_cohort}, while processing
#'   all trials simultaneously using matrix operations. This typically provides
#'   2-5x speedup compared to traditional sequential loop-based approaches.
#'
#'   **Key Performance Features:**
#'   \itemize{
#'     \item Automatic generation of decision tables and boundaries
#'     \item Single random number generation for all trials per cohort
#'     \item Matrix-based decision table lookups via vectorized indexing
#'     \item Batch isotonic regression and MTD selection
#'     \item Early termination of inactive trials to reduce overhead
#'   }
#'
#'   **Safety Stopping Rule:**
#'   When \code{extrasafe = TRUE}, the trial implements an additional safety rule:
#'   if the lowest dose shows excessive toxicity based on posterior probability
#'   Pr(toxicity > target | data) > \code{cutoff_eli - offset}, the entire trial
#'   stops for safety. This rule requires at least 3 patients at the lowest dose.
#'
#'   **boundMTD:**
#'   When \code{boundMTD = TRUE}, the selected MTD must satisfy the condition that
#'   its isotonic-estimated toxicity rate is below the de-escalation boundary.
#'   This provides a more conservative MTD selection, ensuring the selected dose
#'   is not too close to overly toxic doses.
#'
#'   **n_earlystop_rule:**
#'   \itemize{
#'     \item "simple": Stop when n >= n_earlystop (default, backward compatible)
#'     \item "with_stay": Stop when n >= n_earlystop AND next decision = "Stay".
#'       This follows the BOIN standard implementation and ensures the algorithm
#'       has converged before stopping.
#'   }
#'
#'   **Memory Usage:** The function maintains \code{n_trials x n_doses} matrices
#'   for patient counts, DLT counts, and elimination status. For typical
#'   simulations (10000 trials, 5-9 doses), memory usage is negligible on
#'   modern systems (<100MB).
#'
#' @references
#'   Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'   Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @examples
#' \dontrun{
#' # Basic BOIN simulation (backward compatible)
#' result <- sim_boin(
#'   n_trials = 10000,
#'   target = 0.30,
#'   p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
#'   n_doses = 5,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   seed = 123
#' )
#'
#' # With BOIN standard implementation (boundMTD + with_stay)
#' result_standard <- sim_boin(
#'   n_trials = 10000,
#'   target = 0.30,
#'   p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
#'   n_doses = 5,
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   boundMTD = TRUE,
#'   n_earlystop_rule = "with_stay",
#'   seed = 123
#' )
#'
#' # Conservative design with extra safety
#' result_conservative <- sim_boin(
#'   n_trials = 10000,
#'   target = 0.30,
#'   p_true = seq(0.05, 0.45, by = 0.05),
#'   n_doses = 9,
#'   n_cohort = 48,
#'   cohort_size = 3,
#'   extrasafe = TRUE,
#'   boundMTD = TRUE,
#'   n_earlystop_rule = "with_stay",
#'   seed = 123
#' )
#' }
#'
#' @seealso
#'   \code{\link{get_boin_boundary}} for BOIN boundary calculation
#'   \code{\link{get_boin_decision}} for decision table generation
#'   \code{\link{get_boin_stopping_boundaries}} for safety stopping boundaries
#'   \code{\link{summarize_simulation_boin}} for summarizing simulation results
#'
#' @export
sim_boin <- function(
    n_trials = 10000,
    target,
    p_true,
    n_doses,
    n_cohort,
    cohort_size,
    n_earlystop = 18,
    cutoff_eli = 0.95,
    extrasafe = FALSE,
    offset = 0.05,
    min_mtd_sample = 6,
    boundMTD = FALSE,
    n_earlystop_rule = c("simple", "with_stay"),
    seed = 123
) {

  # Validate n_earlystop_rule argument
  n_earlystop_rule <- match.arg(n_earlystop_rule)

  set.seed(seed)

  cat("========================================\n")
  cat("Starting BOIN Simulation\n")
  cat("Number of trials:", n_trials, "\n")
  cat("Target DLT rate:", target * 100, "%\n")
  cat("Number of doses:", n_doses, "\n")
  if (extrasafe) {
    cat("Extra safety: Enabled (offset =", offset, ")\n")
  }
  if (boundMTD) {
    cat("boundMTD: Enabled (conservative MTD selection)\n")
  }
  cat("Early stop rule:", n_earlystop_rule, "\n")
  cat("========================================\n\n")

  # ========== Generate Decision Tables ==========
  cat("Generating BOIN decision tables...\n")

  # Get BOIN boundaries
  boin_bound <- get_boin_boundary(target)

  # Generate decision table
  decision_table <- get_boin_decision(
    target = target,
    lambda_e = boin_bound$lambda_e,
    lambda_d = boin_bound$lambda_d,
    max_sample_size = n_earlystop,
    cutoff_eli = cutoff_eli
  )

  # Generate safety stopping boundaries if extrasafe is enabled
  if (extrasafe) {
    cutoff_stop <- cutoff_eli - offset
    stopping_boundaries <- get_boin_stopping_boundaries(
      target = target,
      max_sample_size = n_earlystop,
      cutoff_stop = cutoff_stop
    )
    cat("Safety stopping boundary generated (cutoff =", cutoff_stop, ")\n")
  } else {
    stopping_boundaries <- NULL
  }

  cat("Decision tables ready.\n\n")

  # ========== Initialization ==========
  # Initialize matrices for all trials (rows = trials, cols = doses)
  n_pts_mat <- matrix(0L, nrow = n_trials, ncol = n_doses)
  n_tox_mat <- matrix(0L, nrow = n_trials, ncol = n_doses)
  current_dose_vec <- rep(1L, n_trials)
  eliminated_mat <- matrix(FALSE, nrow = n_trials, ncol = n_doses)
  active_trials <- rep(TRUE, n_trials)
  cohorts_completed <- rep(0L, n_trials)
  stop_reason <- rep(NA_character_, n_trials)  # Track stopping reason for each trial

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

  # Pre-compute table dimensions
  max_col_decision <- ncol(decision_table)
  if (extrasafe) {
    max_col_stopping <- ncol(stopping_boundaries)
  }

  # ========== Main Cohort Loop ==========
  for (cohort in 1:n_cohort) {

    # Get active trial indices
    active_idx <- which(active_trials)
    n_active <- length(active_idx)

    if (n_active == 0) break

    # Progress reporting
    if (cohort %% max(1, n_cohort %/% 5) == 0 || cohort == 1) {
      cat("Cohort:", cohort, "/", n_cohort,
          "| Active trials:", n_active, "/", n_trials, "\n")
    }

    # Get current cohort size
    current_cohort_size <- cohort_size_vec[cohort]

    # ========== Early Stopping Check ==========
    if (n_earlystop_rule == "simple") {
      # Simple rule: stop if n >= n_earlystop
      early_stop_check <- n_pts_mat[cbind(active_idx, current_dose_vec[active_idx])] >= n_earlystop
      if (any(early_stop_check)) {
        early_stop_trials <- active_idx[early_stop_check]
        active_trials[early_stop_trials] <- FALSE
        cohorts_completed[early_stop_trials] <- cohort - 1L
        stop_reason[early_stop_trials] <- "n_earlystop_reached"

        # Update active indices
        active_idx <- active_idx[!early_stop_check]
        n_active <- length(active_idx)
        if (n_active == 0) break
      }
    } else if (n_earlystop_rule == "with_stay") {
      # BOIN standard: stop if n >= n_earlystop AND decision = "Stay"
      current_doses <- current_dose_vec[active_idx]
      current_n_pts <- n_pts_mat[cbind(active_idx, current_doses)]
      current_n_tox <- n_tox_mat[cbind(active_idx, current_doses)]

      # Check which trials have n >= n_earlystop
      n_check <- current_n_pts >= n_earlystop

      if (any(n_check)) {
        # Get decisions for trials meeting n_earlystop
        check_idx <- which(n_check)
        n_vals <- pmin(current_n_pts[check_idx], max_col_decision)
        y_vals <- pmin(current_n_tox[check_idx], n_vals)

        decisions <- decision_table[cbind(y_vals + 1L, n_vals)]

        # Determine if decision is effectively "Stay"
        # Stay conditions:
        # 1. decision == "S"
        # 2. dose == 1 AND decision %in% c("D", "DE") (can't de-escalate)
        # 3. dose == n_doses AND decision == "E" (can't escalate)
        # 4. dose < n_doses AND next dose eliminated AND decision == "E"

        dose_check <- current_doses[check_idx]
        is_stay <- (decisions == "S") |
          (dose_check == 1 & decisions %in% c("D", "DE")) |
          (dose_check == n_doses & decisions == "E")

        # Check for eliminated next dose
        if (any(!is_stay & dose_check < n_doses & decisions == "E")) {
          elim_check_idx <- which(!is_stay & dose_check < n_doses & decisions == "E")
          for (i in elim_check_idx) {
            trial <- active_idx[check_idx[i]]
            next_dose <- dose_check[i] + 1
            if (eliminated_mat[trial, next_dose]) {
              is_stay[i] <- TRUE
            }
          }
        }

        # Stop trials with Stay decision
        if (any(is_stay)) {
          stop_trials <- active_idx[check_idx[is_stay]]
          active_trials[stop_trials] <- FALSE
          cohorts_completed[stop_trials] <- cohort - 1L
          stop_reason[stop_trials] <- "n_earlystop_with_stay"

          # Update active indices
          active_idx <- setdiff(active_idx, stop_trials)
          n_active <- length(active_idx)
          if (n_active == 0) break
        }
      }
    }

    # ========== Generate DLT Data (Vectorized) ==========
    current_doses <- current_dose_vec[active_idx]
    p_true_current <- p_true[current_doses]
    dlt_counts <- rbinom(n_active, current_cohort_size, p_true_current)

    # Update n_pts and n_tox matrices
    update_idx <- cbind(active_idx, current_doses)
    n_pts_mat[update_idx] <- n_pts_mat[update_idx] + current_cohort_size
    n_tox_mat[update_idx] <- n_tox_mat[update_idx] + dlt_counts

    # ========== Safety Stopping Rule at Lowest Dose ==========
    if (extrasafe) {
      dose1_trials <- active_idx[current_doses == 1]
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

            # Update active indices
            active_idx <- setdiff(active_idx, stopped_trials)
            n_active <- length(active_idx)
            if (n_active == 0) break
          }
        }
      }
    }

    # ========== Dose Decision (Vectorized) ==========
    current_doses <- current_dose_vec[active_idx]
    n_pts_current <- n_pts_mat[cbind(active_idx, current_doses)]
    n_tox_current <- n_tox_mat[cbind(active_idx, current_doses)]

    # Prepare indices for decision table lookup
    n_vals <- pmin(n_pts_current, max_col_decision)
    y_vals <- pmin(n_tox_current, n_vals)

    # Lookup decisions
    decisions <- decision_table[cbind(y_vals + 1L, n_vals)]

    # Process Escalate decisions
    esc_idx <- which(decisions == "E")
    if (length(esc_idx) > 0) {
      esc_trials <- active_idx[esc_idx]
      esc_doses <- current_doses[esc_idx]

      # First check if not at max dose
      not_at_max <- esc_doses < n_doses

      # Initialize can_escalate as FALSE
      can_escalate <- rep(FALSE, length(esc_doses))

      # Then check if next dose is not eliminated (only for those not at max)
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

    # Process De-escalate decisions
    deesc_idx <- which(decisions %in% c("D", "DE"))
    if (length(deesc_idx) > 0) {
      deesc_trials <- active_idx[deesc_idx]
      deesc_doses <- current_doses[deesc_idx]

      # Check for elimination
      elim_idx <- which(decisions[deesc_idx] == "DE")
      if (length(elim_idx) > 0) {
        elim_trials <- deesc_trials[elim_idx]
        elim_doses <- deesc_doses[elim_idx]

        # Process each elimination individually
        for (i in seq_along(elim_idx)) {
          trial_idx <- elim_trials[i]
          dose <- elim_doses[i]

          # Eliminate current dose and all higher doses
          eliminated_mat[trial_idx, dose:n_doses] <- TRUE

          if (dose > 1) {
            # De-escalate
            new_dose <- dose - 1L
            while (new_dose > 1 && eliminated_mat[trial_idx, new_dose]) {
              new_dose <- new_dose - 1L
            }
            current_dose_vec[trial_idx] <- new_dose
          } else {
            # Lowest dose eliminated - stop trial
            active_trials[trial_idx] <- FALSE
            cohorts_completed[trial_idx] <- cohort
            stop_reason[trial_idx] <- "lowest_dose_eliminated"
          }
        }

        # Update active indices
        active_idx <- which(active_trials)
        n_active <- length(active_idx)
        if (n_active == 0) break
      }

      # De-escalate without elimination
      no_elim_idx <- which(decisions[deesc_idx] == "D")
      if (length(no_elim_idx) > 0) {
        no_elim_trials <- deesc_trials[no_elim_idx]
        no_elim_doses <- deesc_doses[no_elim_idx]

        can_deescalate <- no_elim_doses > 1
        if (any(can_deescalate)) {
          de_trials <- no_elim_trials[can_deescalate]
          de_doses <- no_elim_doses[can_deescalate]

          # De-escalate and skip eliminated doses
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

    # --- Stay (S) - no action needed ---

    # Update cohorts completed for still-active trials
    cohorts_completed[active_idx] <- cohort
  }

  cat("\nCohort simulation completed!\n")
  cat("Performing MTD selection for all trials...\n\n")

  # ========== MTD Selection Phase ==========

  # Apply isotonic regression to all trials at once
  iso_est_mat <- .isotonic_regression_batch(
    n_pts_mat, n_tox_mat, eliminated_mat, min_mtd_sample
  )

  # Compute distances from target for all trials
  diffs_mat <- abs(iso_est_mat - target)

  # Select MTD for each trial (with boundMTD support)
  mtd_results <- .select_mtd_batch(
    diffs_mat, iso_est_mat, eliminated_mat,
    cohorts_completed, stop_reason, target,
    boundMTD = boundMTD,
    lambda_d = boin_bound$lambda_d
  )

  # Construct results list
  simulation_results <- vector("list", n_trials)
  for (trial in 1:n_trials) {
    simulation_results[[trial]] <- list(
      n_pts = n_pts_mat[trial, ],
      n_tox = n_tox_mat[trial, ],
      mtd = mtd_results$mtd[trial],
      iso_est = iso_est_mat[trial, ],
      reason = mtd_results$reason[trial],
      cohorts_completed = cohorts_completed[trial]
    )
  }

  cat("MTD selection completed!\n\n")

  # Compute summary statistics
  summary_result <- summarize_simulation_boin(simulation_results, n_doses, p_true)

  return(list(
    detailed_results = simulation_results,
    summary = summary_result
  ))
}
