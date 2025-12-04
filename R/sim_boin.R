#' Run BOIN Simulation with Operating Characteristics
#'
#' @description
#'   Execute multiple BOIN trial simulations and compute operating characteristics.
#'   Combines patient enrollment, toxicity tracking, isotonic regression, and MTD
#'   selection into a single streamlined function.
#'
#' @param n_trials Numeric. Number of trials to simulate. Default is 10000.
#' @param target Numeric. Target toxicity probability (e.g., 0.30 for 30%).
#' @param p_true Numeric vector. True toxicity probabilities for each dose.
#' @param n_cohort Numeric. Maximum number of cohorts per trial.
#' @param cohort_size Numeric vector or scalar. Patients per cohort.
#' @param n_earlystop Numeric. Sample size triggering early stopping. Default is 18.
#' @param cutoff_eli Numeric. Cutoff for dose elimination. Default is 0.95.
#' @param extrasafe Logical. Apply extra safety stopping rule. Default is FALSE.
#' @param offset Numeric. Offset for safety cutoff when extrasafe = TRUE. Default is 0.05.
#' @param n_earlystop_rule Character. Early stopping rule: "with_stay" or "simple". Default is "with_stay".
#' @param titration Logical. Perform accelerated dose titration phase. Default is FALSE.
#' @param min_mtd_sample Numeric. Minimum patients required for MTD consideration. Default is 1.
#' @param boundMTD Logical. Impose boundMTD constraint on MTD selection. Default is FALSE.
#' @param seed Numeric. Random seed for reproducibility. Default is 123.
#'
#' @return List containing detailed_results and summary
#'
#' @examples
#' \dontrun{
#' # Basic BOIN simulation
#' result <- sim_boin(
#'   n_trials = 1000,
#'   target = 0.30,
#'   p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   seed = 123
#' )
#' print(result$summary)
#'
#' # BOIN with boundMTD and with_stay
#' result_standard <- sim_boin(
#'   n_trials = 1000,
#'   target = 0.30,
#'   p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   boundMTD = TRUE,
#'   n_earlystop_rule = "with_stay",
#'   seed = 123
#' )
#' print(result_standard$summary)
#'
#' # BOIN with extra safety
#' result_safe <- sim_boin(
#'   n_trials = 1000,
#'   target = 0.30,
#'   p_true = c(0.10, 0.25, 0.40, 0.55, 0.70),
#'   n_cohort = 10,
#'   cohort_size = 3,
#'   extrasafe = TRUE,
#'   offset = 0.05,
#'   seed = 123
#' )
#' print(result_safe$summary)
#' }
#'
#' @importFrom stats runif
#'
#' @export
sim_boin <- function(n_trials = 10000,
                     target,
                     p_true,
                     n_cohort,
                     cohort_size,
                     n_earlystop = 18,
                     cutoff_eli = 0.95,
                     extrasafe = FALSE,
                     offset = 0.05,
                     n_earlystop_rule = "with_stay",
                     titration = FALSE,
                     min_mtd_sample = 1,
                     boundMTD = FALSE,
                     seed = 123) {

  n_earlystop_rule <- match.arg(n_earlystop_rule, c("with_stay", "simple"))
  n_doses <- length(p_true)

  cat("========================================\n")
  cat("BOIN Simulation\n")
  cat("Trials:", n_trials, "| Target:", target * 100, "%",
      "| Doses:", n_doses, "\n")
  cat("========================================\n\n")

  # Step 1: Generate patient and toxicity data
  cat("Step 1: Generating patient enrollment and toxicity data...\n")

  pts_tox_result <- get_pts_and_tox(
    n_trials = n_trials,
    target = target,
    p_true = p_true,
    n_cohort = n_cohort,
    cohort_size = cohort_size,
    n_earlystop = n_earlystop,
    cutoff_eli = cutoff_eli,
    extrasafe = extrasafe,
    offset = offset,
    n_earlystop_rule = n_earlystop_rule,
    titration = titration,
    seed = seed
  )

  n_pts_mat <- pts_tox_result$n_pts_all
  n_tox_mat <- pts_tox_result$n_tox_all
  eliminated_mat <- pts_tox_result$eliminated_mat
  stop_reason_from_trial <- pts_tox_result$stop_reason
  cohorts_completed <- pts_tox_result$cohorts_completed

  cat("  Completed.\n\n")

  # Step 2: Apply isotonic regression
  cat("Step 2: Applying isotonic regression...\n")

  iso_est_mat <- isotonic_regression(n_pts_mat, n_tox_mat)

  cat("  Completed.\n\n")

  # Step 3: Select MTD for each trial
  cat("Step 3: Selecting MTD for each trial...\n")

  boin_bound <- get_boin_boundary(target)
  lambda_d <- boin_bound$lambda_d

  mtd_results <- select_mtd(
    iso_est_mat = iso_est_mat,
    n_pts_mat = n_pts_mat,
    eliminated_mat = eliminated_mat,
    target = target,
    boundMTD = boundMTD,
    lambda_d = if (boundMTD) lambda_d else NULL,
    min_mtd_sample = min_mtd_sample
  )

  # Use stop_reason from get_pts_and_tox when trial stopped for safety
  final_reason <- mtd_results$reason
  safety_stop_reasons <- c("lowest_dose_eliminated", "lowest_dose_too_toxic")
  override_idx <- stop_reason_from_trial %in% safety_stop_reasons

  if (any(override_idx)) {
    final_reason[override_idx] <- stop_reason_from_trial[override_idx]
    mtd_results$mtd[override_idx] <- NA_integer_
  }

  mtd_vector <- mtd_results$mtd

  cat("  Completed.\n\n")

  # Step 4: Compute operating characteristics
  cat("Step 4: Computing operating characteristics...\n")

  # Construct detailed results
  detailed_results <- lapply(seq_len(n_trials), function(i) {
    list(
      n_pts = n_pts_mat[i, ],
      n_tox = n_tox_mat[i, ],
      mtd = mtd_vector[i],
      iso_est = iso_est_mat[i, ],
      reason = final_reason[i],
      cohorts_completed = cohorts_completed[i]
    )
  })

  # MTD selection percentage
  mtd_selected <- matrix(0, nrow = n_trials, ncol = n_doses)
  valid_mtd <- !is.na(mtd_vector)
  if (any(valid_mtd)) {
    mtd_selected[cbind(which(valid_mtd), mtd_vector[valid_mtd])] <- 1
  }

  mtd_selection_by_dose <- colMeans(mtd_selected) * 100
  percent_no_mtd <- (1 - mean(valid_mtd)) * 100
  mtd_selection_percent <- c(mtd_selection_by_dose, percent_no_mtd)

  # Summary statistics
  summary_obj <- list(
    p_true = p_true,
    mtd_selection_percent = mtd_selection_percent,
    avg_n_pts = colMeans(n_pts_mat),
    avg_n_tox = colMeans(n_tox_mat),
    avg_total_n_pts = mean(rowSums(n_pts_mat)),
    avg_total_n_tox = mean(rowSums(n_tox_mat))
  )

  class(summary_obj) <- c("boin_summary", "list")

  cat("  Completed.\n\n")

  cat("========================================\n")
  cat("Simulation Summary\n")
  cat("========================================\n\n")

  print(summary_obj)

  cat("\n")

  return(list(
    detailed_results = detailed_results,
    summary = summary_obj
  ))
}
