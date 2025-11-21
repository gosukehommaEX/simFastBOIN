#' Generate TITE-BOIN Decision Table
#'
#' @description
#' Generate a comprehensive decision table for the Time-to-Event Bayesian Optimal
#' Interval (TITE-BOIN) design. The table contains dose escalation/de-escalation/stay
#' decisions based on the number of patients treated, observed toxicities, and
#' pending patients (those whose toxicity assessment is incomplete).
#'
#' The decision rules incorporate three scenarios:
#' \itemize{
#'   \item Rule 1: Accrual suspension when pending ratio is high and toxicity is low
#'   \item Rule 2: Standard BOIN decisions when all patients have completed follow-up
#'   \item Rule 3: Modified decisions accounting for pending patients using
#'         Standardized Total Follow-up Time (STFT)
#' }
#'
#' @param lambda_e Lower boundary for dose escalation (from \code{get_boin_boundary()})
#' @param lambda_d Upper boundary for dose de-escalation (from \code{get_boin_boundary()})
#' @param target Target toxicity rate (e.g., 0.30 for 30\%)
#' @param n_max Maximum number of patients in the trial
#' @param cohort_size Number of patients per cohort (default: 3)
#' @param max_accrual_pending_ratio Maximum allowable ratio of pending patients
#'   before suspending accrual (default: 0.51, meaning suspend when >49\% pending)
#' @param v Missing fraction threshold for decision rules (default: 0.25)
#' @param cutoff_eli Posterior probability cutoff for dose elimination (default: 0.95)
#' @param neli Minimum sample size required for dose elimination (default: 3)
#' @param prior_strength Prior strength parameter for toxicity rate estimation (default: 0.5)
#'
#' @return A data frame containing the TITE-BOIN decision table with columns:
#' \describe{
#'   \item{n_patients}{Number of patients treated at the dose}
#'   \item{n_tox}{Number of observed toxicities (completed assessments only)}
#'   \item{n_pending}{Number of patients with pending toxicity assessments}
#'   \item{suspend}{Accrual suspension condition}
#'   \item{escalate}{Dose escalation condition}
#'   \item{stay}{Stay at current dose condition}
#'   \item{de_escalate}{Dose de-escalation condition}
#' }
#'
#' @details
#' The TITE-BOIN design extends the standard BOIN design to handle late-onset
#' toxicities by making real-time decisions while some patients' toxicity outcomes
#' are still pending. The decision table uses the Standardized Total Follow-up Time
#' (STFT) to quantify the amount of follow-up information available from pending
#' patients.
#'
#' @references
#' Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K.E. (2018).
#' Time-to-Event Bayesian Optimal Interval Design to Accelerate Phase I Trials.
#' Clinical Cancer Research, 24(20): 4921-4930.
#'
#' Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I
#' Clinical Trials. Journal of the Royal Statistical Society: Series C, 64, 507-523.
#'
#' @importFrom dplyr %>% rename filter arrange mutate select case_when if_else
#' @importFrom stats pbeta
#'
#' @examples
#' # Generate TITE-BOIN decision table for a trial with target DLT rate of 30%
#' target <- 0.30
#' boundary <- get_boin_boundary(target)
#'
#' decision_table <- get_titeboin_decision(
#'   lambda_e = boundary$lambda_e,
#'   lambda_d = boundary$lambda_d,
#'   target = target,
#'   n_max = 18,
#'   cohort_size = 3,
#'   max_accrual_pending_ratio = 0.51,
#'   cutoff_eli = 0.95
#' )
#'
#' # View the first few rows
#' head(decision_table)
#'
#' @export
get_titeboin_decision <- function(
    lambda_e, lambda_d, target, n_max, cohort_size = 3,
    max_accrual_pending_ratio = 0.51, v = 0.25,
    cutoff_eli = 0.95, neli = 3, prior_strength = 0.5
) {
  
  n_seq <- seq(cohort_size, n_max, by = cohort_size)
  
  # Create all combinations of (n_patients, n_tox, n_pending)
  dt <- expand.grid(
    n = n_seq,
    y = 0:n_max,
    c = 0:n_max
  ) %>%
    rename(n_patients = n, n_tox = y, n_pending = c) %>%
    filter(n_tox <= n_patients, n_pending <= n_patients - n_tox) %>%
    arrange(n_patients, n_tox, n_pending)
  
  # Calculate basic statistics
  dt <- dt %>%
    mutate(
      obs_rate = n_tox / n_patients,
      mf = n_pending / n_patients,
      n_completed = n_patients - n_pending,
      p_hat = (n_tox + prior_strength * target) / (n_completed + 1)
    )
  
  # Determine which rule applies
  # RULE 1: Accrual suspension if pending ratio is high AND obs_rate is low (safe)
  # RULE 2: No pending patients (all completed)
  # RULE 3: Other cases with pending patients
  dt <- dt %>%
    mutate(
      rule1_flag = (mf > (1 - max_accrual_pending_ratio)) & (obs_rate < lambda_d),
      rule2_flag = (n_pending == 0),
      rule3_flag = !rule1_flag & !rule2_flag,
      rule3_esc = rule3_flag & (obs_rate < target),
      rule3_deesc = rule3_flag & (obs_rate > target),
      rule3_stay = rule3_flag & (obs_rate == target)
    )
  
  # Initialize decision columns
  dt <- dt %>%
    mutate(
      suspend = "No",
      escalate = "No",
      stay = "No",
      de_escalate = "No"
    )
  
  # ========== RULE 1: Accrual suspension ==========
  dt <- dt %>%
    mutate(
      suspend = if_else(rule1_flag, "Yes", suspend),
      escalate = if_else(rule1_flag, "No", escalate),
      stay = if_else(rule1_flag, "No", stay),
      de_escalate = if_else(rule1_flag, "No", de_escalate)
    )
  
  # ========== RULE 2: No pending patients ==========
  dt <- dt %>%
    mutate(
      prob_exceed_rule2 = case_when(
        rule2_flag ~ 1 - pbeta(target, n_tox + 1, n_patients - n_tox + 1),
        TRUE ~ NA_real_),
      eliminate_flag_rule2 = case_when(
        rule2_flag & (n_patients >= neli) & (prob_exceed_rule2 > cutoff_eli) ~ TRUE,
        TRUE ~ FALSE)
    )
  
  dt <- dt %>%
    mutate(
      escalate = case_when(
        rule2_flag & (obs_rate < lambda_e) ~ "Yes",
        TRUE ~ escalate),
      de_escalate = case_when(
        rule2_flag & (obs_rate > lambda_d) & eliminate_flag_rule2 ~ "Yes & Eliminate",
        rule2_flag & (obs_rate > lambda_d) & !eliminate_flag_rule2 ~ "Yes",
        TRUE ~ de_escalate),
      stay = case_when(
        rule2_flag & (obs_rate >= lambda_e) & (obs_rate <= lambda_d) ~ "Yes",
        TRUE ~ stay)
    )
  
  # ========== RULE 3a: Escalation with pending patients ==========
  dt <- dt %>%
    mutate(
      stft_e = case_when(
        rule3_esc ~ n_pending - (1 - p_hat) / p_hat * (n_patients * lambda_e - n_tox),
        TRUE ~ NA_real_)
    )
  
  dt <- dt %>%
    mutate(
      suspend = case_when(
        rule3_esc & (stft_e >= 0) ~ paste0("MF<", v, " & STFT>=", round(stft_e, 2)),
        rule3_esc & (stft_e < 0) ~ paste0("MF<", v),
        TRUE ~ suspend),
      escalate = case_when(
        rule3_esc & (stft_e >= 0) ~ paste0("MF>=", v, " & STFT>=", round(stft_e, 2)),
        rule3_esc & (stft_e < 0) ~ paste0("MF>=", v),
        TRUE ~ escalate),
      stay = case_when(
        rule3_esc & (stft_e >= 0) ~ paste0("STFT<", round(stft_e, 2)),
        rule3_esc & (stft_e < 0) ~ "No",
        TRUE ~ stay),
      de_escalate = case_when(
        rule3_esc ~ "No",
        TRUE ~ de_escalate)
    )
  
  # ========== RULE 3b: De-escalation with pending patients ==========
  dt <- dt %>%
    mutate(
      stft_d = case_when(
        rule3_deesc ~ n_pending - (1 - p_hat) / p_hat * (n_patients * lambda_d - n_tox),
        TRUE ~ NA_real_),
      prob_exceed_rule3 = case_when(
        rule3_deesc ~ 1 - pbeta(target, n_tox + 1, n_patients - n_tox + 1),
        TRUE ~ NA_real_),
      eliminate_flag_rule3 = case_when(
        rule3_deesc & (n_patients >= neli) & (prob_exceed_rule3 > cutoff_eli) ~ TRUE,
        TRUE ~ FALSE)
    )
  
  dt <- dt %>%
    mutate(
      suspend = case_when(
        rule3_deesc ~ "No",
        TRUE ~ suspend),
      escalate = case_when(
        rule3_deesc ~ "No",
        TRUE ~ escalate),
      stay = case_when(
        rule3_deesc & eliminate_flag_rule3 ~ "No",
        rule3_deesc & !eliminate_flag_rule3 & (stft_d >= n_pending) ~ "No",
        rule3_deesc & !eliminate_flag_rule3 & (stft_d >= 0) & (stft_d < n_pending) ~ paste0("STFT>", round(stft_d, 2)),
        rule3_deesc & !eliminate_flag_rule3 & (stft_d < 0) ~ "No",
        TRUE ~ stay),
      de_escalate = case_when(
        rule3_deesc & eliminate_flag_rule3 ~ "Yes & Eliminate",
        rule3_deesc & !eliminate_flag_rule3 & (stft_d >= n_pending) ~ "Yes",
        rule3_deesc & !eliminate_flag_rule3 & (stft_d >= 0) & (stft_d < n_pending) ~ paste0("STFT<=", round(stft_d, 2)),
        rule3_deesc & !eliminate_flag_rule3 & (stft_d < 0) ~ "Yes",
        TRUE ~ de_escalate)
    )
  
  # ========== RULE 3c: Stay with pending patients ==========
  dt <- dt %>%
    mutate(
      suspend = case_when(
        rule3_stay ~ "No",
        TRUE ~ suspend),
      escalate = case_when(
        rule3_stay ~ "No",
        TRUE ~ escalate),
      stay = case_when(
        rule3_stay ~ "Yes",
        TRUE ~ stay),
      de_escalate = case_when(
        rule3_stay ~ "No",
        TRUE ~ de_escalate)
    )
  
  # Clean up and return
  dt <- dt %>%
    select(n_patients, n_tox, n_pending, suspend, escalate, stay, de_escalate) %>%
    arrange(n_patients, n_tox, n_pending)
  
  return(dt)
}