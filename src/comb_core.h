#ifndef SIMFASTBOIN_COMB_CORE_H
#define SIMFASTBOIN_COMB_CORE_H

// Core algorithms of the BOIN drug-combination design of Lin and Yin (2017),
// written in plain C++ so that the same code can be compiled both inside the
// package (against R's numerical routines) and inside a standalone
// verification harness.
//
// The caller supplies post_prob_above() and beta_cdf().
//
// A dose matrix with J levels of drug A (rows) and K levels of drug B
// (columns) is stored column-major, as R stores a matrix: the combination in
// 0-based row a and column b is element a + J * b.
//
// Integer boundary vectors are indexed by the number of patients minus one.
// A value of 0 in b_elim means that no elimination boundary exists at that
// sample size, which corresponds to NA in the reference implementation.

#include <vector>
#include <algorithm>
#include <cmath>

namespace simfastboin {

// Pr(p > target | y DLTs out of n patients) under a uniform Beta(1, 1) prior.
// Defined by the caller.
double post_prob_above(double target, int y, int n);

// Distribution function of the beta distribution. Defined by the caller.
double beta_cdf(double x, double shape1, double shape2);

// Stop codes returned by comb_simulate_one().
enum CombStopCode {
  COMB_STOP_LOWEST_ELIMINATED = 0,
  COMB_STOP_LOWEST_TOO_TOXIC = 1,
  COMB_STOP_N_EARLYSTOP = 2,
  COMB_STOP_MAX_COHORTS = 3
};

// Reason codes returned by comb_choose_mtd().
enum CombMtdReason {
  COMB_MTD_SELECTED = 0,
  COMB_MTD_LOWEST_ELIMINATED = 1,
  COMB_MTD_NO_DOSE_BELOW_LAMBDA_D = 2
};

// Score of a neighboring combination considered for the next cohort: the
// posterior probability that its DLT probability lies in (lambda_e, lambda_d)
// under a Beta(0.5, 0.5) prior, plus 0.0005 per patient already treated there,
// which is how the reference implementation favors the better informed of
// two otherwise equal candidates.
inline double comb_candidate_score(int y, int n, double lambda_e,
                                   double lambda_d) {
  const double shape1 = y + 0.5;
  const double shape2 = (n - y) + 0.5;
  const double prob = beta_cdf(lambda_d, shape1, shape2) -
    beta_cdf(lambda_e, shape1, shape2);
  const double bonus = n * 5e-4;
  return prob + bonus;
}

// Choose among up to two candidate moves. Candidates that are not available
// have a score of zero. When no candidate is available the dose stays and no
// random variate is drawn; otherwise one variate is drawn to break ties, even
// when there is no tie, as in the reference implementation.
template <class Rng>
inline int comb_pick(const double* score, Rng& rng) {
  const double best = std::max(score[0], score[1]);
  if (best == 0.0) return -1;
  int ties[2];
  int n_ties = 0;
  for (int k = 0; k < 2; ++k) {
    if (score[k] == best) ties[n_ties++] = k;
  }
  const double u = rng();
  const int pick = static_cast<int>(u * n_ties + 1.0) - 1;
  return ties[pick];
}

// Actions returned by comb_decide().
enum CombAction {
  COMB_ACTION_STAY = 0,
  COMB_ACTION_ESCALATE = 1,
  COMB_ACTION_DEESCALATE = 2,
  COMB_ACTION_STOP = 3
};

// The decision after a cohort at the current combination (row a, column b,
// 0-based), once elimination and the safety rule have been applied: stop for
// n_earlystop, escalate, de-escalate or stay. For an escalation or a
// de-escalation, score[0] and score[1] receive the scores of the move by one
// level of drug A and of drug B, zero for a move that is not available; when
// both are zero the dose stays. Shared by the simulation and by the
// trial-conduct function, so that the two apply the same rule.
inline int comb_decide(const int* n, const int* y, const char* elim,
                       int n_rows, int n_cols, int a, int b, int n_earlystop,
                       double lambda_e, double lambda_d,
                       const std::vector<int>& b_esc,
                       const std::vector<int>& b_deesc, double* score) {
  const int J = n_rows;
  const int K = n_cols;
  const int max_n = static_cast<int>(b_esc.size());
  const int idx = a + J * b;
  int nc = n[idx];
  if (nc > max_n) nc = max_n;
  const int ix = nc - 1;
  score[0] = 0.0;
  score[1] = 0.0;

  // Early stopping on the number of patients at the current combination:
  // the design would not escalate (or cannot), and would not de-escalate
  // (or cannot, at the lowest combination).
  if (n[idx] >= n_earlystop) {
    const bool no_escalation = y[idx] > b_esc[ix];
    const bool at_highest = (a == J - 1) && (b == K - 1);
    const bool blocked_last_row =
      (a == J - 1) && (b < K - 1) && elim[a + J * (b + 1)] == 1;
    const bool blocked_last_col =
      (a < J - 1) && (b == K - 1) && elim[(a + 1) + J * b] == 1;
    const bool blocked_both =
      (a < J - 1) && (b < K - 1) &&
      elim[(a + 1) + J * b] == 1 && elim[a + J * (b + 1)] == 1;
    const bool no_deescalation = (y[idx] < b_deesc[ix]) || (a == 0 && b == 0);
    if ((no_escalation || at_highest || blocked_last_row ||
         blocked_last_col || blocked_both) && no_deescalation) {
      return COMB_ACTION_STOP;
    }
  }

  if (y[idx] <= b_esc[ix]) {
    // Candidates: one level up of drug A, then one level up of drug B.
    if (a + 1 < J && elim[(a + 1) + J * b] == 0) {
      const int c = (a + 1) + J * b;
      score[0] = comb_candidate_score(y[c], n[c], lambda_e, lambda_d);
    }
    if (b + 1 < K && elim[a + J * (b + 1)] == 0) {
      const int c = a + J * (b + 1);
      score[1] = comb_candidate_score(y[c], n[c], lambda_e, lambda_d);
    }
    return COMB_ACTION_ESCALATE;
  }
  if (y[idx] >= b_deesc[ix]) {
    // Candidates: one level down of drug A, then one level down of drug B.
    if (a - 1 >= 0) {
      const int c = (a - 1) + J * b;
      score[0] = comb_candidate_score(y[c], n[c], lambda_e, lambda_d);
    }
    if (b - 1 >= 0) {
      const int c = a + J * (b - 1);
      score[1] = comb_candidate_score(y[c], n[c], lambda_e, lambda_d);
    }
    return COMB_ACTION_DEESCALATE;
  }
  return COMB_ACTION_STAY;
}

// Simulate a single trial of the BOIN drug-combination design.
//
// rng() must return one uniform variate on (0, 1) per call. Variates are drawn
// in this order: during titration, one per patient and one for the direction
// of each move that has two possible directions; afterwards one per patient,
// in enrollment order, and one for each escalation or de-escalation that has
// at least one candidate. This reproduces the random number consumption of
// BOIN::get.oc.comb().
template <class Rng>
inline void comb_simulate_one(const std::vector<double>& p_true,
                              int n_rows,
                              int n_cols,
                              const std::vector<int>& cohort_size,
                              int start_row,
                              int start_col,
                              int n_earlystop,
                              bool titration,
                              bool extrasafe,
                              double target,
                              double cutoff_eli,
                              double offset,
                              double lambda_e,
                              double lambda_d,
                              const std::vector<int>& b_esc,
                              const std::vector<int>& b_deesc,
                              const std::vector<int>& b_elim,
                              Rng& rng,
                              std::vector<int>& n_pts,
                              std::vector<int>& n_tox,
                              std::vector<char>& eliminated,
                              int& cohorts_used,
                              int& stop_code) {

  const int J = n_rows;
  const int K = n_cols;
  const int n_cohort = static_cast<int>(cohort_size.size());
  const int max_n = static_cast<int>(b_esc.size());

  std::fill(n_pts.begin(), n_pts.end(), 0);
  std::fill(n_tox.begin(), n_tox.end(), 0);
  std::fill(eliminated.begin(), eliminated.end(), static_cast<char>(0));

  int* n = &n_pts[0];
  int* y = &n_tox[0];
  char* elim = &eliminated[0];

  int a = start_row;
  int b = start_col;
  stop_code = COMB_STOP_MAX_COHORTS;
  cohorts_used = 0;
  bool first_titration_cohort = true;

  // Titration: one patient per combination, moving up one level of either
  // drug at random, until the first DLT or the highest combination.
  if (titration) {
    int idx = a + J * b;
    y[idx] = (rng() < p_true[idx]) ? 1 : 0;
    n[idx] = 1;
    int n_dlt = y[idx];
    for (;;) {
      if (a == J - 1 && b == K - 1) break;
      if (n_dlt == 1) break;
      if (a < J - 1 && b < K - 1) {
        // rbinom(1, 1, 0.5) is 0 exactly when its single variate is below 0.5.
        if (rng() < 0.5) a = a + 1; else b = b + 1;
      } else if (a == J - 1) {
        b = b + 1;
      } else {
        a = a + 1;
      }
      idx = a + J * b;
      y[idx] = (rng() < p_true[idx]) ? 1 : 0;
      n[idx] = 1;
      n_dlt += y[idx];
    }
  }

  for (int i = 0; i < n_cohort; ++i) {
    cohorts_used = i + 1;
    const int cs = cohort_size[i];
    const int idx = a + J * b;

    if (titration && n[idx] < cs && first_titration_cohort) {
      // Fill the current combination up to the nominal cohort size.
      first_titration_cohort = false;
      const int add = cs - 1;
      int dlt = 0;
      for (int k = 0; k < add; ++k) {
        if (rng() < p_true[idx]) ++dlt;
      }
      y[idx] += dlt;
      n[idx] += add;
    } else {
      int dlt = 0;
      for (int k = 0; k < cs; ++k) {
        if (rng() < p_true[idx]) ++dlt;
      }
      y[idx] += dlt;
      n[idx] += cs;
    }

    int nc = n[idx];
    if (nc > max_n) nc = max_n;
    const int ix = nc - 1;

    // Elimination of the combination and every combination with at least the
    // same level of both drugs, and the extra safety rule, which is nested in
    // the elimination branch as in the reference implementation.
    if (b_elim[ix] > 0) {
      if (y[idx] >= b_elim[ix]) {
        for (int bb = b; bb < K; ++bb) {
          for (int aa = a; aa < J; ++aa) elim[aa + J * bb] = 1;
        }
        if (a == 0 && b == 0) {
          stop_code = COMB_STOP_LOWEST_ELIMINATED;
          break;
        }
      }
      if (extrasafe && a == 0 && b == 0 && n[0] >= 3) {
        if (post_prob_above(target, y[0], n[0]) > cutoff_eli - offset) {
          stop_code = COMB_STOP_LOWEST_TOO_TOXIC;
          break;
        }
      }
    }

    // Early stopping and the move to the next combination.
    double score[2] = {0.0, 0.0};
    const int action = comb_decide(n, y, elim, J, K, a, b, n_earlystop,
                                   lambda_e, lambda_d, b_esc, b_deesc, score);
    if (action == COMB_ACTION_STOP) {
      stop_code = COMB_STOP_N_EARLYSTOP;
      break;
    }
    if (action == COMB_ACTION_ESCALATE) {
      const int k = comb_pick(score, rng);
      if (k == 0) a = a + 1;
      if (k == 1) b = b + 1;
    } else if (action == COMB_ACTION_DEESCALATE) {
      const int k = comb_pick(score, rng);
      if (k == 0) a = a - 1;
      if (k == 1) b = b - 1;
    }
  }
}

// Elimination implied by the data of a trial: a combination whose DLT count
// has reached its elimination boundary eliminates itself and every combination
// with at least the same level of both drugs. During a simulated trial the
// combinations eliminated are exactly these, because elimination is checked
// after every cohort and an eliminated combination is never treated again.
inline void comb_data_elimination(const int* n, const int* y, int n_rows,
                                  int n_cols, const std::vector<int>& b_elim,
                                  char* elim) {
  const int J = n_rows;
  const int K = n_cols;
  const int max_n = static_cast<int>(b_elim.size());
  std::fill(elim, elim + J * K, static_cast<char>(0));
  for (int j = 0; j < K; ++j) {
    for (int i = 0; i < J; ++i) {
      const int idx = i + J * j;
      if (n[idx] <= 0) continue;
      int nc = n[idx];
      if (nc > max_n) nc = max_n;
      const int bound = b_elim[nc - 1];
      if (bound > 0 && y[idx] >= bound) {
        for (int bb = j; bb < K; ++bb) {
          for (int aa = i; aa < J; ++aa) elim[aa + J * bb] = 1;
        }
      }
    }
  }
}

// The decision for the next cohort during a trial, from the data so far and
// the combination (row a, column b, 0-based) that treated the last cohort,
// with the rules of comb_simulate_one(). Elimination is derived from the data
// with comb_data_elimination(). Returns a stop code, or -1 when the trial
// continues, in which case action and score describe the move as returned by
// comb_decide(); no random variate is drawn.
inline int comb_next(const int* n, const int* y, int n_rows, int n_cols,
                     int a, int b, int n_earlystop, bool extrasafe,
                     double target, double cutoff_eli, double offset,
                     double lambda_e, double lambda_d,
                     const std::vector<int>& b_esc,
                     const std::vector<int>& b_deesc,
                     const std::vector<int>& b_elim,
                     char* elim, int& action, double* score) {
  const int J = n_rows;
  const int K = n_cols;
  const int max_n = static_cast<int>(b_esc.size());
  action = COMB_ACTION_STAY;
  score[0] = 0.0;
  score[1] = 0.0;

  comb_data_elimination(n, y, J, K, b_elim, elim);
  if (elim[0] == 1) return COMB_STOP_LOWEST_ELIMINATED;

  int nc = n[a + J * b];
  if (nc > max_n) nc = max_n;
  const int ix = nc - 1;
  if (a == 0 && b == 0 && b_elim[ix] > 0 && extrasafe && n[0] >= 3 &&
      post_prob_above(target, y[0], n[0]) > cutoff_eli - offset) {
    return COMB_STOP_LOWEST_TOO_TOXIC;
  }

  action = comb_decide(n, y, elim, J, K, a, b, n_earlystop, lambda_e,
                       lambda_d, b_esc, b_deesc, score);
  if (action == COMB_ACTION_STOP) return COMB_STOP_N_EARLYSTOP;
  return -1;
}

// Elimination used for the MTD selection, re-derived from the final data as
// the reference implementation does. Within each row, the first combination
// (lowest level of drug B) that meets the elimination criterion eliminates
// itself, the rest of its row and the rest of its column; the scan then moves
// to the next row. The extra safety rule eliminates every combination.
inline void comb_selection_elimination(const int* n, const int* y,
                                       int n_rows, int n_cols,
                                       double target, double cutoff_eli,
                                       bool extrasafe, double offset,
                                       char* elim) {
  const int J = n_rows;
  const int K = n_cols;
  std::fill(elim, elim + J * K, static_cast<char>(0));

  if (extrasafe && n[0] >= 3 &&
      post_prob_above(target, y[0], n[0]) > cutoff_eli - offset) {
    std::fill(elim, elim + J * K, static_cast<char>(1));
  }

  for (int i = 0; i < J; ++i) {
    for (int j = 0; j < K; ++j) {
      const int idx = i + J * j;
      if (n[idx] >= 3 && post_prob_above(target, y[idx], n[idx]) > cutoff_eli) {
        for (int ii = i; ii < J; ++ii) elim[ii + J * j] = 1;
        for (int jj = j; jj < K; ++jj) elim[i + J * jj] = 1;
        break;
      }
    }
  }
}

// Choose the MTD from the bivariate isotonic estimates of a single trial.
//
// iso holds the isotonic estimates of (y + 0.05) / (n + 0.1) with weights
// n + 0.1, computed by the caller. Eliminated combinations get the value 1.1,
// untreated ones are excluded, and 1e-5 * (row + column) (1-based) is added to
// break ties in favor of lower combinations. The combination closest to the
// target is selected; among exact ties the first in column-major order wins.
// Returns the reason code and writes the 1-based row and column, or 0 and 0
// when no combination is selected.
inline int comb_choose_mtd(const double* iso, const int* n, const char* elim,
                           int n_rows, int n_cols, double target,
                           bool bound_mtd, double lambda_d,
                           int& sel_row, int& sel_col) {
  const int J = n_rows;
  const int K = n_cols;
  sel_row = 0;
  sel_col = 0;

  if (elim[0] == 1) return COMB_MTD_LOWEST_ELIMINATED;

  std::vector<double> phat(J * K);
  for (int j = 0; j < K; ++j) {
    for (int i = 0; i < J; ++i) {
      const int idx = i + J * j;
      const double value = (elim[idx] == 1) ? 1.1 : iso[idx];
      const double masked = value * ((n[idx] != 0) ? 1.0 : 0.0);
      const double shift = 1e-5 * static_cast<double>((i + 1) + (j + 1));
      phat[idx] = masked + shift;
    }
  }

  if (bound_mtd) {
    bool all_above = true;
    for (int idx = 0; idx < J * K; ++idx) {
      if (n[idx] != 0 && !(phat[idx] >= lambda_d)) { all_above = false; break; }
    }
    if (all_above) return COMB_MTD_NO_DOSE_BELOW_LAMBDA_D;
    for (int idx = 0; idx < J * K; ++idx) {
      if (phat[idx] >= lambda_d) phat[idx] = 10.0;
    }
  }

  for (int idx = 0; idx < J * K; ++idx) {
    if (n[idx] == 0) phat[idx] = 10.0;
  }

  double best = std::fabs(phat[0] - target);
  for (int idx = 1; idx < J * K; ++idx) {
    const double dist = std::fabs(phat[idx] - target);
    if (dist < best) best = dist;
  }
  for (int idx = 0; idx < J * K; ++idx) {
    if (std::fabs(phat[idx] - target) == best) {
      sel_row = idx % J + 1;
      sel_col = idx / J + 1;
      break;
    }
  }
  return COMB_MTD_SELECTED;
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_COMB_CORE_H
