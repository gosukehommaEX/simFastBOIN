#ifndef SIMFASTBOIN_BOIN_CORE_H
#define SIMFASTBOIN_BOIN_CORE_H

// Core algorithms of the BOIN design, written in plain C++ so that the same
// translation unit can be compiled both inside the package (against R's
// numerical routines) and inside a standalone verification harness.
//
// The only external dependency is post_prob_above(), which the caller supplies.
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

// Stop codes returned by simulate_one().
enum StopCode {
  STOP_LOWEST_ELIMINATED = 0,
  STOP_LOWEST_TOO_TOXIC = 1,
  STOP_N_EARLYSTOP = 2,
  STOP_MAX_SAMPLE_SIZE = 3,
  STOP_MAX_COHORTS = 4
};

// Reason codes returned by select_mtd_one().
enum MtdReason {
  MTD_SELECTED = 0,
  MTD_LOWEST_ELIMINATED = 1,
  MTD_NO_ADMISSIBLE_DOSE = 2,
  MTD_NO_DOSE_BELOW_LAMBDA_D = 3
};

// Pool adjacent violators algorithm for weighted isotonic regression.
//
// The pooling loop reproduces, step for step, the implementation used by the
// reference BOIN package: repeatedly pool the leftmost pair of level sets that
// violates monotonicity, recomputing the pooled value as a weighted mean over
// the individual elements in index order. The result is the unique weighted
// isotonic regression estimate; keeping the accumulation order identical makes
// the output bit-for-bit reproducible against the reference.
inline void pava(std::vector<double>& x, const std::vector<double>& w) {
  const int n = static_cast<int>(x.size());
  if (n <= 1) return;

  std::vector<int> lvl(n);
  for (int i = 0; i < n; ++i) lvl[i] = i;

  for (;;) {
    int viol = -1;
    for (int k = 0; k + 1 < n; ++k) {
      if (x[k + 1] - x[k] < 0.0) { viol = k; break; }
    }
    if (viol < 0) break;

    const int lvl1 = lvl[viol];
    const int lvl2 = lvl[viol + 1];

    double sum_xw = 0.0;
    double sum_w = 0.0;
    for (int k = 0; k < n; ++k) {
      if (lvl[k] == lvl1 || lvl[k] == lvl2) {
        sum_xw += x[k] * w[k];
        sum_w += w[k];
      }
    }
    const double pooled = sum_xw / sum_w;
    for (int k = 0; k < n; ++k) {
      if (lvl[k] == lvl1 || lvl[k] == lvl2) {
        x[k] = pooled;
        lvl[k] = lvl1;
      }
    }
  }
}

// Posterior mean estimate and inverse-variance weight used for the isotonic fit.
inline void phat_and_weight(int y, int n, double& phat, double& weight) {
  const double yy = y + 0.05;
  const double nn = n + 0.1;
  const double zz = n - y + 0.05;
  phat = yy / nn;
  weight = 1.0 / (yy * zz / (nn * nn * (nn + 1.0)));
}

// Simulate a single BOIN trial.
//
// rng() must return one uniform variate on (0, 1) per call. Exactly one variate
// is consumed per patient, drawn in enrollment order, and the titration phase
// draws one variate per dose level unconditionally. This reproduces the random
// number consumption of the reference implementation.
template <class Rng>
inline void simulate_one(const std::vector<double>& p_true,
                         const std::vector<int>& cohort_size,
                         int start_dose,
                         int n_earlystop,
                         bool early_stop_simple,
                         bool titration,
                         bool extrasafe,
                         double target,
                         double cutoff_eli,
                         double offset,
                         const std::vector<int>& b_esc,
                         const std::vector<int>& b_deesc,
                         const std::vector<int>& b_elim,
                         int max_total_pts,
                         Rng& rng,
                         std::vector<int>& n_pts,
                         std::vector<int>& n_tox,
                         std::vector<char>& eliminated,
                         int& cohorts_used,
                         int& stop_code) {

  const int n_doses = static_cast<int>(p_true.size());
  const int n_cohort = static_cast<int>(cohort_size.size());
  const int max_n = static_cast<int>(b_esc.size());

  std::fill(n_pts.begin(), n_pts.end(), 0);
  std::fill(n_tox.begin(), n_tox.end(), 0);
  std::fill(eliminated.begin(), eliminated.end(), static_cast<char>(0));

  int* n = &n_pts[0];
  int* y = &n_tox[0];
  char* elim = &eliminated[0];

  int d = start_dose - 1;
  int total = 0;
  stop_code = STOP_MAX_COHORTS;
  cohorts_used = 0;
  bool first_titration_cohort = true;

  int max_cs = 1;
  for (int i = 0; i < n_cohort; ++i) max_cs = std::max(max_cs, cohort_size[i]);
  std::vector<int> cohort_dlt(max_cs, 0);

  // Titration phase: one patient per dose level until the first DLT.
  if (titration) {
    int first_dlt = -1;
    for (int j = 0; j < n_doses; ++j) {
      const double u = rng();
      if (first_dlt < 0 && u < p_true[j]) first_dlt = j;
    }
    if (first_dlt < 0) {
      d = n_doses - 1;
      for (int j = 0; j < n_doses; ++j) n[j] = 1;
      total = n_doses;
    } else {
      d = first_dlt;
      for (int j = 0; j <= first_dlt; ++j) n[j] = 1;
      y[first_dlt] = 1;
      total = first_dlt + 1;
    }
  }

  for (int i = 0; i < n_cohort; ++i) {
    cohorts_used = i + 1;
    const int cs = cohort_size[i];

    if (titration && n[d] < cs && first_titration_cohort) {

      // Fill the current dose up to the nominal cohort size.
      first_titration_cohort = false;
      const int add = cs - 1;
      int dlt = 0;
      for (int k = 0; k < add; ++k) {
        if (rng() < p_true[d]) ++dlt;
      }
      y[d] += dlt;
      n[d] += add;
      total += add;

    } else {

      for (int k = 0; k < cs; ++k) {
        cohort_dlt[k] = (rng() < p_true[d]) ? 1 : 0;
      }

      if (total + cs >= max_total_pts) {
        int n_remain = max_total_pts - total;
        if (n_remain < 0) n_remain = 0;
        if (n_remain > cs) n_remain = cs;
        int dlt = 0;
        for (int k = 0; k < n_remain; ++k) dlt += cohort_dlt[k];
        y[d] += dlt;
        n[d] += n_remain;
        total += n_remain;
        stop_code = STOP_MAX_SAMPLE_SIZE;
        break;
      }

      int dlt = 0;
      for (int k = 0; k < cs; ++k) dlt += cohort_dlt[k];
      y[d] += dlt;
      n[d] += cs;
      total += cs;
    }

    int nd = n[d];
    if (nd > max_n) nd = max_n;
    const int idx = nd - 1;

    // Dose elimination and the extra safety stopping rule.
    // The safety rule is nested inside the elimination branch. This mirrors the
    // reference implementation: when no elimination boundary exists at the
    // current sample size, the safety rule is not evaluated either.
    if (b_elim[idx] > 0) {
      if (y[d] >= b_elim[idx]) {
        for (int j = d; j < n_doses; ++j) elim[j] = 1;
        if (d == 0) { stop_code = STOP_LOWEST_ELIMINATED; break; }
      }
      if (extrasafe && d == 0 && n[0] >= 3) {
        if (post_prob_above(target, y[0], n[0]) > cutoff_eli - offset) {
          stop_code = STOP_LOWEST_TOO_TOXIC;
          break;
        }
      }
    }

    // Early stopping on the number of patients at the current dose. This is
    // evaluated before the dose transition, at the dose that was just treated.
    bool stop_now = false;
    if (n[d] >= n_earlystop) {
      if (early_stop_simple) {
        stop_now = true;
      } else {
        const int be = b_esc[idx];
        const int bd = b_deesc[idx];
        const bool stay = (y[d] > be) && (y[d] < bd);
        const bool at_lowest = (d == 0) && (y[d] >= bd);
        const bool at_highest =
          ((d == n_doses - 1) || (elim[d + 1] == 1)) && (y[d] <= be);
        stop_now = stay || at_lowest || at_highest;
      }
    }
    if (stop_now) { stop_code = STOP_N_EARLYSTOP; break; }

    // Dose transition.
    if (y[d] <= b_esc[idx] && d != n_doses - 1) {
      if (elim[d + 1] == 0) d = d + 1;
    } else if (y[d] >= b_deesc[idx] && d != 0) {
      d = d - 1;
    }
  }
}

// Select the MTD from the final data of a single trial.
//
// Dose elimination is re-derived from the final data rather than carried over
// from the trial, which is what the reference implementation does.
// Returns the selected dose (1-based) or 0 when no dose is selected.
inline int select_mtd_one(const int* n, const int* y, int n_doses,
                          double target, double cutoff_eli,
                          bool extrasafe, double offset,
                          bool bound_mtd, double lambda_d,
                          int min_mtd_sample, int& reason) {

  std::vector<char> elim(n_doses, 0);

  for (int j = 0; j < n_doses; ++j) {
    if (n[j] >= 3 && post_prob_above(target, y[j], n[j]) > cutoff_eli) {
      for (int k = j; k < n_doses; ++k) elim[k] = 1;
      break;
    }
  }
  if (extrasafe && n[0] >= 3 &&
      post_prob_above(target, y[0], n[0]) > cutoff_eli - offset) {
    std::fill(elim.begin(), elim.end(), static_cast<char>(1));
  }

  if (elim[0] == 1) {
    reason = MTD_LOWEST_ELIMINATED;
    return 0;
  }

  std::vector<int> adm;
  for (int j = 0; j < n_doses; ++j) {
    if (elim[j] == 0 && n[j] > 0 && n[j] >= min_mtd_sample) adm.push_back(j);
  }
  if (adm.empty()) {
    reason = MTD_NO_ADMISSIBLE_DOSE;
    return 0;
  }

  const int k_adm = static_cast<int>(adm.size());
  std::vector<double> phat(k_adm), weight(k_adm);
  for (int k = 0; k < k_adm; ++k) {
    phat_and_weight(y[adm[k]], n[adm[k]], phat[k], weight[k]);
  }
  pava(phat, weight);
  for (int k = 0; k < k_adm; ++k) phat[k] += (k + 1) * 1e-10;

  int n_keep = k_adm;
  if (bound_mtd) {
    // The isotonic estimates are non-decreasing, so the doses satisfying the
    // constraint always form a prefix of the admissible set.
    n_keep = 0;
    for (int k = 0; k < k_adm; ++k) {
      if (phat[k] <= lambda_d) ++n_keep; else break;
    }
    if (n_keep == 0) {
      reason = MTD_NO_DOSE_BELOW_LAMBDA_D;
      return 0;
    }
  }

  int best = 0;
  double best_dist = std::fabs(phat[0] - target);
  for (int k = 1; k < n_keep; ++k) {
    const double dist = std::fabs(phat[k] - target);
    if (dist < best_dist) { best_dist = dist; best = k; }
  }

  reason = MTD_SELECTED;
  return adm[best] + 1;
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_BOIN_CORE_H
