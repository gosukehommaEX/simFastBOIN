#ifndef SIMFASTBOIN_IBOIN_CORE_H
#define SIMFASTBOIN_IBOIN_CORE_H

// Core algorithm of the iBOIN design of Zhou et al. (2021), written in plain
// C++ so that it can be compiled inside the package and inside a standalone
// verification harness.
//
// The trial runs as a BOIN trial without titration, see simulate_one() in
// boin_core.h, but the escalation and de-escalation boundaries depend on the
// dose as well as on the number of patients treated there. The boundaries
// b_esc and b_deesc hold one block of max_n values per dose, the value for
// dose j (0-based) and n patients at index j * max_n + n - 1. The elimination
// boundary b_elim does not depend on the dose and is indexed by the number
// of patients minus one, 0 meaning that no boundary exists.

#include "boin_core.h"

namespace simfastboin {

// Simulate a single trial whose boundaries depend on the dose.
//
// rng() must return one uniform variate on (0, 1) per call. Exactly one
// variate is consumed per patient, drawn in enrollment order, and a cohort
// cut short by the maximum sample size still draws a full cohort's worth, as
// in simulate_one(). With the same boundaries at every dose the trials are
// those of simulate_one() without titration.
template <class Rng>
inline void simulate_one_by_dose(const std::vector<double>& p_true,
                                 const std::vector<int>& cohort_size,
                                 int start_dose,
                                 int n_earlystop,
                                 bool early_stop_simple,
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
  const int max_n = static_cast<int>(b_elim.size());

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

  int max_cs = 1;
  for (int i = 0; i < n_cohort; ++i) max_cs = std::max(max_cs, cohort_size[i]);
  std::vector<int> cohort_dlt(max_cs, 0);

  for (int i = 0; i < n_cohort; ++i) {
    cohorts_used = i + 1;
    const int cs = cohort_size[i];

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

    int nd = n[d];
    if (nd > max_n) nd = max_n;
    const int idx = nd - 1;
    const int be = b_esc[d * max_n + idx];
    const int bd = b_deesc[d * max_n + idx];

    // Dose elimination and the extra safety stopping rule, both under the
    // uniform prior, nested as in simulate_one().
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

    // Early stopping on the number of patients at the current dose.
    bool stop_now = false;
    if (n[d] >= n_earlystop) {
      if (early_stop_simple) {
        stop_now = true;
      } else {
        const bool stay = (y[d] > be) && (y[d] < bd);
        const bool at_lowest = (d == 0) && (y[d] >= bd);
        const bool at_highest =
          ((d == n_doses - 1) || (elim[d + 1] == 1)) && (y[d] <= be);
        stop_now = stay || at_lowest || at_highest;
      }
    }
    if (stop_now) { stop_code = STOP_N_EARLYSTOP; break; }

    // Dose transition.
    if (y[d] <= be && d != n_doses - 1) {
      if (elim[d + 1] == 0) d = d + 1;
    } else if (y[d] >= bd && d != 0) {
      d = d - 1;
    }
  }
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_IBOIN_CORE_H
