#ifndef SIMFASTBOIN_UBOIN_CORE_H
#define SIMFASTBOIN_UBOIN_CORE_H

// Core algorithms of the U-BOIN design of Zhou, Lee and Yuan (2019) for a
// binary toxicity and a binary efficacy outcome observed before the next
// cohort is treated, with the pick-the-winner allocation. The code is plain
// C++ so that it can be compiled both inside the package and inside a
// standalone verification harness.
//
// The outcome of a patient is one of four categories, in the order of the
// article: 0 (no response, DLT), 1 (no response, no DLT), 2 (response, DLT)
// and 3 (response, no DLT). Counts are stored per dose, nk[4 * j + k].
//
// Every posterior probability that decides whether a dose is eliminated or
// admissible is evaluated in R and passed in as a table, indexed by the
// number of patients n (0 to M) and the number of events c (0 to M):
//   elim[n * (M + 1) + c]    1 when c DLTs in n patients eliminate the dose
//                            in stage I (only for n >= 3);
//   toxic[n * (M + 1) + c]   1 when c DLTs in n patients make it toxic;
//   futile[n * (M + 1) + c]  1 when c responses in n patients make it futile.
// Boundary vectors are indexed by the number of patients minus one.

#include <vector>
#include <algorithm>
#include <cmath>

namespace simfastboin {

// Stop codes returned by uboin_simulate_one().
enum UbStop {
  UB_STOP_MAX_SAMPLE_SIZE = 0,
  UB_STOP_N_STOP = 1,
  UB_STOP_LOWEST_ELIMINATED = 2,
  UB_STOP_NO_ADMISSIBLE = 3
};

// The settings of a U-BOIN trial that do not change from trial to trial.
struct UbDesign {
  int n_doses;
  std::vector<double> p_outcome;   // 4 * n_doses, category probabilities
  int n_cohort;
  int cohort_size;
  int start_dose;                  // 0-based
  int n_stage1;                    // s1
  int n_stop;                      // s2
  std::vector<double> utility;     // 4 values, one per category
  double prior_a;                  // Dirichlet parameter of each category
  int max_n;                       // M, the largest number at one dose
  std::vector<int> elim;
  std::vector<int> toxic;
  std::vector<int> futile;
  std::vector<int> b_esc;
  std::vector<int> b_deesc;
};

// Rounds to ten decimals. Utilities are compared after this rounding, so
// that doses whose utilities are equal as fractions tie whatever the
// rounding error of the platform.
inline double ub_round10(double x) {
  return std::floor(x * 1e10 + 0.5) / 1e10;
}

// Posterior mean utility of a dose with the counts nk[0..3]: the utilities
// weighted by the Dirichlet posterior means (prior_a + n_k) / (4 prior_a + n).
inline double ub_utility(const int* nk, const UbDesign& des) {
  double num = 0.0;
  int n = 0;
  for (int k = 0; k < 4; ++k) {
    num += des.utility[k] * (des.prior_a + nk[k]);
    n += nk[k];
  }
  return num / (4.0 * des.prior_a + n);
}

inline int ub_index(const UbDesign& des, int n, int c) {
  return n * (des.max_n + 1) + c;
}

// Whether a dose with the counts nk[0..3] is admissible: neither toxic nor
// futile.
inline bool ub_admissible(const int* nk, const UbDesign& des) {
  const int n = nk[0] + nk[1] + nk[2] + nk[3];
  const int dlt = nk[0] + nk[2];
  const int resp = nk[2] + nk[3];
  return des.toxic[ub_index(des, n, dlt)] == 0 &&
    des.futile[ub_index(des, n, resp)] == 0;
}

// The dose (0-based) with the largest posterior mean utility among the doses
// that have been tried, are not eliminated and are admissible, the lowest of
// them among ties; -1 when there is none.
inline int ub_best_dose(const int* nk, const char* eliminated,
                        const UbDesign& des) {
  int best = -1;
  double best_u = 0.0;
  for (int j = 0; j < des.n_doses; ++j) {
    const int* c = nk + 4 * j;
    const int n = c[0] + c[1] + c[2] + c[3];
    if (n == 0 || eliminated[j] == 1) continue;
    if (!ub_admissible(c, des)) continue;
    const double u = ub_round10(ub_utility(c, des));
    if (best < 0 || u > best_u) {
      best = j;
      best_u = u;
    }
  }
  return best;
}

// The category (0 to 3) of a patient at dose d from one uniform variate u:
// the first category whose cumulative probability exceeds u.
inline int ub_category(const UbDesign& des, int d, double u) {
  const double* p = &des.p_outcome[4 * d];
  double cum = 0.0;
  for (int k = 0; k < 3; ++k) {
    cum += p[k];
    if (u < cum) return k;
  }
  return 3;
}

// Simulate a single U-BOIN trial.
//
// rng() must return one uniform variate on (0, 1) per call; one is drawn per
// patient, in enrollment order. nk receives the counts of every dose and
// category, eliminated the doses eliminated in stage I.
//
// Stage I is the BOIN design: after each cohort, a dose whose DLTs reach the
// elimination table is eliminated together with every higher dose (the trial
// stops when it is the lowest), and the next cohort escalates (unless the
// next dose is eliminated), stays or de-escalates by the boundaries. Once a
// dose has n_stage1 patients, the next cohorts follow stage II: rule B1
// escalates to the dose above the highest dose tried while the DLT rate
// observed there is at most the escalation boundary (and that dose is not
// eliminated); otherwise rule B2 assigns the dose of ub_best_dose(), and the
// trial stops when there is none. The trial ends at the maximum sample size
// or when a dose has n_stop patients.
template <class Rng>
inline void uboin_simulate_one(const UbDesign& des, Rng& rng,
                               std::vector<int>& nk,
                               std::vector<char>& eliminated,
                               int& cohorts_used, int& stage2,
                               int& stop_code) {
  const int J = des.n_doses;
  std::fill(nk.begin(), nk.end(), 0);
  std::fill(eliminated.begin(), eliminated.end(), static_cast<char>(0));
  cohorts_used = 0;
  stage2 = 0;
  stop_code = UB_STOP_MAX_SAMPLE_SIZE;

  int d = des.start_dose;
  for (int cohort = 0; cohort < des.n_cohort; ++cohort) {
    cohorts_used = cohort + 1;
    for (int i = 0; i < des.cohort_size; ++i) {
      const int k = ub_category(des, d, rng());
      ++nk[4 * d + k];
    }
    const int* c = &nk[4 * d];
    const int n = c[0] + c[1] + c[2] + c[3];
    const int dlt = c[0] + c[2];

    if (stage2 == 0 && des.elim[ub_index(des, n, dlt)] == 1) {
      for (int j = d; j < J; ++j) eliminated[j] = 1;
      if (d == 0) {
        stop_code = UB_STOP_LOWEST_ELIMINATED;
        return;
      }
    }

    if (cohort == des.n_cohort - 1) break;
    int largest = 0;
    for (int j = 0; j < J; ++j) {
      const int nj = nk[4 * j] + nk[4 * j + 1] + nk[4 * j + 2] +
        nk[4 * j + 3];
      if (nj > largest) largest = nj;
    }
    if (largest >= des.n_stop) {
      stop_code = UB_STOP_N_STOP;
      return;
    }

    if (stage2 == 0) {
      if (dlt <= des.b_esc[n - 1] && d < J - 1 && eliminated[d + 1] == 0) {
        d = d + 1;
      } else if (dlt >= des.b_deesc[n - 1] && d > 0) {
        d = d - 1;
      }
      if (largest < des.n_stage1) continue;
      stage2 = 1;
    }

    // Stage II, rule B1.
    int jstar = -1;
    for (int j = 0; j < J; ++j) {
      if (nk[4 * j] + nk[4 * j + 1] + nk[4 * j + 2] + nk[4 * j + 3] > 0) {
        jstar = j;
      }
    }
    const int* cs = &nk[4 * jstar];
    const int ns = cs[0] + cs[1] + cs[2] + cs[3];
    const int dlts = cs[0] + cs[2];
    if (jstar < J - 1 && dlts <= des.b_esc[ns - 1] &&
        eliminated[jstar + 1] == 0) {
      d = jstar + 1;
      continue;
    }

    // Rule B2.
    const int best = ub_best_dose(&nk[0], &eliminated[0], des);
    if (best < 0) {
      stop_code = UB_STOP_NO_ADMISSIBLE;
      return;
    }
    d = best;
  }
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_UBOIN_CORE_H
