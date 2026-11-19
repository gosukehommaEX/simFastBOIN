#ifndef SIMFASTBOIN_BOIN_ETC_CORE_H
#define SIMFASTBOIN_BOIN_ETC_CORE_H

// Core algorithms of the BOIN-ETC design (Kakizume, Takeda, Taguri and
// Morita, 2024) for two drugs given together, written in plain C++ like
// boinet_core.h so that they can be verified outside of R as well.
//
// A dose matrix with J levels of drug A (rows a) and K levels of drug B
// (columns b) is stored column-major: the combination (a, b), 0-based, is
// element c = a + J * b. J and K are at least 2.
//
// The matrix is divided into subtrials, each a path along which the DLT
// probability increases:
//   - the first subtrial S_L: (0, 0), (1, 0), ..., (J - 1, 0), (J - 1, 1),
//     ..., (J - 1, K - 1), started at (0, 0);
//   - ETC_ETC1 and ETC_ETC2: then the rows r = J - 2, ..., 0, each the
//     combinations (r, 1), ..., (r, K - 1). When the candidate of S_L is in
//     its first column, (a*, 0), the rows start at max(a* - 1, 0) and the
//     rows above are skipped. ETC_ETC2 starts every row at (r, 1); ETC_ETC1
//     at the column of the previous candidate (ETC_START_BELOW) or one column
//     to its left (ETC_START_DIAGONAL), but not before column 1, and at
//     (r, 1) after a candidate in column 0 or a subtrial without candidate;
//   - ETC_ETC3: then the L-shapes j = 1, ..., min(J, K) - 1, each column j
//     from (0, j) to (J - 1 - j, j) followed by row J - 1 - j from column
//     j + 1 to K - 1, started at (0, j).
//
// Within a subtrial, after each cohort at its dose d (local index i):
//   - a dose whose DLT data fail (table tox_ok) is excluded together with
//     every higher dose of the subtrial; when d is the lowest dose of the
//     subtrial, the subtrial ends without a candidate (the trial ends when
//     this happens in S_L);
//   - the decision of the BOIN-ET table dec at d: escalation to i + 1
//     (staying at the highest dose and below an excluded dose), stay,
//     de-escalation to i - 1 (staying at the lowest dose), or rule 4;
//   - rule 4 (ETC_RULE4_PAPER): i + 1 when it has not been treated;
//     otherwise the treated, not excluded dose of {i - 1, i, i + 1} with the
//     highest response rate, ties broken with one uniform variate u, the
//     dose floor(u m) of the m tied doses in increasing order;
//   - rule 4 (ETC_RULE4_AUTHORS), the rule of the code of the authors
//     (FUNC_BOINETC.subtrial.R): i + 1 when it has not been treated;
//     otherwise, with D the doses of {i - 1, i, i + 1} inside the subtrial
//     that are not excluded and R the response rates of the treated ones of
//     D, in order, the position p of the highest rate in R selects the dose
//     D[p] (when i - 1 is untreated the choice moves one dose down); when
//     the highest rate occurs nn > 1 times, one uniform variate u selects
//     D[ceil(u nn) - 1];
//   - the subtrial ends when the next dose has n_stop patients (that dose is
//     the candidate) and the trial ends when all cohorts are enrolled.
//
// Each patient draws two uniform variates in enrollment order, the first for
// the DLT and the second for the response (as simulate_et_one()); a tie of
// rule 4 draws one more after the cohort.
//
// Selection at the end of a trial (etc_select_one()):
//   - iso: unweighted least-squares fit of (y + 0.05) / (n + 0.1),
//     non-decreasing in a and in b (Dykstra's algorithm until no estimate
//     changes by ETC_ISO_TOL in a cycle);
//   - in every row a, the MTD column is the highest column whose
//     |iso - target| is within ETC_MTD_TOL of the smallest in the row (tried
//     or not);
//   - the candidates are the tried, not excluded combinations (a, b) with b
//     at most the MTD column of row a; the utility estimate is
//     100 (q + alpha) / (n + alpha + beta), with the quasi-binomial number of
//     events q = (100 x + psi00 (n - y - x) + both_weight z) / 100;
//   - the OBD is the candidate with the highest estimate; all candidates
//     within ETC_UTIL_TOL of it are tied, and the one with the lowest a + b
//     and then the first in column-major order is returned.

#include <vector>
#include <cmath>
#include <algorithm>
#include "boinet_core.h"
#include "comb_boin12_core.h"

namespace simfastboin {

enum EtcSubtrials {
  ETC_ETC1 = 0,
  ETC_ETC2 = 1,
  ETC_ETC3 = 2
};

enum EtcStartRule {
  ETC_START_DIAGONAL = 0,
  ETC_START_BELOW = 1
};

enum EtcRule4 {
  ETC_RULE4_PAPER = 0,
  ETC_RULE4_AUTHORS = 1
};

// Stop codes returned by etc_simulate_one().
enum EtcStop {
  ETC_STOP_LOWEST_EXCLUDED = 0,  // the lowest combination (0, 0) is excluded
  ETC_STOP_COMPLETED = 1,        // every subtrial was run
  ETC_STOP_MAX_COHORTS = 2       // all cohorts were enrolled
};

const double ETC_MTD_TOL = 1e-8;
const double ETC_UTIL_TOL = 1e-6;
const double ETC_ISO_TOL = 1e-13;
const int ETC_ISO_MAX_CYCLES = 100000;

// The settings of a BOIN-ETC trial that do not change from trial to trial.
struct EtcDesign {
  int n_rows;                      // J
  int n_cols;                      // K
  std::vector<double> p_true;      // column-major
  std::vector<double> p_resp;
  std::vector<double> q0;          // Pr(response | no DLT), binary link only
  std::vector<double> q1;          // Pr(response | DLT), binary link only
  int cor_type;                    // EtCorType
  double resp_cor;                 // copula correlation
  std::vector<int> cohort_size;    // one entry per cohort
  int subtrials;                   // EtcSubtrials
  int start_rule;                  // EtcStartRule (ETC_ETC1 only)
  int rule4;                       // EtcRule4
  int n_stop;                      // patients at the next dose that end a
                                   // subtrial
  int max_n;                       // largest number of patients in the tables
  std::vector<int> dec;            // EtDecision table, as in EtDesign
  std::vector<int> tox_ok;         // tox_ok[n * (max_n + 1) + y]
};

// Combinations of the subtrials (column-major indices).
inline std::vector<int> etc_first_subtrial(int J, int K) {
  std::vector<int> out;
  for (int a = 0; a < J; ++a) out.push_back(a);
  for (int b = 1; b < K; ++b) out.push_back(J - 1 + J * b);
  return out;
}

inline std::vector<int> etc_row_subtrial(int J, int K, int r) {
  std::vector<int> out;
  for (int b = 1; b < K; ++b) out.push_back(r + J * b);
  return out;
}

inline std::vector<int> etc_lshape_subtrial(int J, int K, int j) {
  std::vector<int> out;
  for (int a = 0; a <= J - 1 - j; ++a) out.push_back(a + J * j);
  for (int b = j + 1; b < K; ++b) out.push_back(J - 1 - j + J * b);
  return out;
}

// One patient at combination c: two uniform variates, DLT and response.
template <class Rng>
inline void etc_patient(const EtcDesign& des, int c, Rng& rng, bool& dlt,
                        bool& resp) {
  const double u = rng();
  const double v = rng();
  dlt = u < des.p_true[c];
  if (des.cor_type == ET_COR_BINARY) {
    resp = v < (dlt ? des.q1[c] : des.q0[c]);
  } else {
    resp = response_variate(u, v, des.resp_cor) < des.p_resp[c];
  }
}

// Rule 4 at local dose i of a subtrial with doses dose[0..m-1].
template <class Rng>
inline int etc_rule4(int i, const std::vector<int>& dose,
                     const std::vector<char>& excl, const int* n,
                     const int* x, int rule, Rng& rng) {
  const int m = static_cast<int>(dose.size());
  if (rule == ETC_RULE4_AUTHORS) {
    if (i + 1 < m && n[dose[i + 1]] == 0) return i + 1;
    std::vector<int> dd;
    if (i == 0) {
      dd.push_back(0);
      if (m > 1) dd.push_back(1);
    } else if (i == m - 1) {
      dd.push_back(i - 1);
      dd.push_back(i);
    } else {
      dd.push_back(i - 1);
      dd.push_back(i);
      dd.push_back(i + 1);
    }
    std::vector<int> elig;
    std::vector<double> rate;
    for (std::size_t k = 0; k < dd.size(); ++k) {
      if (excl[dd[k]]) continue;
      elig.push_back(dd[k]);
      const int c = dose[dd[k]];
      if (n[c] > 0) rate.push_back(static_cast<double>(x[c]) / n[c]);
    }
    if (rate.empty()) return i;
    const double best = *std::max_element(rate.begin(), rate.end());
    std::vector<int> pos;
    for (std::size_t k = 0; k < rate.size(); ++k) {
      if (rate[k] == best) pos.push_back(static_cast<int>(k));
    }
    if (pos.size() == 1) return elig[pos[0]];
    const int nn = static_cast<int>(pos.size());
    const double u = rng();
    int pick = static_cast<int>(std::ceil(u * nn)) - 1;
    if (pick < 0) pick = 0;
    if (pick >= nn) pick = nn - 1;
    return elig[pick];
  }
  if (i + 1 < m && n[dose[i + 1]] == 0 && !excl[i + 1]) return i + 1;
  std::vector<int> cand;
  std::vector<double> rate;
  for (int k = i - 1; k <= i + 1; ++k) {
    if (k < 0 || k >= m || excl[k] || n[dose[k]] == 0) continue;
    cand.push_back(k);
    rate.push_back(static_cast<double>(x[dose[k]]) / n[dose[k]]);
  }
  if (cand.empty()) return i;
  const double best = *std::max_element(rate.begin(), rate.end());
  std::vector<int> ties;
  for (std::size_t k = 0; k < cand.size(); ++k) {
    if (rate[k] == best) ties.push_back(cand[k]);
  }
  if (ties.size() == 1) return ties[0];
  const int mt = static_cast<int>(ties.size());
  const double u = rng();
  int pick = static_cast<int>(u * mt);
  if (pick >= mt) pick = mt - 1;
  return ties[pick];
}

// Run one subtrial from local dose start. cohort is the index of the next
// cohort and is advanced. Returns the candidate (column-major index), or -1
// when the lowest dose of the subtrial was excluded; full is set when all
// cohorts have been enrolled.
template <class Rng>
inline int etc_run_subtrial(const EtcDesign& des,
                            const std::vector<int>& dose, int start,
                            Rng& rng, std::vector<int>& n, std::vector<int>& y,
                            std::vector<int>& x, std::vector<int>& z,
                            std::vector<char>& excluded, int& cohort,
                            bool& full) {
  const int m = static_cast<int>(dose.size());
  const int n_cohorts = static_cast<int>(des.cohort_size.size());
  const int w = des.max_n + 1;
  std::vector<char> excl(m, 0);
  int i = start;
  full = false;
  while (true) {
    const int c = dose[i];
    for (int p = 0; p < des.cohort_size[cohort]; ++p) {
      bool dlt, resp;
      etc_patient(des, c, rng, dlt, resp);
      n[c] += 1;
      if (dlt) y[c] += 1;
      if (resp) x[c] += 1;
      if (dlt && resp) z[c] += 1;
    }
    ++cohort;
    if (des.tox_ok[n[c] * w + y[c]] == 0) {
      for (int k = i; k < m; ++k) {
        excl[k] = 1;
        excluded[dose[k]] = 1;
      }
      if (i == 0) {
        full = cohort >= n_cohorts;
        return -1;
      }
    }
    const int code = des.dec[(n[c] - 1) * w * w + y[c] * w + x[c]];
    int nxt = i;
    if (code == ET_ESCALATE) {
      if (i + 1 < m && !excl[i + 1]) nxt = i + 1;
    } else if (code == ET_DEESCALATE) {
      if (i > 0) nxt = i - 1;
    } else if (code == ET_ADMISSIBLE_SET) {
      nxt = etc_rule4(i, dose, excl, &n[0], &x[0], des.rule4, rng);
    }
    i = nxt;
    if (cohort >= n_cohorts) {
      full = true;
      return dose[i];
    }
    if (n[dose[i]] >= des.n_stop) return dose[i];
  }
}

// Simulate one BOIN-ETC trial.
template <class Rng>
inline void etc_simulate_one(const EtcDesign& des, Rng& rng,
                             std::vector<int>& n, std::vector<int>& y,
                             std::vector<int>& x, std::vector<int>& z,
                             std::vector<char>& excluded, int& cohorts_used,
                             int& n_subtrials, int& stop_code) {
  const int J = des.n_rows;
  const int K = des.n_cols;
  std::fill(n.begin(), n.end(), 0);
  std::fill(y.begin(), y.end(), 0);
  std::fill(x.begin(), x.end(), 0);
  std::fill(z.begin(), z.end(), 0);
  std::fill(excluded.begin(), excluded.end(), 0);
  int cohort = 0;
  bool full = false;
  n_subtrials = 1;

  int cand = etc_run_subtrial(des, etc_first_subtrial(J, K), 0, rng, n, y, x,
                              z, excluded, cohort, full);
  cohorts_used = cohort;
  if (cand < 0) {
    stop_code = ETC_STOP_LOWEST_EXCLUDED;
    return;
  }
  if (full) {
    stop_code = ETC_STOP_MAX_COHORTS;
    return;
  }

  if (des.subtrials == ETC_ETC3) {
    for (int j = 1; j < std::min(J, K); ++j) {
      etc_run_subtrial(des, etc_lshape_subtrial(J, K, j), 0, rng, n, y, x, z,
                       excluded, cohort, full);
      ++n_subtrials;
      cohorts_used = cohort;
      if (full) {
        stop_code = ETC_STOP_MAX_COHORTS;
        return;
      }
    }
    stop_code = ETC_STOP_COMPLETED;
    return;
  }

  const int a0 = cand % J;
  const int b0 = cand / J;
  int row = (b0 == 0) ? std::max(a0 - 1, 0) : J - 2;
  int prev = cand;
  for (int r = row; r >= 0; --r) {
    int start = 0;
    if (des.subtrials == ETC_ETC1 && prev >= 0 && prev / J > 0) {
      int col = prev / J;
      if (des.start_rule == ETC_START_DIAGONAL) col -= 1;
      if (col < 1) col = 1;
      start = col - 1;
    }
    prev = etc_run_subtrial(des, etc_row_subtrial(J, K, r), start, rng, n, y,
                            x, z, excluded, cohort, full);
    ++n_subtrials;
    cohorts_used = cohort;
    if (full) {
      stop_code = ETC_STOP_MAX_COHORTS;
      return;
    }
  }
  stop_code = ETC_STOP_COMPLETED;
}

// Unweighted bivariate isotonic estimates of the DLT rates (y + 0.05) /
// (n + 0.1), non-decreasing in a and in b. Returns false when
// ETC_ISO_MAX_CYCLES cycles are not enough.
inline bool etc_biviso(const int* n, const int* y, int J, int K,
                       double* iso) {
  const int m = J * K;
  std::vector<double> wt(m, 1.0), p(m, 0.0), q(m, 0.0), v(m), r(m);
  for (int c = 0; c < m; ++c) iso[c] = (y[c] + 0.05) / (n[c] + 0.1);
  for (int cycle = 0; cycle < ETC_ISO_MAX_CYCLES; ++cycle) {
    for (int c = 0; c < m; ++c) v[c] = iso[c] + p[c];
    for (int a = 0; a < J; ++a) cb12_pava(&v[a], &wt[a], K, J, &r[a]);
    for (int c = 0; c < m; ++c) {
      p[c] = v[c] - r[c];
      v[c] = r[c] + q[c];
    }
    for (int b = 0; b < K; ++b) {
      cb12_pava(&v[J * b], &wt[J * b], J, 1, &r[J * b]);
    }
    double change = 0.0;
    for (int c = 0; c < m; ++c) {
      q[c] = v[c] - r[c];
      change = std::max(change, std::fabs(r[c] - iso[c]));
      iso[c] = r[c];
    }
    if (change < ETC_ISO_TOL) return true;
  }
  return false;
}

// Select the OBD of one trial from its final data and the excluded flags.
// Returns its index (0-based) or -1, sets reason (EtObdReason), the MTD
// column of every row (mtd_col, size J), the tied candidates (tied, size
// J * K, 1 when within ETC_UTIL_TOL of the best) and the utility estimates
// (0 to 100) of the tried combinations (util; the others are left
// untouched).
inline int etc_select_one(const int* n, const int* y, const int* x,
                          const int* z, const char* excluded,
                          const double* iso, bool early_stop, int J, int K,
                          double target, double psi00, double both_weight,
                          double alpha, double beta, int& reason,
                          int* mtd_col, char* tied, double* util) {
  const int m = J * K;
  for (int c = 0; c < m; ++c) tied[c] = 0;
  for (int a = 0; a < J; ++a) mtd_col[a] = -1;
  if (early_stop) {
    reason = ET_OBD_EARLY_STOP;
    return -1;
  }
  for (int a = 0; a < J; ++a) {
    double best = -1.0;
    for (int b = 0; b < K; ++b) {
      const double d = std::fabs(iso[a + J * b] - target);
      if (best < 0.0 || d < best) best = d;
    }
    for (int b = 0; b < K; ++b) {
      if (std::fabs(iso[a + J * b] - target) <= best + ETC_MTD_TOL) {
        mtd_col[a] = b;
      }
    }
  }
  bool any = false;
  double best = 0.0;
  std::vector<char> cand(m, 0);
  for (int c = 0; c < m; ++c) {
    if (n[c] == 0) continue;
    const double q = (100.0 * x[c] + psi00 * (n[c] - y[c] - x[c]) +
                      both_weight * z[c]) / 100.0;
    util[c] = 100.0 * (q + alpha) / (n[c] + alpha + beta);
    if (excluded[c] || c / J > mtd_col[c % J]) continue;
    cand[c] = 1;
    if (!any || util[c] > best) best = util[c];
    any = true;
  }
  if (!any) {
    reason = ET_OBD_NO_ADMISSIBLE;
    return -1;
  }
  int obd = -1;
  int best_s = 0;
  for (int c = 0; c < m; ++c) {
    if (!cand[c] || util[c] < best - ETC_UTIL_TOL) continue;
    tied[c] = 1;
    const int s = c % J + c / J;
    if (obd < 0 || s < best_s) {
      obd = c;
      best_s = s;
    }
  }
  reason = ET_OBD_SELECTED;
  return obd;
}

}  // namespace simfastboin

#endif
