#ifndef SIMFASTBOIN_BOIN12_CORE_H
#define SIMFASTBOIN_BOIN12_CORE_H

// Core algorithms of the BOIN12 design (Lin, Zhou, Yan, Li and Yuan, 2020),
// written in plain C++ like boinet_core.h so that they can be verified
// outside of R as well.
//
// The utility of a patient is 100 for a response without DLT, psi00 for
// neither, psi11 for both and 0 for a DLT without response. The utility of a
// dose, standardized to (0, 1), has the Beta(alpha + x, beta + n - x)
// posterior of the quasi-binomial model, with x = sum of the utilities / 100,
// and the desirability of a dose is the posterior probability that its
// utility exceeds the benchmark u_benchmark. The posterior probabilities come
// from the function object pb(q, a, b) = Pr(Beta(a, b) > q), which is
// R::pbeta in the package.
//
// Rules for the next cohort at the current dose d, with phat = y / n:
//   - a dose is admissible when its DLT and response data both pass (tables
//     tox_ok and resp_ok); a dose excluded by its DLT data also excludes every
//     higher dose; with no admissible dose the trial stops;
//   - exploration: n >= n_explore, phat < lambda_d and the dose above is
//     untried and admissible (with explore_admissible_only, the current dose
//     must be admissible as well): escalate;
//   - phat >= lambda_d: the highest admissible dose below d, or d itself;
//   - lambda_e < phat < lambda_d and n >= n_star: the more desirable of d - 1
//     and d; otherwise the most desirable of d - 1, d and d + 1; only
//     admissible doses count and ties go to the higher dose, or with
//     tie_random to one of the tied doses drawn at random;
//   - when none of these doses is admissible, the admissible dose nearest to
//     d, ties to the lower dose.
// With stop_no_candidate, the rule of the simulations of Lin et al. (2020),
// the trial stops instead in the last case and when phat >= lambda_d with no
// admissible dose below d > 1; at the lowest dose it still stays when that
// dose is admissible. explore_admissible_only is the rule of exploration of
// these simulations.
// The trial stops early when n_earlystop > 0, the current dose has
// n_earlystop patients and the next dose is the current one (or, with
// early_stop_simple, in any case).
//
// Table layouts, with M the largest number of patients at one dose:
//   tox_ok[n * (M + 1) + y], n = 0..M: 1 when y DLTs do not exclude the dose;
//   resp_ok[n * (M + 1) + x], n = 0..M: 1 when x responses do not.

#include <vector>
#include <cmath>
#include <algorithm>
#include "boinet_core.h"

namespace simfastboin {

// Stop codes returned by b12_simulate_one().
enum B12Stop {
  B12_STOP_NO_ADMISSIBLE = 0,  // every dose is excluded
  B12_STOP_N_EARLYSTOP = 1,    // the current dose has n_earlystop patients
  B12_STOP_MAX_COHORTS = 2,    // all cohorts were enrolled
  B12_STOP_NO_CANDIDATE = 3    // no admissible dose to go to, with
                               // stop_no_candidate
};

// Value of b12_next_dose() and tb12_next_dose() when stop_no_candidate stops
// the trial.
const int B12_NEXT_STOP = -1;

// Estimates of the utility used to select the OBD.
enum B12Est {
  B12_EST_POSTERIOR_MEAN = 0,  // (x + alpha) / (n + alpha + beta)
  B12_EST_OBSERVED = 1,        // x / n
  B12_EST_TAIL = 2             // the desirability Pr(u > u_benchmark | data)
};

// The settings of a BOIN12 trial that do not change from trial to trial.
struct B12Design {
  std::vector<double> p_true;      // DLT probability per dose
  std::vector<double> p_resp;      // response probability per dose
  std::vector<double> q0;          // Pr(response | no DLT), binary link only
  std::vector<double> q1;          // Pr(response | DLT), binary link only
  int cor_type;                    // EtCorType
  double resp_cor;                 // copula correlation
  std::vector<int> cohort_size;    // one entry per cohort
  int start_dose;                  // 1-based
  int n_earlystop;                 // 0 for no early stopping
  int early_stop_simple;
  int n_star;
  int n_explore;
  int tie_random;                  // 1: break ties of desirability at random
  int stop_no_candidate;           // 1: stop when no candidate is admissible
  int explore_admissible_only;     // 1: explore only from an admissible dose
  double lambda_e;
  double lambda_d;
  double psi00;
  double psi11;
  double alpha;
  double beta;
  double u_benchmark;              // on the (0, 1) scale
  int max_n;                       // M of the tables
  std::vector<int> tox_ok;
  std::vector<int> resp_ok;
};

// Quasi-binomial number of events of n patients with y DLTs, x responses and
// z patients with both.
inline double b12_quasi_events(double n, double y, double x, double z,
                               double psi00, double psi11) {
  return (100.0 * (x - z) + psi00 * (n - y - x + z) + psi11 * z) / 100.0;
}

// Desirability of a dose with n patients and quasi-events q.
template <class Pb>
inline double b12_desirability(double n, double q, const B12Design& des,
                               const Pb& pb) {
  return pb(des.u_benchmark, des.alpha + q, des.beta + n - q);
}

// One of the n_tied doses tied[0] < tied[1] < ... with the highest
// desirability (0-based): the highest one, or with tie_random the dose
// tied[floor(u * n_tied)] for one uniform variate u = rng(), which is drawn
// only when n_tied > 1.
template <class Rng>
inline int b12_break_tie(const int* tied, int n_tied, const B12Design& des,
                         Rng& rng) {
  if (n_tied == 1 || !des.tie_random) return tied[n_tied - 1];
  int i = static_cast<int>(rng() * n_tied);
  if (i >= n_tied) i = n_tied - 1;
  return tied[i];
}

// Next dose (1-based) from the current dose cur, 0 when no dose is admissible
// and B12_NEXT_STOP when stop_no_candidate stops the trial. adm receives the
// admissibility of every dose; rng() gives the variate of a random tie break.
template <class Pb, class Rng>
inline int b12_next_dose(int cur, const std::vector<int>& nn,
                         const std::vector<int>& yy,
                         const std::vector<int>& xx,
                         const std::vector<int>& zz, const B12Design& des,
                         const Pb& pb, Rng& rng, std::vector<char>& adm) {
  const int n_doses = static_cast<int>(nn.size());
  if (!et_admissible_set(nn, yy, xx, des.tox_ok, des.resp_ok, des.max_n,
                         false, adm)) {
    return 0;
  }
  const int d = cur - 1;
  const double phat = nn[d] > 0 ? static_cast<double>(yy[d]) / nn[d] : 0.0;

  if (nn[d] >= des.n_explore && phat < des.lambda_d && d + 1 < n_doses &&
      nn[d + 1] == 0 && adm[d + 1] &&
      (!des.explore_admissible_only || adm[d])) {
    return cur + 1;
  }

  if (phat >= des.lambda_d) {
    for (int j = d - 1; j >= 0; --j) {
      if (adm[j]) return j + 1;
    }
    if (des.stop_no_candidate && d > 0) return B12_NEXT_STOP;
    if (adm[d]) return cur;
  } else {
    const bool pair = phat > des.lambda_e && nn[d] >= des.n_star;
    const int top = pair ? d : d + 1;
    int tied[3];
    int n_tied = 0;
    double best_value = -1.0;
    for (int j = d - 1; j <= top; ++j) {
      if (j < 0 || j >= n_doses || !adm[j]) continue;
      const double q = b12_quasi_events(nn[j], yy[j], xx[j], zz[j], des.psi00,
                                        des.psi11);
      const double value = b12_desirability(nn[j], q, des, pb);
      if (value > best_value) {
        best_value = value;
        tied[0] = j;
        n_tied = 1;
      } else if (value == best_value) {
        tied[n_tied++] = j;
      }
    }
    if (n_tied > 0) return b12_break_tie(tied, n_tied, des, rng) + 1;
  }

  if (des.stop_no_candidate) return B12_NEXT_STOP;
  for (int k = 1; k < n_doses; ++k) {
    if (d - k >= 0 && adm[d - k]) return d - k + 1;
    if (d + k < n_doses && adm[d + k]) return d + k + 1;
  }
  return 0;
}

// Draw the DLT and the response of one patient at dose d (0-based) from two
// uniform variates u and v.
inline void b12_patient(const B12Design& des, int d, double u, double v,
                        bool& dlt, bool& resp) {
  dlt = u < des.p_true[d];
  if (des.cor_type == ET_COR_BINARY) {
    resp = v < (dlt ? des.q1[d] : des.q0[d]);
  } else {
    resp = response_variate(u, v, des.resp_cor) < des.p_resp[d];
  }
}

// Simulate a single BOIN12 trial.
//
// rng() must return one uniform variate on (0, 1) per call. Each patient draws
// two variates in enrollment order, the first for the DLT and the second for
// the response. With tie_random, a decision between tied doses then draws one
// more variate.
template <class Rng, class Pb>
inline void b12_simulate_one(const B12Design& des, Rng& rng, const Pb& pb,
                             std::vector<int>& nn, std::vector<int>& yy,
                             std::vector<int>& xx, std::vector<int>& zz,
                             int& cohorts_used, int& stop_code) {
  const int n_doses = static_cast<int>(des.p_true.size());
  std::fill(nn.begin(), nn.end(), 0);
  std::fill(yy.begin(), yy.end(), 0);
  std::fill(xx.begin(), xx.end(), 0);
  std::fill(zz.begin(), zz.end(), 0);
  std::vector<char> adm(n_doses);

  int cur = des.start_dose;
  cohorts_used = 0;
  stop_code = B12_STOP_MAX_COHORTS;

  for (std::size_t i = 0; i < des.cohort_size.size(); ++i) {
    ++cohorts_used;
    const int d = cur - 1;
    for (int k = 0; k < des.cohort_size[i]; ++k) {
      const double u = rng();
      const double v = rng();
      bool dlt, resp;
      b12_patient(des, d, u, v, dlt, resp);
      nn[d] += 1;
      if (dlt) yy[d] += 1;
      if (resp) xx[d] += 1;
      if (dlt && resp) zz[d] += 1;
    }

    const int nxt = b12_next_dose(cur, nn, yy, xx, zz, des, pb, rng, adm);
    if (nxt == 0) {
      stop_code = B12_STOP_NO_ADMISSIBLE;
      return;
    }
    if (nxt == B12_NEXT_STOP) {
      stop_code = B12_STOP_NO_CANDIDATE;
      return;
    }
    if (des.n_earlystop > 0 && nn[d] >= des.n_earlystop &&
        (des.early_stop_simple || nxt == cur)) {
      stop_code = B12_STOP_N_EARLYSTOP;
      return;
    }
    cur = nxt;
  }
}

// Select the OBD of one trial from its final data (n, y, x, z per dose).
// target is the DLT limit of the MTD. Returns the OBD (1-based) or 0, sets
// reason (EtObdReason) and writes the utility estimates of the doses tried to
// util (on the scale 0 to 100; the other doses are left untouched).
//
// The DLT rates of the doses tried are made monotone by isotonic regression
// with equal weights, and the MTD is the dose tried whose estimate is closest
// to target, ties (within 1e-10) to the higher dose. The OBD is the
// admissible dose tried, up to the MTD, with the highest utility estimate,
// ties (within 1e-12) to the lower dose.
template <class Pb>
inline int b12_select_obd_one(const int* n, const int* y, const int* x,
                              const int* z, int n_doses, bool early_stop,
                              const B12Design& des, double target, int est,
                              const Pb& pb, int& reason, double* util) {
  if (early_stop) {
    reason = ET_OBD_EARLY_STOP;
    return 0;
  }
  std::vector<int> nn(n, n + n_doses), yy(y, y + n_doses), xx(x, x + n_doses);
  std::vector<char> adm(n_doses);
  if (!et_admissible_set(nn, yy, xx, des.tox_ok, des.resp_ok, des.max_n,
                         false, adm)) {
    reason = ET_OBD_NO_ADMISSIBLE;
    return 0;
  }

  std::vector<int> tried;
  std::vector<double> rate;
  for (int j = 0; j < n_doses; ++j) {
    if (n[j] > 0) {
      tried.push_back(j);
      rate.push_back(static_cast<double>(y[j]) / n[j]);
    }
  }
  pava_equal(rate);
  int mtd = tried[0];
  double best_gap = std::fabs(rate[0] - target);
  for (std::size_t k = 1; k < tried.size(); ++k) {
    const double gap = std::fabs(rate[k] - target);
    if (gap <= best_gap + 1e-10) {
      if (gap < best_gap) best_gap = gap;
      mtd = tried[k];
    }
  }

  int obd = -1;
  double best = 0.0;
  for (std::size_t k = 0; k < tried.size(); ++k) {
    const int j = tried[k];
    const double q = b12_quasi_events(n[j], y[j], x[j], z[j], des.psi00,
                                      des.psi11);
    double value;
    if (est == B12_EST_OBSERVED) {
      value = q / n[j];
    } else if (est == B12_EST_TAIL) {
      value = b12_desirability(n[j], q, des, pb);
    } else {
      value = (q + des.alpha) / (n[j] + des.alpha + des.beta);
    }
    util[j] = 100.0 * value;
    if (j > mtd || !adm[j]) continue;
    if (obd < 0 || value > best + 1e-12) {
      obd = j;
      best = value;
    }
  }
  if (obd < 0) {
    reason = ET_OBD_NO_ADMISSIBLE;
    return 0;
  }
  reason = ET_OBD_SELECTED;
  return obd + 1;
}

}  // namespace simfastboin

#endif
