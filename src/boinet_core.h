#ifndef SIMFASTBOIN_BOINET_CORE_H
#define SIMFASTBOIN_BOINET_CORE_H

// Core algorithms of the BOIN-ET design (Takeda, Taguri and Morita, 2018),
// written in plain C++ like boin_core.h so that they can be verified outside
// of R as well.
//
// The rules follow the boinet package (version 1.6.0), with two changes that
// keep to Takeda et al. (2018), who do not allow dose skipping during
// escalation: a dose excluded by its DLT data also excludes every higher dose,
// as in Takeda et al. (2025) and the bfboinet package, and an escalation from
// the highest dose ends the trial when that dose is excluded. With compat set,
// the rules are exactly those of boinet 1.6.0. Every comparison of an
// observed rate with a boundary, and every posterior probability that decides
// whether a dose is admissible, is evaluated in R and passed in as a table, so
// the tables carry exactly the floating point results of the reference.
//
// Table layouts, with M the largest number of patients at one dose:
//   dec[(n - 1) * (M + 1) * (M + 1) + y * (M + 1) + x], n = 1..M
//       decision for y DLTs and x responses among n patients (EtDecision);
//   tox_ok[n * (M + 1) + y], n = 0..M
//       1 when the DLT data do not exclude the dose;
//   resp_ok[n * (M + 1) + x], n = 0..M
//       1 when the response data do not exclude the dose.

#include <vector>
#include <cmath>
#include <algorithm>
#include "backfill_core.h"

namespace simfastboin {

// Decision for the next cohort from the data at the current dose.
enum EtDecision {
  ET_ESCALATE = 0,
  ET_STAY = 1,
  ET_DEESCALATE = 2,
  ET_ADMISSIBLE_SET = 3   // rule 4: choose from the doses around the current one
};

// Stop codes returned by simulate_et_one().
enum EtStop {
  ET_STOP_NO_ADMISSIBLE = 0,  // every dose is excluded
  ET_STOP_NO_LOWER = 1,       // de-escalation finds no admissible dose
  ET_STOP_N_EARLYSTOP = 2,    // a dose has reached the stopping sample size
  ET_STOP_MAX_COHORTS = 3,    // all cohorts were enrolled
  ET_STOP_TOP_EXCLUDED = 4    // escalation from the highest dose, which is
                              // excluded (not an early stop; compat == 0 only)
};

// How the response of a patient is linked to the DLT.
enum EtCorType {
  ET_COR_COPULA = 0,   // Gaussian copula with correlation resp_cor
  ET_COR_BINARY = 1    // conditional response probabilities q0 and q1
};

// Methods to select the OBD among the admissible doses.
enum EtObdMethod {
  ET_OBD_MAX_RESP = 0,
  ET_OBD_UTILITY_WEIGHTED = 1,
  ET_OBD_UTILITY_TRUNCATED_LINEAR = 2,
  ET_OBD_UTILITY_SCORING = 3
};

// Reason codes returned by select_obd_et_one().
enum EtObdReason {
  ET_OBD_SELECTED = 0,
  ET_OBD_EARLY_STOP = 1,
  ET_OBD_NO_ADMISSIBLE = 2
};

// The settings of a BOIN-ET trial that do not change from trial to trial.
struct EtDesign {
  std::vector<double> p_true;      // DLT probability per dose
  std::vector<double> p_resp;      // response probability per dose
  std::vector<double> q0;          // Pr(response | no DLT), binary link only
  std::vector<double> q1;          // Pr(response | DLT), binary link only
  int cor_type;                    // EtCorType
  double resp_cor;                 // copula correlation
  std::vector<int> cohort_size;    // one entry per cohort
  int start_dose;                  // 1-based
  int n_stop;                      // stopping sample size at any dose
  int max_n;                       // M of the tables
  std::vector<int> dec;
  std::vector<int> tox_ok;
  std::vector<int> resp_ok;
  int compat;                      // 1 for the rules of boinet 1.6.0
};

// Whether a dose with n patients, y DLTs and x responses is admissible.
inline bool et_admissible(int n, int y, int x, const std::vector<int>& tox_ok,
                          const std::vector<int>& resp_ok, int max_n) {
  const int w = max_n + 1;
  return tox_ok[n * w + y] == 1 && resp_ok[n * w + x] == 1;
}

// Admissibility of every dose from its patients nn, DLTs yy and responses xx.
// Unless compat is set, a dose excluded by its DLT data also excludes every
// higher dose. Returns whether any dose is admissible.
inline bool et_admissible_set(const std::vector<int>& nn,
                              const std::vector<int>& yy,
                              const std::vector<int>& xx,
                              const std::vector<int>& tox_ok,
                              const std::vector<int>& resp_ok, int max_n,
                              bool compat, std::vector<char>& adm) {
  const int n_doses = static_cast<int>(nn.size());
  const int w = max_n + 1;
  bool tox_out = false;
  bool any_adm = false;
  for (int j = 0; j < n_doses; ++j) {
    const bool t_ok = tox_ok[nn[j] * w + yy[j]] == 1;
    if (!t_ok && !compat) tox_out = true;
    const bool ok = t_ok && !tox_out && resp_ok[nn[j] * w + xx[j]] == 1;
    adm[j] = ok ? 1 : 0;
    if (ok) any_adm = true;
  }
  return any_adm;
}

// Move from the current dose cur (1-based) to the dose nxt chosen by the
// decision rule, given the admissibility adm of every dose. Returns false,
// with stop_code set, when the trial stops. When nxt is not admissible, an
// escalation moves to the lowest admissible dose above it (the next cohort
// stays when there is none), and a stay or a de-escalation to the highest
// admissible dose below it (the trial stops when there is none). An
// escalation from the highest dose stays there; unless compat is set, the
// trial stops instead when that dose is not admissible.
inline bool et_move(int& cur, int nxt, const std::vector<char>& adm,
                    bool compat, int& stop_code) {
  const int n_doses = static_cast<int>(adm.size());
  if (nxt == 0) {
    if (adm[0]) {
      cur = 1;
      return true;
    }
    stop_code = ET_STOP_NO_LOWER;
    return false;
  }
  if (nxt == n_doses + 1) {
    if (compat || adm[n_doses - 1]) {
      cur = n_doses;
      return true;
    }
    stop_code = ET_STOP_TOP_EXCLUDED;
    return false;
  }
  if (adm[nxt - 1]) {
    cur = nxt;
    return true;
  }
  if (cur < nxt) {
    for (int j = nxt; j <= n_doses; ++j) {
      if (adm[j - 1]) {
        cur = j;
        break;
      }
    }
    return true;
  }
  for (int j = nxt; j >= 1; --j) {
    if (adm[j - 1]) {
      cur = j;
      return true;
    }
  }
  stop_code = ET_STOP_NO_LOWER;
  return false;
}

// Next dose (1-based, from 0 to n_doses + 1) chosen by the decision rule.
//
// Under rule 4 the dose with the highest observed response rate among the
// current dose and its neighbors is chosen; a dose not yet tried counts as a
// rate of zero, as in the reference. Ties are broken with one uniform variate.
template <class Rng>
inline int et_next_dose(int cur, const std::vector<int>& nn,
                        const std::vector<int>& yy, const std::vector<int>& xx,
                        const EtDesign& des, Rng& rng) {
  const int n_doses = static_cast<int>(nn.size());
  const int w = des.max_n + 1;
  const int c = cur - 1;
  const int code = des.dec[(nn[c] - 1) * w * w + yy[c] * w + xx[c]];
  if (code == ET_ESCALATE) return cur + 1;
  if (code == ET_STAY) return cur;
  if (code == ET_DEESCALATE) return cur - 1;

  std::vector<int> three;
  if (cur == n_doses) {
    if (n_doses > 1) three.push_back(cur - 1);
    three.push_back(cur);
  } else if (nn[cur] == 0) {
    return cur + 1;
  } else if (cur == 1) {
    three.push_back(cur);
    three.push_back(cur + 1);
  } else {
    three.push_back(cur - 1);
    three.push_back(cur);
    three.push_back(cur + 1);
  }

  std::vector<double> pe(three.size());
  double best = -1.0;
  for (std::size_t k = 0; k < three.size(); ++k) {
    const int d = three[k] - 1;
    pe[k] = (nn[d] > 0) ? static_cast<double>(xx[d]) / static_cast<double>(nn[d])
                        : 0.0;
    if (pe[k] > best) best = pe[k];
  }
  std::vector<int> ties;
  for (std::size_t k = 0; k < three.size(); ++k) {
    if (pe[k] == best) ties.push_back(three[k]);
  }
  if (ties.size() == 1) return ties[0];
  const int k = static_cast<int>(ties.size());
  int pick = static_cast<int>(rng() * k);
  if (pick >= k) pick = k - 1;
  return ties[pick];
}

// Simulate a single BOIN-ET trial.
//
// rng() must return one uniform variate on (0, 1) per call. Each patient draws
// two variates in enrollment order, the first for the DLT and the second for
// the response, and a tie under rule 4 draws one more.
template <class Rng>
inline void simulate_et_one(const EtDesign& des, Rng& rng,
                            std::vector<int>& nn, std::vector<int>& yy,
                            std::vector<int>& xx, int& cohorts_used,
                            int& stop_code) {
  const int n_doses = static_cast<int>(des.p_true.size());
  std::fill(nn.begin(), nn.end(), 0);
  std::fill(yy.begin(), yy.end(), 0);
  std::fill(xx.begin(), xx.end(), 0);
  std::vector<char> adm(n_doses);

  int cur = des.start_dose;
  cohorts_used = 0;
  stop_code = ET_STOP_MAX_COHORTS;

  for (std::size_t i = 0; i < des.cohort_size.size(); ++i) {
    ++cohorts_used;
    const int d = cur - 1;
    for (int k = 0; k < des.cohort_size[i]; ++k) {
      const double u = rng();
      const double v = rng();
      const bool dlt = u < des.p_true[d];
      bool resp;
      if (des.cor_type == ET_COR_BINARY) {
        resp = v < (dlt ? des.q1[d] : des.q0[d]);
      } else {
        resp = response_variate(u, v, des.resp_cor) < des.p_resp[d];
      }
      nn[d] += 1;
      if (dlt) yy[d] += 1;
      if (resp) xx[d] += 1;
    }

    const int nxt = et_next_dose(cur, nn, yy, xx, des, rng);

    const bool any_adm = et_admissible_set(nn, yy, xx, des.tox_ok, des.resp_ok,
                                           des.max_n, des.compat != 0, adm);
    if (!any_adm) {
      stop_code = ET_STOP_NO_ADMISSIBLE;
      return;
    }
    bool reached = false;
    for (int j = 0; j < n_doses; ++j) {
      if (nn[j] >= des.n_stop) reached = true;
    }
    if (reached) {
      stop_code = ET_STOP_N_EARLYSTOP;
      return;
    }

    if (!et_move(cur, nxt, adm, des.compat != 0, stop_code)) return;
  }
}

// Isotonic regression with equal weights (pool adjacent violators).
inline void pava_equal(std::vector<double>& x) {
  const int n = static_cast<int>(x.size());
  std::vector<double> sum;
  std::vector<int> count;
  sum.reserve(n);
  count.reserve(n);
  for (int i = 0; i < n; ++i) {
    sum.push_back(x[i]);
    count.push_back(1);
    while (sum.size() > 1) {
      const std::size_t b = sum.size() - 1;
      if (sum[b - 1] / count[b - 1] > sum[b] / count[b]) {
        sum[b - 1] += sum[b];
        count[b - 1] += count[b];
        sum.pop_back();
        count.pop_back();
      } else {
        break;
      }
    }
  }
  int pos = 0;
  for (std::size_t b = 0; b < sum.size(); ++b) {
    const double value = sum[b] / count[b];
    for (int k = 0; k < count[b]; ++k) x[pos++] = value;
  }
}

// The parameters of the OBD selection.
struct EtObdRule {
  int method;          // EtObdMethod
  double target;       // target DLT probability
  double tox_upper;    // DLT rate penalised by the weighted utility
  double w1;
  double w2;
  double t_low;        // truncated linear utility
  double t_upp;
  double e_low;
  double e_upp;
  double psi00;        // scoring utility
  double psi11;
};

// Index, among the evaluated doses, of the OBD chosen by rule from the
// isotonic DLT estimates probt, the response estimates probe and the
// admissibility flags cand. At least one dose must be admissible.
inline int et_obd_choose(const std::vector<double>& probt,
                         const std::vector<double>& probe,
                         const std::vector<char>& cand, const EtObdRule& rule) {
  const int k_doses = static_cast<int>(probt.size());

  if (rule.method == ET_OBD_MAX_RESP) {
    // The MTD is the last evaluated dose whose estimate is as close to the
    // target as the closest admissible one; the OBD is the first admissible
    // dose up to the MTD with the highest response estimate.
    double mdif = 0.0;
    bool first = true;
    for (int k = 0; k < k_doses; ++k) {
      if (!cand[k]) continue;
      const double dist = std::fabs(probt[k] - rule.target);
      if (first || dist < mdif) {
        mdif = dist;
        first = false;
      }
    }
    int mtd = 0;
    for (int k = 0; k < k_doses; ++k) {
      if (std::fabs(probt[k] - rule.target) == mdif) mtd = k;
    }
    double meff = 0.0;
    first = true;
    for (int k = 0; k <= mtd; ++k) {
      if (!cand[k]) continue;
      if (first || probe[k] > meff) {
        meff = probe[k];
        first = false;
      }
    }
    for (int k = 0; k < k_doses; ++k) {
      if (cand[k] && probe[k] == meff) return k;
    }
    return -1;
  }

  std::vector<double> fu(k_doses);
  for (int k = 0; k < k_doses; ++k) {
    const double pt = probt[k];
    const double pe = probe[k];
    if (rule.method == ET_OBD_UTILITY_WEIGHTED) {
      const double over = (pt > rule.tox_upper) ? 1.0 : 0.0;
      fu[k] = (pe - rule.w1 * pt) - rule.w2 * pt * over;
    } else if (rule.method == ET_OBD_UTILITY_TRUNCATED_LINEAR) {
      double fpe, fpt;
      if (pe <= rule.e_low) {
        fpe = 0.0;
      } else if (pe >= rule.e_upp) {
        fpe = 1.0;
      } else {
        fpe = (pe - rule.e_low) / (rule.e_upp - rule.e_low);
      }
      if (pt <= rule.t_low) {
        fpt = 1.0;
      } else if (pt >= rule.t_upp) {
        fpt = 0.0;
      } else {
        fpt = 1.0 - (pt - rule.t_low) / (rule.t_upp - rule.t_low);
      }
      fu[k] = fpe * fpt;
    } else {
      fu[k] = 0.0 * (1.0 - pe) * pt + rule.psi00 * (1.0 - pe) * (1.0 - pt) +
        rule.psi11 * pe * pt + 100.0 * pe * (1.0 - pt);
    }
  }
  double best = 0.0;
  bool first = true;
  for (int k = 0; k < k_doses; ++k) {
    if (!cand[k]) continue;
    if (first || fu[k] > best) {
      best = fu[k];
      first = false;
    }
  }
  for (int k = 0; k < k_doses; ++k) {
    if (cand[k] && fu[k] == best) return k;
  }
  return -1;
}

// Select the OBD from the final data of a single trial.
//
// probe holds the response estimate of every dose; only the doses with at
// least one patient are read. Unless compat is set, a dose excluded by its DLT
// data also excludes every higher dose. Returns the OBD (1-based) or 0 when no
// dose is selected, and sets reason to an EtObdReason.
inline int select_obd_et_one(const int* n, const int* y, const int* x,
                             const double* probe_all, int n_doses,
                             bool early_stop, const std::vector<int>& tox_ok,
                             const std::vector<int>& resp_ok, int max_n,
                             bool compat, const EtObdRule& rule, int& reason) {
  if (early_stop) {
    reason = ET_OBD_EARLY_STOP;
    return 0;
  }
  std::vector<int> eva;
  for (int j = 0; j < n_doses; ++j) {
    if (n[j] > 0) eva.push_back(j);
  }
  const int k_doses = static_cast<int>(eva.size());
  std::vector<int> nv(n, n + n_doses), yv(y, y + n_doses), xv(x, x + n_doses);
  std::vector<char> adm_all(n_doses);
  et_admissible_set(nv, yv, xv, tox_ok, resp_ok, max_n, compat, adm_all);
  std::vector<char> cand(k_doses);
  bool any_cand = false;
  for (int k = 0; k < k_doses; ++k) {
    cand[k] = adm_all[eva[k]];
    if (cand[k]) any_cand = true;
  }
  if (!any_cand) {
    reason = ET_OBD_NO_ADMISSIBLE;
    return 0;
  }
  if (k_doses == 1) {
    reason = ET_OBD_SELECTED;
    return eva[0] + 1;
  }

  std::vector<double> probt(k_doses), probe(k_doses);
  for (int k = 0; k < k_doses; ++k) {
    const int j = eva[k];
    probt[k] = static_cast<double>(y[j]) / static_cast<double>(n[j]);
    probe[k] = probe_all[j];
  }
  pava_equal(probt);

  const int k = et_obd_choose(probt, probe, cand, rule);
  if (k < 0) {
    reason = ET_OBD_NO_ADMISSIBLE;
    return 0;
  }
  reason = ET_OBD_SELECTED;
  return eva[k] + 1;
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_BOINET_CORE_H
