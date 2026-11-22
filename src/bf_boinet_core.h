#ifndef SIMFASTBOIN_BF_BOINET_CORE_H
#define SIMFASTBOIN_BF_BOINET_CORE_H

// Core algorithm of the BF-BOIN-ET design (Takeda et al., 2025), written in
// plain C++ like boinet_core.h so that it can be verified outside of R as
// well. The rules follow get.oc.backboinet() of the bfboinet package (version
// 1.1.0):
//   - a cohort is treated at the current dose; the next cohort starts with
//     the first patient who arrives after the DLT outcomes of the cohort are
//     known (all patients followed for the DLT window or all DLTs observed),
//     or, with a suspension rule, after enough of them are known; after the
//     last cohort, backfilling continues for a further period;
//   - the patients arriving in between are backfilled, one by one, to the
//     dose with the highest utility among the backfill doses that have fewer
//     than n_stop patients, the DLT rates being made monotone first;
//   - the decision uses the time-to-event DLT rate at the current dose with
//     the boundaries of BOIN, a stay being kept after 1 DLT in 3 patients; the
//     backfill doses are those from the lowest dose with a response rate of
//     at least eta1 up to the current dose (escalation), the dose below it
//     (stay) or two below it (de-escalation), among the admissible doses;
//   - once a dose has n_elim patients it is excluded, with every higher dose,
//     by its DLT data, or alone by its response data, for the rest of the
//     trial;
//   - the OBD is the admissible dose up to the estimated MTD with the highest
//     utility, from the data at the last decision.

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "tite_boinet_core.h"

namespace simfastboin {

// Stop codes returned by simulate_bf_et_one().
enum BfEtStop {
  BF_STOP_NO_ADMISSIBLE = 0,     // every dose is excluded (an early stop)
  BF_STOP_LOWEST_EXCLUDED = 1,   // de-escalation from an excluded lowest dose
  BF_STOP_TOP_EXCLUDED = 2,      // escalation from an excluded highest dose
  BF_STOP_NO_LOWER = 3,          // no admissible dose at or below the target
  BF_STOP_N_EARLYSTOP = 4,       // n_stop patients at the current dose
  BF_STOP_MAX_COHORTS = 5        // all cohorts were enrolled
};

// The settings of a BF-BOIN-ET trial that do not change from trial to trial.
struct BfEtDesign {
  std::vector<double> p_true;
  std::vector<double> p_resp;
  double resp_cor;
  std::vector<int> cohort_size;
  int start_dose;                // 1-based
  double lambda1;
  double lambda2;
  double eta1;
  int n_elim;                    // patients needed before a dose is excluded
  int n_stop;                    // patients at a dose that end the trial
  double suspend;                // share of pending DLT outcomes allowed
  double extension;              // backfilling after the last cohort
  std::vector<int> tox_ok;       // tables as in boinet_core.h
  std::vector<int> resp_ok;
  int max_n;
  double target;                 // target DLT probability (OBD selection)
  double psi00, psi01, psi10, psi11;  // scores of the utility
};

// Utility of a dose with response rate pe and DLT rate pt, rounded to ten
// decimals (et_round10() of boinet_core.h), so that doses with the same
// utility as fractions tie on every platform, in the choice of the backfill
// dose and in the selection of the OBD.
inline double bf_utility(double pe, double pt, const BfEtDesign& des) {
  return et_round10(des.psi01 * (1.0 - pe) * pt +
                    des.psi00 * (1.0 - pe) * (1.0 - pt) +
                    des.psi11 * pe * pt + des.psi10 * pe * (1.0 - pt));
}

// Sample quantile of type 7 (the default of R) of x at probability prob.
inline double quantile_type7(std::vector<double> x, double prob) {
  std::sort(x.begin(), x.end());
  const int n = static_cast<int>(x.size());
  const double h = (n - 1) * prob;
  const int lo = static_cast<int>(std::floor(h));
  if (lo >= n - 1) return x[n - 1];
  return x[lo] + (h - lo) * (x[lo + 1] - x[lo]);
}

// OBD among the evaluated doses (0-based positions) from the isotonic DLT
// estimates probt, the response estimates probe and the admissibility cand,
// as obd.select() inside get.oc.backboinet(); -1 when none.
inline int bf_obd_choose(const std::vector<double>& probt,
                         const std::vector<double>& probe,
                         const std::vector<char>& cand, const BfEtDesign& des) {
  const int k_doses = static_cast<int>(probt.size());
  double mdif = 0.0;
  bool first = true;
  for (int k = 0; k < k_doses; ++k) {
    if (!cand[k]) continue;
    const double dist = std::fabs(probt[k] - des.target);
    if (first || dist < mdif) {
      mdif = dist;
      first = false;
    }
  }
  if (first) return -1;
  int mtd = 0;
  for (int k = 0; k < k_doses; ++k) {
    if (std::fabs(probt[k] - des.target) == mdif) mtd = k;
  }
  double best = 0.0;
  first = true;
  for (int k = 0; k <= mtd; ++k) {
    if (!cand[k]) continue;
    const double u = bf_utility(probe[k], probt[k], des);
    if (first || u > best) {
      best = u;
      first = false;
    }
  }
  if (first) return -1;
  for (int k = 0; k <= mtd; ++k) {
    if (cand[k] && bf_utility(probe[k], probt[k], des) == best) return k;
  }
  return -1;
}

// Simulate a single BF-BOIN-ET trial and select its OBD.
//
// nn, nb, yy and xx receive the patients, the backfilled patients, and the
// DLTs and responses observed at the last decision at every dose; yb and xb
// the DLTs and responses among the backfilled patients. obd is 1-based, or 0
// when none is selected; early_stop is the early stop of the reference, when
// every dose is excluded. duration is the time at which the last patient
// finishes both assessments.
//
// For every cohort, rng() is called for the times between the arrivals of its
// patients, twice for each of its patients in arrival order, then for the
// arrivals until the next cohort, then twice for each backfilled patient.
template <class Rng>
inline void simulate_bf_et_one(const BfEtDesign& des, const EtTiming& tm,
                               Rng& rng, std::vector<int>& nn,
                               std::vector<int>& nb, std::vector<int>& yy,
                               std::vector<int>& xx, std::vector<int>& yb,
                               std::vector<int>& xb, int& cohorts_used,
                               int& stop_code, bool& early_stop, int& obd,
                               double& duration) {
  const int n_doses = static_cast<int>(des.p_true.size());
  const int n_cohorts = static_cast<int>(des.cohort_size.size());
  const int w = des.max_n + 1;
  std::fill(nn.begin(), nn.end(), 0);
  std::fill(nb.begin(), nb.end(), 0);
  std::fill(yy.begin(), yy.end(), 0);
  std::fill(xx.begin(), xx.end(), 0);
  std::fill(yb.begin(), yb.end(), 0);
  std::fill(xb.begin(), xb.end(), 0);
  std::vector<double> pt(n_doses, 0.0), pe(n_doses, 0.0);
  std::vector<char> elim_tox(n_doses, 0), elim_eff(n_doses, 0), adm(n_doses);
  std::vector<int> bdose;        // backfill doses (0-based), ascending
  std::vector<EtPatient> pts;

  int n_esc_max = 0;
  for (int i = 0; i < n_cohorts; ++i) n_esc_max += des.cohort_size[i];
  int n_esc = 0;
  double t_curr = 0.0;
  int cur = des.start_dose;
  cohorts_used = 0;
  stop_code = BF_STOP_MAX_COHORTS;
  early_stop = false;

  for (int i = 0; i < n_cohorts; ++i) {
    ++cohorts_used;
    const int m = des.cohort_size[i];
    const int d = cur - 1;

    // The dose-escalation cohort.
    std::vector<double> enter(m);
    enter[0] = t_curr;
    for (int k = 1; k < m; ++k) enter[k] = enter[k - 1] + et_gap(tm, rng);
    double max_finish = -std::numeric_limits<double>::infinity();
    double max_onset = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < m; ++k) {
      EtPatient p = et_new_patient(d, enter[k], des.p_true, des.p_resp,
                                   des.resp_cor, tm, rng);
      max_finish = std::max(max_finish, p.enter + tm.window);
      max_onset = std::max(max_onset, p.enter + p.t_tox);
      pts.push_back(p);
    }
    double q50 = std::numeric_limits<double>::infinity();
    if (des.suspend != 0.0) {
      std::vector<double> endt;
      for (std::size_t q = 0; q < pts.size(); ++q) {
        if (pts[q].dose == d) endt.push_back(pts[q].end_tox);
      }
      const double qv = quantile_type7(endt, 1.0 - des.suspend);
      q50 = -std::numeric_limits<double>::infinity();
      for (std::size_t q = 0; q < endt.size(); ++q) {
        if (endt[q] <= qv && endt[q] > q50) q50 = endt[q];
      }
    }
    nn[d] += m;
    n_esc += m;

    // Arrivals until the next cohort.
    t_curr = -std::numeric_limits<double>::infinity();
    for (std::size_t q = 0; q < pts.size(); ++q) {
      t_curr = std::max(t_curr, pts[q].enter);
    }
    double bench = std::min(max_finish, max_onset);
    if (des.suspend != 0.0) bench = std::min(bench, q50);
    if (i == n_cohorts - 1) bench += des.extension;
    std::vector<double> queue;
    do {
      t_curr += et_gap(tm, rng);
      queue.push_back(t_curr);
    } while (bench > t_curr);
    queue.pop_back();

    // Backfilling.
    if (cur != 1 && !bdose.empty()) {
      for (std::size_t qi = 0; qi < queue.size(); ++qi) {
        const double tq = queue[qi];
        std::vector<int> eva;
        for (std::size_t b = 0; b < bdose.size(); ++b) {
          const int j = bdose[b];
          if (nn[j] != 0 && nn[j] < des.n_stop) eva.push_back(j);
        }
        if (eva.empty()) break;
        if (qi != 0) et_tite_rates(pts, tq, tm, nn, yy, xx, pt, pe);
        std::vector<double> est(eva.size());
        for (std::size_t k = 0; k < eva.size(); ++k) est[k] = pt[eva[k]];
        pava_equal(est);
        double best = 0.0;
        int s = -1;
        for (std::size_t k = 0; k < eva.size(); ++k) {
          const double u = bf_utility(pe[eva[k]], est[k], des);
          if (s < 0 || u >= best) {
            best = u;
            s = eva[k];
          }
        }
        EtPatient p = et_new_patient(s, tq, des.p_true, des.p_resp,
                                     des.resp_cor, tm, rng);
        p.backfill = 1;
        pts.push_back(p);
        nn[s] += 1;
        nb[s] += 1;
      }
    }

    // Decision for the next cohort.
    et_tite_rates(pts, t_curr, tm, nn, yy, xx, pt, pe);
    std::fill(yb.begin(), yb.end(), 0);
    std::fill(xb.begin(), xb.end(), 0);
    for (std::size_t q = 0; q < pts.size(); ++q) {
      const EtPatient& p = pts[q];
      if (!p.backfill) continue;
      if (p.end_tox <= t_curr && p.dlt) yb[p.dose] += 1;
      if (p.end_resp <= t_curr && p.resp) xb[p.dose] += 1;
    }
    const double t = pt[d];
    int nxt = cur;
    int top = cur - 1;           // highest backfill dose, 1-based
    if (t <= des.lambda1) {
      nxt = cur + 1;
      top = cur;
    } else if ((t > des.lambda1 && t <= des.lambda2) ||
               (yy[d] == 1 && nn[d] == 3)) {
      nxt = cur;
      top = cur - 1;
    } else if (t > des.lambda2) {
      nxt = cur - 1;
      top = cur - 2;
    }
    int jmin = 0;
    for (int j = 0; j < n_doses; ++j) {
      if (pe[j] >= des.eta1) {
        jmin = j + 1;
        break;
      }
    }

    // Exclusion, kept for the rest of the trial.
    for (int j = 0; j < n_doses; ++j) {
      if (nn[j] >= des.n_elim && des.tox_ok[nn[j] * w + yy[j]] != 1) {
        for (int k = j; k < n_doses; ++k) elim_tox[k] = 1;
        break;
      }
    }
    for (int j = 0; j < n_doses; ++j) {
      if (nn[j] >= des.n_elim && des.resp_ok[nn[j] * w + xx[j]] != 1) {
        elim_eff[j] = 1;
      }
    }
    bool any_adm = false;
    for (int j = 0; j < n_doses; ++j) {
      adm[j] = (elim_tox[j] == 0 && elim_eff[j] == 0) ? 1 : 0;
      if (adm[j]) any_adm = true;
    }

    if (!any_adm) {
      early_stop = true;
      stop_code = BF_STOP_NO_ADMISSIBLE;
      break;
    }
    if (n_esc == n_esc_max) {
      stop_code = BF_STOP_MAX_COHORTS;
      break;
    }
    if (nn[d] >= des.n_stop &&
        ((t > des.lambda1 && t <= des.lambda2) ||
         (cur == 1 && t >= des.lambda2) ||
         ((cur == n_doses || !adm[cur]) && t <= des.lambda1))) {
      stop_code = BF_STOP_N_EARLYSTOP;
      break;
    }
    if (nxt == 0) {
      if (!adm[0]) {
        stop_code = BF_STOP_LOWEST_EXCLUDED;
        break;
      }
      cur = 1;
    } else if (nxt == n_doses + 1) {
      if (!adm[n_doses - 1]) {
        stop_code = BF_STOP_TOP_EXCLUDED;
        break;
      }
      cur = n_doses;
    } else if (adm[nxt - 1]) {
      cur = nxt;
    } else if (cur < nxt) {
      for (int j = nxt; j <= n_doses; ++j) {
        if (adm[j - 1]) {
          cur = j;
          break;
        }
      }
    } else {
      int lower = 0;
      for (int j = nxt; j >= 1; --j) {
        if (adm[j - 1]) {
          lower = j;
          break;
        }
      }
      if (lower == 0) {
        stop_code = BF_STOP_NO_LOWER;
        break;
      }
      cur = lower;
    }
    bdose.clear();
    if (jmin >= 1) {
      for (int j = jmin; j <= top; ++j) {
        if (adm[j - 1]) bdose.push_back(j - 1);
      }
    }
  }

  duration = -std::numeric_limits<double>::infinity();
  for (std::size_t q = 0; q < pts.size(); ++q) {
    duration = std::max(duration, std::max(pts[q].end_tox, pts[q].end_resp));
  }

  // OBD from the data at the last decision.
  obd = 0;
  if (static_cast<int>(pts.size()) < n_esc_max && early_stop) return;
  std::vector<int> eva;
  for (int j = 0; j < n_doses; ++j) {
    if (nn[j] != 0) eva.push_back(j);
  }
  const int k_doses = static_cast<int>(eva.size());
  std::vector<char> cand(k_doses);
  bool any_cand = false;
  for (int k = 0; k < k_doses; ++k) {
    const int j = eva[k];
    cand[k] = (des.tox_ok[nn[j] * w + yy[j]] == 1 &&
               des.resp_ok[nn[j] * w + xx[j]] == 1) ? 1 : 0;
    if (cand[k]) any_cand = true;
  }
  if (k_doses == 1) {
    if (cand[0]) obd = eva[0] + 1;
    return;
  }
  if (!any_cand) return;
  std::vector<double> probt(k_doses), probe(k_doses);
  for (int k = 0; k < k_doses; ++k) {
    const int j = eva[k];
    probt[k] = static_cast<double>(yy[j]) / static_cast<double>(nn[j]);
    probe[k] = static_cast<double>(xx[j]) / static_cast<double>(nn[j]);
  }
  pava_equal(probt);
  const int k = bf_obd_choose(probt, probe, cand, des);
  if (k >= 0) obd = eva[k] + 1;
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_BF_BOINET_CORE_H
