#ifndef SIMFASTBOIN_TITE_COMB_BOIN12_CORE_H
#define SIMFASTBOIN_TITE_COMB_BOIN12_CORE_H

// Core algorithm of the TITE-Comb-BOIN12 design (Lu, Zhang, Yuan and Lin,
// 2025): the rules of Comb-BOIN12 (comb_boin12_core.h) applied at the
// decision times of TITE-BOIN12 (tite_boin12_core.h), with the interim
// quantities of TITE-BOIN12 (tb12_interim()) at every combination.
//
// Patients arrive and are followed as in tb12_simulate_one(): the first
// patient of a cohort arrives at the time of the previous decision and the
// next decision is due when the next patient arrives; while fewer than a
// share min_share of the patients at the current combination have finished
// both assessments, the decision waits for the next assessment to end; after
// the last cohort it waits until every patient has finished both
// assessments.
//
// At a decision, with the interim quantities of every combination:
//   - a combination is admissible when Pr(p_T > target) < cutoff_eli (the
//     cutoff of tox_ok_lowest, cutoff_lowest, at (0, 0)) under
//     Beta(1 + DLTs, 1 + completed without DLT + sum of w) and Pr(p_E <
//     target_resp) < cutoff_eli_resp under the same posterior of the
//     responses; a combination excluded by its DLT data also excludes every
//     combination with both levels at least its own;
//   - the desirability is that of tb12_desirability();
//   - the next combination follows cb12_decide() with phat = pi_hat of the
//     DLTs at the current combination and the patients treated.
// With pending_use = TCB12_PENDING_UTILITY (the rules as Sections 2.2 and
// 3.1 of Lu et al. state them) the pending outcomes enter only the
// quasi-events of the desirability: the admissibility uses Beta(1 + DLTs,
// 1 + n - DLTs) and Beta(1 + responses, 1 + n - responses), n the patients
// treated, and phat is the observed DLTs over the patients treated, so that
// a pending outcome counts as no event.
// With every outcome observed these are exactly the quantities of
// Comb-BOIN12, so that with min_share = 1 and fixed arrivals the trials are
// those of cb12_simulate_one().
//
// The trial stops at once when no combination is admissible, or when
// stop_no_candidate stops it. When it stops because the current combination
// has n_earlystop patients, or after the last cohort, it waits until every
// patient has finished both assessments; the duration is that time, and the
// counts returned are the final outcomes of all patients.

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "comb_boin12_core.h"
#include "tite_boin12_core.h"

namespace simfastboin {

// Use of the pending outcomes.
enum Tcb12PendingUse {
  TCB12_PENDING_ALL = 0,       // estimates, admissibility and utility
  TCB12_PENDING_UTILITY = 1    // the quasi-events of the utility only
};

// The settings of a TITE-Comb-BOIN12 trial. The admissibility tables of cb
// are not used; cutoff_lowest is the cutoff of the DLT rule at (0, 0).
struct Tcb12Design {
  Cb12Design cb;
  double target;
  double target_resp;
  double cutoff_eli;
  double cutoff_lowest;
  double cutoff_eli_resp;
  int pending;                   // Tb12Pending
  int pending_use;               // Tcb12PendingUse
};

// Admissibility of every combination from the interim quantities. Returns
// whether any combination is admissible.
template <class Pb>
inline bool tcb12_admissible_set(const std::vector<Tb12Interim>& est,
                                 const Tcb12Design& des, const Pb& pb,
                                 char* adm) {
  const int J = des.cb.n_rows;
  const int K = des.cb.n_cols;
  bool any_adm = false;
  std::vector<char> blocked(J * K, 0);
  for (int b = 0; b < K; ++b) {
    for (int a = 0; a < J; ++a) {
      const int c = a + J * b;
      const Tb12Interim& e = est[c];
      const double cut = c == 0 ? des.cutoff_lowest : des.cutoff_eli;
      const bool all = des.pending_use == TCB12_PENDING_ALL;
      const double rest_t = all ? e.ne_t + e.s_t : e.n - e.ev_t;
      const double rest_e = all ? e.ne_e + e.s_e : e.n - e.ev_e;
      bool blk = !(pb(des.target, 1.0 + e.ev_t, 1.0 + rest_t, false) < cut);
      if (a > 0 && blocked[c - 1]) blk = true;
      if (b > 0 && blocked[c - J]) blk = true;
      blocked[c] = blk ? 1 : 0;
      const bool r_ok = pb(des.target_resp, 1.0 + e.ev_e, 1.0 + rest_e,
                           true) < des.cutoff_eli_resp;
      const bool ok = !blk && r_ok;
      adm[c] = ok ? 1 : 0;
      if (ok) any_adm = true;
    }
  }
  return any_adm;
}

// Simulate a single TITE-Comb-BOIN12 trial.
//
// nn, yy, xx and zz (J * K each) receive the patients, DLTs, responses and
// patients with both at every combination, from the final outcomes of all
// patients. stop_code is a B12Stop code and duration the end of the trial.
//
// rng() is called as in tb12_simulate_one(): for every cohort, for the times
// between the arrivals of its patients, then for the time to the next
// decision (except after the last cohort), then twice for each patient in
// arrival order, and once more at a decision when a tie remains.
template <class Rng, class Pb>
inline void tcb12_simulate_one(const Tcb12Design& des, const EtTiming& tm,
                               Rng& rng, const Pb& pb,
                               std::vector<Tb12Patient>& pts,
                               std::vector<int>& nn, std::vector<int>& yy,
                               std::vector<int>& xx, std::vector<int>& zz,
                               int& cohorts_used, int& stop_code,
                               double& duration) {
  const Cb12Design& cb = des.cb;
  const B12Design& b = cb.b12;
  const int n_cells = cb.n_rows * cb.n_cols;
  const int n_cohorts = static_cast<int>(b.cohort_size.size());
  std::fill(nn.begin(), nn.end(), 0);
  pts.clear();
  std::vector<Tb12Interim> est(n_cells);
  std::vector<double> desir(n_cells);
  std::vector<char> adm(n_cells);

  double t_dec = 0.0;
  int cur = cb.start_row + cb.n_rows * cb.start_col;
  cohorts_used = 0;
  stop_code = B12_STOP_MAX_COHORTS;
  bool wait_all = true;

  for (int i = 0; i < n_cohorts; ++i) {
    ++cohorts_used;
    const int m = b.cohort_size[i];
    const int d = cur;
    const bool last = i == n_cohorts - 1;

    std::vector<double> enter(m);
    enter[0] = t_dec;
    for (int k = 1; k < m; ++k) enter[k] = enter[k - 1] + et_gap(tm, rng);
    if (!last) t_dec = enter[m - 1] + et_gap(tm, rng);
    for (int k = 0; k < m; ++k) {
      pts.push_back(tb12_new_patient(b, tm, d, enter[k], rng));
    }
    nn[d] += m;
    if (last) t_dec = std::max(enter[m - 1], tb12_all_done(pts));

    // Wait until a share min_share of the patients at the current
    // combination have finished both assessments.
    for (;;) {
      int done = 0;
      double next = std::numeric_limits<double>::infinity();
      for (std::size_t q = 0; q < pts.size(); ++q) {
        const Tb12Patient& p = pts[q];
        if (p.dose != d) continue;
        if (p.end_tox <= t_dec && p.end_resp <= t_dec) ++done;
        if (p.end_tox > t_dec && p.end_tox < next) next = p.end_tox;
        if (p.end_resp > t_dec && p.end_resp < next) next = p.end_resp;
      }
      if (static_cast<double>(done) / nn[d] >= tm.min_share) break;
      t_dec = next;
    }

    for (int c = 0; c < n_cells; ++c) {
      tb12_interim(pts, c, t_dec, tm, b, des.pending, est[c]);
    }
    if (!tcb12_admissible_set(est, des, pb, &adm[0])) {
      stop_code = B12_STOP_NO_ADMISSIBLE;
      wait_all = false;
      break;
    }
    for (int c = 0; c < n_cells; ++c) {
      desir[c] = tb12_desirability(est[c], b, pb);
    }
    auto value = [&](int c) { return desir[c]; };
    const double phat = des.pending_use == TCB12_PENDING_ALL ? est[d].pi_t :
      est[d].ev_t / nn[d];
    const int nxt = cb12_decide(cur, &nn[0], phat, &adm[0], value, cb, rng);
    if (nxt == CB12_NEXT_STOP) {
      stop_code = B12_STOP_NO_CANDIDATE;
      wait_all = false;
      break;
    }
    if (b.n_earlystop > 0 && nn[d] >= b.n_earlystop &&
        (b.early_stop_simple || nxt == cur)) {
      stop_code = B12_STOP_N_EARLYSTOP;
      break;
    }
    cur = nxt;
  }

  duration = wait_all ? std::max(t_dec, tb12_all_done(pts)) : t_dec;
  std::fill(yy.begin(), yy.end(), 0);
  std::fill(xx.begin(), xx.end(), 0);
  std::fill(zz.begin(), zz.end(), 0);
  for (std::size_t i = 0; i < pts.size(); ++i) {
    const Tb12Patient& p = pts[i];
    if (p.dlt) yy[p.dose] += 1;
    if (p.resp) xx[p.dose] += 1;
    if (p.dlt && p.resp) zz[p.dose] += 1;
  }
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_TITE_COMB_BOIN12_CORE_H
