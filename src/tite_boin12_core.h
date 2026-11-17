#ifndef SIMFASTBOIN_TITE_BOIN12_CORE_H
#define SIMFASTBOIN_TITE_BOIN12_CORE_H

// Core algorithm of the TITE-BOIN12 design (Zhou, Lin, Lee, Li, Wang, Li and
// Yuan, 2022) with the approximated likelihood, written in plain C++ like
// boin12_core.h so that it can be verified outside of R as well.
//
// Patients arrive over time as in tite_boinet_core.h: the first patient of a
// cohort arrives at the time of the previous decision and the next decision
// is due when the next patient arrives. While fewer than a share min_share of
// the patients at the current dose have finished both assessments, the
// decision waits for the next assessment to end. After the last cohort the
// decision waits until every patient has finished both assessments.
//
// At a decision at time t, for every dose with patients:
//   - an outcome is observed once its event has occurred or its window has
//     ended; otherwise it is pending, and w = (t - arrival) / window;
//   - pi_hat = observed events / (observed events + completed without the
//     event + sum of w over the pending patients) (Supplement S.7, Eq. 14);
//   - a pending outcome counts as E(Y = 1 | pending), which is
//     pi_hat (1 - w) / (1 - pi_hat w) with pending = conditional (Eq. 11) or
//     pi_hat itself with pending = marginal; the DLT and the response of a
//     patient are taken as independent for this purpose only;
//   - the quasi-events x are the sum of the expected utilities / 100, the
//     patients with both outcomes observed counting as in BOIN12, and the
//     desirability is Pr(Beta(alpha + x, beta + n - x) > u_benchmark) with n
//     the patients treated;
//   - the dose is admissible when Pr(p_T > target) < cutoff_eli under
//     Beta(1 + DLTs, 1 + completed without DLT + sum of w) and
//     Pr(p_E < target_resp) < cutoff_eli_resp under the same posterior of the
//     responses; a dose excluded by its DLT data also excludes every higher
//     dose.
// The next dose follows b12_next_dose() with phat = pi_hat of the DLTs at the
// current dose. With every outcome observed, all of these are exactly the
// quantities of BOIN12.
//
// The trial stops at once when no dose is admissible, or when
// stop_no_candidate stops it (see b12_next_dose()). When it stops because
// the current dose has n_earlystop patients, or after the last cohort, it
// waits until every patient has finished both assessments; the duration is
// that time, and the counts returned are the final outcomes of all patients.
//
// Each patient draws two uniform variates u and v exactly as in
// b12_patient(): the DLT occurs when u < p_true and the response follows from
// u and v by the copula or by the binary link. The time to a DLT is the
// quantile u of its distribution; the time to a response is the quantile of
// the copula variate, or, with the binary link, of (v / q) p_resp, which is
// uniform given a response.

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "boin12_core.h"
#include "tite_boinet_core.h"

namespace simfastboin {

// Value that replaces a pending outcome.
enum Tb12Pending {
  TB12_PENDING_CONDITIONAL = 0,  // pi_hat (1 - w) / (1 - pi_hat w)
  TB12_PENDING_MARGINAL = 1      // pi_hat
};

// The settings of a TITE-BOIN12 trial. The admissibility tables of b12 are
// not used.
struct Tb12Design {
  B12Design b12;
  double target;
  double target_resp;
  double cutoff_eli;
  double cutoff_eli_resp;
  int pending;                   // Tb12Pending
};

// A patient of a TITE-BOIN12 trial.
struct Tb12Patient {
  int dose;          // 0-based
  double enter;      // arrival time
  double end_tox;    // time at which the DLT is observed or the window ends
  double end_resp;   // the same for the response
  char dlt;
  char resp;
};

// Interim quantities of one dose.
struct Tb12Interim {
  int n;             // patients treated
  double ev_t;       // DLTs observed
  double ne_t;       // patients who completed the DLT window without DLT
  double s_t;        // sum of w over the patients with a pending DLT outcome
  double ev_e;       // the same for the response
  double ne_e;
  double s_e;
  double pi_t;
  double pi_e;
  double q;          // quasi-events
};

// Draw a patient arriving at time t at dose d (0-based).
template <class Rng>
inline Tb12Patient tb12_new_patient(const B12Design& des, const EtTiming& tm,
                                    int d, double t, Rng& rng) {
  const double u = rng();
  const double v = rng();
  bool dlt, resp;
  b12_patient(des, d, u, v, dlt, resp);
  Tb12Patient pt;
  pt.dose = d;
  pt.enter = t;
  pt.dlt = dlt ? 1 : 0;
  pt.resp = resp ? 1 : 0;
  double tt = tm.window;
  if (dlt) {
    tt = std::min(et_event_time(u, des.p_true[d], tm.shape_t[d], tm.scale_t[d],
                                tm.window, tm.event_time), tm.window);
  }
  double te = tm.resp_window;
  if (resp) {
    double q;
    if (des.cor_type == ET_COR_BINARY) {
      const double c = dlt ? des.q1[d] : des.q0[d];
      q = (v / c) * des.p_resp[d];
    } else {
      q = response_variate(u, v, des.resp_cor);
    }
    te = std::min(et_event_time(q, des.p_resp[d], tm.shape_e[d], tm.scale_e[d],
                                tm.resp_window, tm.event_time), tm.resp_window);
  }
  pt.end_tox = t + tt;
  pt.end_resp = t + te;
  return pt;
}

// E(Y = 1 | pending) for an estimate pi and a followed fraction w.
inline double tb12_pending_value(double pi, double w, int pending) {
  if (pending == TB12_PENDING_MARGINAL) return pi;
  const double den = 1.0 - pi * w;
  if (den <= 0.0) return pi;
  return pi * (1.0 - w) / den;
}

// Interim quantities of dose j (0-based) at time t from the patients pts.
inline void tb12_interim(const std::vector<Tb12Patient>& pts, int j, double t,
                         const EtTiming& tm, const B12Design& des, int pending,
                         Tb12Interim& est) {
  est.n = 0;
  est.ev_t = 0.0; est.ne_t = 0.0; est.s_t = 0.0;
  est.ev_e = 0.0; est.ne_e = 0.0; est.s_e = 0.0;
  for (std::size_t i = 0; i < pts.size(); ++i) {
    const Tb12Patient& p = pts[i];
    if (p.dose != j) continue;
    ++est.n;
    if (p.end_tox <= t) {
      if (p.dlt) est.ev_t += 1.0; else est.ne_t += 1.0;
    } else {
      est.s_t += (t - p.enter) / tm.window;
    }
    if (p.end_resp <= t) {
      if (p.resp) est.ev_e += 1.0; else est.ne_e += 1.0;
    } else {
      est.s_e += (t - p.enter) / tm.resp_window;
    }
  }
  const double ess_t = est.ev_t + est.ne_t + est.s_t;
  const double ess_e = est.ev_e + est.ne_e + est.s_e;
  est.pi_t = ess_t > 0.0 ? est.ev_t / ess_t : 0.0;
  est.pi_e = ess_e > 0.0 ? est.ev_e / ess_e : 0.0;

  // Patients with both outcomes observed count as in BOIN12; the others add
  // their expected utility. The sum is divided by 100 only at the end, so
  // that utilities that are whole numbers add up exactly and doses with equal
  // data have equal desirabilities.
  int nc = 0, yc = 0, xc = 0, zc = 0;
  double extra = 0.0;
  for (std::size_t i = 0; i < pts.size(); ++i) {
    const Tb12Patient& p = pts[i];
    if (p.dose != j) continue;
    const bool obs_t = p.end_tox <= t;
    const bool obs_e = p.end_resp <= t;
    if (obs_t && obs_e) {
      ++nc;
      if (p.dlt) ++yc;
      if (p.resp) ++xc;
      if (p.dlt && p.resp) ++zc;
      continue;
    }
    const double pt1 = obs_t ? static_cast<double>(p.dlt) :
      tb12_pending_value(est.pi_t, (t - p.enter) / tm.window, pending);
    const double pe1 = obs_e ? static_cast<double>(p.resp) :
      tb12_pending_value(est.pi_e, (t - p.enter) / tm.resp_window, pending);
    extra += 100.0 * (1.0 - pt1) * pe1 + des.psi00 * (1.0 - pt1) * (1.0 - pe1) +
      des.psi11 * pt1 * pe1;
  }
  // The expression of b12_quasi_events(), with the pending patients added
  // before the division; with none pending it is exactly that of BOIN12.
  const double n_c = nc, y_c = yc, x_c = xc, z_c = zc;
  est.q = (100.0 * (x_c - z_c) + des.psi00 * (n_c - y_c - x_c + z_c) +
           des.psi11 * z_c + extra) / 100.0;
}

// Desirability of a dose from its interim quantities. pb(q, a, b, lower) is
// the beta distribution function, lower or upper tail.
template <class Pb>
inline double tb12_desirability(const Tb12Interim& est, const B12Design& des,
                                const Pb& pb) {
  return pb(des.u_benchmark, des.alpha + est.q, des.beta + est.n - est.q,
            false);
}

// Admissibility of every dose from the interim quantities. Returns whether
// any dose is admissible.
template <class Pb>
inline bool tb12_admissible_set(const std::vector<Tb12Interim>& est,
                                const Tb12Design& des, const Pb& pb,
                                std::vector<char>& adm) {
  const int n_doses = static_cast<int>(est.size());
  bool tox_out = false;
  bool any_adm = false;
  for (int j = 0; j < n_doses; ++j) {
    const Tb12Interim& e = est[j];
    const bool t_ok = pb(des.target, 1.0 + e.ev_t, 1.0 + e.ne_t + e.s_t,
                         false) < des.cutoff_eli;
    if (!t_ok) tox_out = true;
    const bool r_ok = pb(des.target_resp, 1.0 + e.ev_e, 1.0 + e.ne_e + e.s_e,
                         true) < des.cutoff_eli_resp;
    const bool ok = t_ok && !tox_out && r_ok;
    adm[j] = ok ? 1 : 0;
    if (ok) any_adm = true;
  }
  return any_adm;
}

// Next dose (1-based) from the current dose cur by the rules of
// b12_next_dose(), with the DLT rate phat of the current dose and the
// desirability of every dose; 0 when no dose is admissible and B12_NEXT_STOP
// when stop_no_candidate stops the trial. rng() gives the variate of a random
// tie break.
template <class Rng>
inline int tb12_next_dose(int cur, const std::vector<int>& nn, double phat,
                          const std::vector<double>& desir,
                          const std::vector<char>& adm,
                          const B12Design& des, Rng& rng) {
  const int n_doses = static_cast<int>(nn.size());
  const int d = cur - 1;

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
      if (desir[j] > best_value) {
        best_value = desir[j];
        tied[0] = j;
        n_tied = 1;
      } else if (desir[j] == best_value) {
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

// Time at which every patient has finished both assessments.
inline double tb12_all_done(const std::vector<Tb12Patient>& pts) {
  double t = 0.0;
  for (std::size_t i = 0; i < pts.size(); ++i) {
    t = std::max(t, std::max(pts[i].end_tox, pts[i].end_resp));
  }
  return t;
}

// Simulate a single TITE-BOIN12 trial.
//
// nn, yy, xx and zz receive the patients, DLTs, responses and patients with
// both at every dose, from the final outcomes of all patients. stop_code is a
// B12Stop code and duration the end of the trial.
//
// For every cohort, rng() is called for the times between the arrivals of its
// patients, then for the time to the next decision (except after the last
// cohort), then twice for each patient in arrival order, as in
// simulate_tite_et_one(), and with tie_random once more at a decision between
// tied doses. With arrivals fixed no variate is drawn for the times, so that
// with min_share = 1 the trials are those of b12_simulate_one().
template <class Rng, class Pb>
inline void tb12_simulate_one(const Tb12Design& des, const EtTiming& tm,
                              Rng& rng, const Pb& pb,
                              std::vector<Tb12Patient>& pts,
                              std::vector<int>& nn, std::vector<int>& yy,
                              std::vector<int>& xx, std::vector<int>& zz,
                              int& cohorts_used, int& stop_code,
                              double& duration) {
  const B12Design& b = des.b12;
  const int n_doses = static_cast<int>(b.p_true.size());
  const int n_cohorts = static_cast<int>(b.cohort_size.size());
  std::fill(nn.begin(), nn.end(), 0);
  pts.clear();
  std::vector<Tb12Interim> est(n_doses);
  std::vector<double> desir(n_doses);
  std::vector<char> adm(n_doses);

  double t_dec = 0.0;
  int cur = b.start_dose;
  cohorts_used = 0;
  stop_code = B12_STOP_MAX_COHORTS;
  bool wait_all = true;

  for (int i = 0; i < n_cohorts; ++i) {
    ++cohorts_used;
    const int m = b.cohort_size[i];
    const int d = cur - 1;
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

    // Wait until a share min_share of the patients at the current dose have
    // finished both assessments.
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

    for (int j = 0; j < n_doses; ++j) {
      tb12_interim(pts, j, t_dec, tm, b, des.pending, est[j]);
    }
    if (!tb12_admissible_set(est, des, pb, adm)) {
      stop_code = B12_STOP_NO_ADMISSIBLE;
      wait_all = false;
      break;
    }
    for (int j = 0; j < n_doses; ++j) {
      desir[j] = tb12_desirability(est[j], b, pb);
    }
    const int nxt = tb12_next_dose(cur, nn, est[d].pi_t, desir, adm, b, rng);
    if (nxt == 0) {
      stop_code = B12_STOP_NO_ADMISSIBLE;
      wait_all = false;
      break;
    }
    if (nxt == B12_NEXT_STOP) {
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

#endif  // SIMFASTBOIN_TITE_BOIN12_CORE_H
