#ifndef SIMFASTBOIN_TITE_BOINET_CORE_H
#define SIMFASTBOIN_TITE_BOINET_CORE_H

// Core algorithm of the TITE-BOIN-ET design (Takeda, Morita and Taguri,
// 2020), written in plain C++ like boinet_core.h so that it can be verified
// outside of R as well.
//
// The rules follow tite.boinet() of the boinet package (version 1.6.0), with
// the two changes of boinet_core.h unless compat is set:
//   - the first patient of a cohort arrives at the time of the previous
//     decision, and the next decision is due when the next patient arrives,
//     or, after the last cohort, once every patient has been followed for both
//     assessment windows;
//   - while fewer than half of the patients at the current dose have finished
//     both assessments, the decision waits for the next assessment to end;
//   - the DLT and response rates of every dose count a pending patient by the
//     fraction of the window followed so far, and the decision rule of
//     BOIN-ET is applied to these rates;
//   - the admissibility of a dose counts the pending patients as patients
//     without the event, and the trial moves between doses as in
//     simulate_et_one().
// The event times come from Weibull or uniform distributions, linked by a
// Gaussian copula, so that the probability of an event within its window is
// the DLT or the response probability of the dose.

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "boinet_core.h"

namespace simfastboin {

// Distribution of the times between arrivals.
enum EtArrival {
  ET_ARRIVAL_UNIFORM = 0,      // uniform on (0, 2 / accrual_rate)
  ET_ARRIVAL_EXPONENTIAL = 1,  // exponential with rate accrual_rate
  ET_ARRIVAL_FIXED = 2         // exactly 1 / accrual_rate
};

// Distribution of the time to an event.
enum EtEventTime {
  ET_TIME_WEIBULL = 0,
  ET_TIME_UNIFORM = 1          // uniform on (0, window / probability)
};

// The settings of the time-to-event part of a trial.
struct EtTiming {
  double window;               // DLT assessment window
  double resp_window;          // response assessment window
  double accrual_rate;         // patients per unit of time
  int arrival;                 // EtArrival
  int event_time;              // EtEventTime
  std::vector<double> shape_t; // Weibull shape and scale of the DLT time
  std::vector<double> scale_t;
  std::vector<double> shape_e; // and of the response time
  std::vector<double> scale_e;
};

// Time between two arrivals. Draws one uniform variate unless arrivals are
// fixed.
template <class Rng>
inline double et_gap(const EtTiming& tm, Rng& rng) {
  if (tm.arrival == ET_ARRIVAL_FIXED) return 1.0 / tm.accrual_rate;
  const double u = rng();
  if (tm.arrival == ET_ARRIVAL_EXPONENTIAL) return -std::log(u) / tm.accrual_rate;
  return u * 2.0 / tm.accrual_rate;
}

// Time to an event of probability p within the window from the uniform
// variate u; infinite when p is zero.
inline double et_event_time(double u, double p, double shape, double scale,
                            double window, int model) {
  if (p <= 0.0) return std::numeric_limits<double>::infinity();
  if (model == ET_TIME_UNIFORM) return u * (window / p);
  return scale * std::pow(-std::log(1.0 - u), 1.0 / shape);
}

// A patient of a time-to-event trial.
struct EtPatient {
  int dose;          // 0-based
  double enter;      // arrival time
  double t_tox;      // time to DLT from arrival (may exceed the window)
  double t_resp;     // time to response from arrival
  double end_tox;    // end of the DLT assessment
  double end_resp;   // end of the response assessment
  char dlt;
  char resp;
  char backfill;
};

// Draw the outcomes of a patient arriving at time t at dose d with two
// uniform variates, the first for the DLT and the second for the response.
template <class Rng>
inline EtPatient et_new_patient(int d, double t, const std::vector<double>& p_true,
                                const std::vector<double>& p_resp,
                                double resp_cor, const EtTiming& tm, Rng& rng) {
  const double u = rng();
  const double v = rng();
  const double w = response_variate(u, v, resp_cor);
  EtPatient pt;
  pt.dose = d;
  pt.enter = t;
  pt.t_tox = et_event_time(u, p_true[d], tm.shape_t[d], tm.scale_t[d],
                           tm.window, tm.event_time);
  pt.t_resp = et_event_time(w, p_resp[d], tm.shape_e[d], tm.scale_e[d],
                            tm.resp_window, tm.event_time);
  pt.dlt = (pt.t_tox <= tm.window) ? 1 : 0;
  pt.resp = (pt.t_resp <= tm.resp_window) ? 1 : 0;
  pt.end_tox = t + std::min(pt.t_tox, tm.window);
  pt.end_resp = t + std::min(pt.t_resp, tm.resp_window);
  pt.backfill = 0;
  return pt;
}

// Rates of every dose with patients at time t. A pending patient counts by the
// fraction of the window followed so far. yy and xx receive the numbers of
// DLTs and responses observed by time t, and pt and pe the rates; doses
// without patients keep their values.
inline void et_tite_rates(const std::vector<EtPatient>& pts, double t,
                          const EtTiming& tm, const std::vector<int>& nn,
                          std::vector<int>& yy, std::vector<int>& xx,
                          std::vector<double>& pt, std::vector<double>& pe) {
  const int n_doses = static_cast<int>(nn.size());
  for (int j = 0; j < n_doses; ++j) {
    if (nn[j] == 0) continue;
    int x_t = 0, c_t = 0, x_e = 0, c_e = 0;
    double s_t = 0.0, s_e = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
      const EtPatient& p = pts[i];
      if (p.dose != j) continue;
      if (p.end_tox <= t) {
        if (p.dlt) ++x_t; else ++c_t;
      } else {
        s_t += t - p.enter;
      }
      if (p.end_resp <= t) {
        if (p.resp) ++x_e; else ++c_e;
      } else {
        s_e += t - p.enter;
      }
    }
    const double n_t = static_cast<double>(x_t) + static_cast<double>(c_t) +
      s_t / tm.window;
    const double n_e = static_cast<double>(x_e) + static_cast<double>(c_e) +
      s_e / tm.resp_window;
    yy[j] = x_t;
    xx[j] = x_e;
    pt[j] = static_cast<double>(x_t) / n_t;
    pe[j] = static_cast<double>(x_e) / n_e;
  }
}

// Decision of BOIN-ET from the rates of every dose, as in et_next_dose() but
// with rates that need not be ratios of counts. nn tells which doses have
// been tried; ties under rule 4 draw one uniform variate.
template <class Rng>
inline int et_next_dose_rates(int cur, const std::vector<int>& nn,
                              const std::vector<double>& pt,
                              const std::vector<double>& pe, double lambda1,
                              double lambda2, double eta1, Rng& rng) {
  const int n_doses = static_cast<int>(nn.size());
  const double t = pt[cur - 1];
  const double e = pe[cur - 1];
  if (t <= lambda1 && e <= eta1) return cur + 1;
  if (t < lambda2 && e > eta1) return cur;
  if (t >= lambda2) return cur - 1;

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
  double best = -1.0;
  for (std::size_t k = 0; k < three.size(); ++k) {
    const double r = pe[three[k] - 1];
    if (r > best) best = r;
  }
  std::vector<int> ties;
  for (std::size_t k = 0; k < three.size(); ++k) {
    if (pe[three[k] - 1] == best) ties.push_back(three[k]);
  }
  if (ties.size() == 1) return ties[0];
  const int k = static_cast<int>(ties.size());
  int pick = static_cast<int>(rng() * k);
  if (pick >= k) pick = k - 1;
  return ties[pick];
}

// Simulate a single TITE-BOIN-ET trial.
//
// des supplies the probabilities, the cohorts, the stopping sample size, the
// admissibility tables and compat; its decision table is not used. nn
// receives the patients treated at every dose, and yy and xx the DLTs and
// responses observed at the last decision. duration is the time of the last
// decision.
//
// For every cohort, rng() is called for the times between the arrivals of its
// patients, then for the time to the next decision (except after the last
// cohort), then twice for each patient in arrival order, and once more for a
// tie under rule 4.
template <class Rng>
inline void simulate_tite_et_one(const EtDesign& des, const EtTiming& tm,
                                 double lambda1, double lambda2, double eta1,
                                 Rng& rng, std::vector<int>& nn,
                                 std::vector<int>& yy, std::vector<int>& xx,
                                 int& cohorts_used, int& stop_code,
                                 double& duration) {
  const int n_doses = static_cast<int>(des.p_true.size());
  const int n_cohorts = static_cast<int>(des.cohort_size.size());
  std::fill(nn.begin(), nn.end(), 0);
  std::fill(yy.begin(), yy.end(), 0);
  std::fill(xx.begin(), xx.end(), 0);
  std::vector<double> pt(n_doses, 0.0), pe(n_doses, 0.0);
  std::vector<char> adm(n_doses);
  std::vector<EtPatient> pts;

  double t_dec = 0.0;
  int cur = des.start_dose;
  cohorts_used = 0;
  stop_code = ET_STOP_MAX_COHORTS;

  for (int i = 0; i < n_cohorts; ++i) {
    ++cohorts_used;
    const int m = des.cohort_size[i];
    const int d = cur - 1;

    std::vector<double> enter(m);
    enter[0] = t_dec;
    for (int k = 1; k < m; ++k) enter[k] = enter[k - 1] + et_gap(tm, rng);
    if (i == n_cohorts - 1) {
      t_dec = enter[m - 1] + std::max(tm.window, tm.resp_window);
    } else {
      t_dec = enter[m - 1] + et_gap(tm, rng);
    }
    for (int k = 0; k < m; ++k) {
      pts.push_back(et_new_patient(d, enter[k], des.p_true, des.p_resp,
                                   des.resp_cor, tm, rng));
    }
    nn[d] += m;

    // Wait until at least half of the patients at the current dose have
    // finished both assessments.
    for (;;) {
      int done = 0;
      double next = std::numeric_limits<double>::infinity();
      for (std::size_t q = 0; q < pts.size(); ++q) {
        const EtPatient& p = pts[q];
        if (p.dose != d) continue;
        if (p.end_tox <= t_dec && p.end_resp <= t_dec) ++done;
        if (p.end_tox > t_dec && p.end_tox < next) next = p.end_tox;
        if (p.end_resp > t_dec && p.end_resp < next) next = p.end_resp;
      }
      if (2 * done >= nn[d]) break;
      t_dec = next;
    }

    et_tite_rates(pts, t_dec, tm, nn, yy, xx, pt, pe);
    const int nxt = et_next_dose_rates(cur, nn, pt, pe, lambda1, lambda2,
                                       eta1, rng);

    const bool any_adm = et_admissible_set(nn, yy, xx, des.tox_ok, des.resp_ok,
                                           des.max_n, des.compat != 0, adm);
    if (!any_adm) {
      stop_code = ET_STOP_NO_ADMISSIBLE;
      break;
    }
    bool reached = false;
    for (int j = 0; j < n_doses; ++j) {
      if (nn[j] >= des.n_stop) reached = true;
    }
    if (reached) {
      stop_code = ET_STOP_N_EARLYSTOP;
      break;
    }
    if (!et_move(cur, nxt, adm, des.compat != 0, stop_code)) break;
  }
  duration = t_dec;
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_TITE_BOINET_CORE_H
