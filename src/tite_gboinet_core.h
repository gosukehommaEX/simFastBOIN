#ifndef SIMFASTBOIN_TITE_GBOINET_CORE_H
#define SIMFASTBOIN_TITE_GBOINET_CORE_H

// Core algorithm of the TITE-gBOIN-ET design (Takeda, Yamaguchi, Taguri and
// Morita, 2023), the time-to-event version of gBOIN-ET for toxicity and
// efficacy outcomes with several categories that each carry a score,
// written in plain C++ like gboinet_core.h and tite_boinet_core.h so that it
// can be verified outside of R as well.
//
// The rules follow tite.gboinet() of the boinet package (version 1.6.0),
// with the two changes of boinet_core.h unless compat is set:
//   - the patients arrive, the decisions are due and the decisions wait as
//     in simulate_tite_et_one();
//   - the mean score of a dose is the total of the scores observed so far
//     divided by the effective number of patients, the patients whose
//     assessment has ended plus, for every patient still under assessment,
//     the fraction of the window followed so far, and the decision rule of
//     gBOIN-ET is applied to these means;
//   - the admissibility of a dose uses the total of the scores observed so
//     far and the number of patients treated, so that a patient still under
//     assessment counts as a patient without an event, and the trial moves
//     between doses as in simulate_get_one().
// A patient has an event when the event time falls within the window; the
// category of a patient with an event is drawn as in simulate_get_one() and
// is observed at the time of the event.
//
// Floating point. The mean score of a dose is one division of the total of
// the scores as given (raw_score_total() of gboin_core.h) by the largest
// score times the effective number of patients. When no patient of the dose
// is still under assessment the effective number is the number of patients
// and the mean is that of get_mean() of gboinet_core.h, so that doses whose
// mean scores are equal as fractions compare equal on every platform.
//
// The only external dependency is beta_cdf(), which the caller supplies.

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "tite_boinet_core.h"
#include "gboinet_core.h"

namespace simfastboin {

// A patient of a time-to-event trial with categories.
struct TgetPatient {
  int dose;          // 0-based
  double enter;      // arrival time
  double end_tox;    // end of the toxicity assessment
  double end_resp;   // end of the efficacy assessment
  int cat_tox;       // toxicity category, 0 without an event
  int cat_resp;      // efficacy category, 0 without an event
};

// Mean normalized score from the total raw_total of the scores as given and
// the effective number of patients n_eff: one division, by the largest
// score times the effective number.
inline double tget_mean(double raw_total, double largest, double n_eff) {
  volatile double denom = largest * n_eff;
  return raw_total / denom;
}

// Draw the outcomes of a patient arriving at time t at dose d (0-based).
// The variates are drawn in the order of simulate_get_one(): the toxicity
// event, the efficacy event, the toxicity category when the toxicity
// outcome has more than two categories and the efficacy category when the
// efficacy outcome has more than two.
template <class Rng>
inline TgetPatient tget_new_patient(int d, double t, const GetDesign& des,
                                    const EtTiming& tm, Rng& rng) {
  const int n_cat_t = static_cast<int>(des.score_tox.size());
  const int n_cat_e = static_cast<int>(des.score_resp.size());
  const double u = rng();
  const double v = rng();
  const double g = (n_cat_t > 2) ? rng() : 0.0;
  const double h = (n_cat_e > 2) ? rng() : 0.0;
  const double w = response_variate(u, v, des.resp_cor);
  const double* cum_t = &des.cum_tox[d * n_cat_t];
  const double* cum_e = &des.cum_resp[d * n_cat_e];
  const double t_tox = et_event_time(u, cum_t[1], tm.shape_t[d], tm.scale_t[d],
                                     tm.window, tm.event_time);
  const double t_resp = et_event_time(w, cum_e[1], tm.shape_e[d],
                                      tm.scale_e[d], tm.resp_window,
                                      tm.event_time);
  TgetPatient pt;
  pt.dose = d;
  pt.enter = t;
  // get_category() takes an event when its first variate is below the
  // probability of an event.
  pt.cat_tox = (t_tox <= tm.window) ? get_category(cum_t, n_cat_t, -1.0, g)
                                    : 0;
  pt.cat_resp = (t_resp <= tm.resp_window)
    ? get_category(cum_e, n_cat_e, -1.0, h) : 0;
  pt.end_tox = t + std::min(t_tox, tm.window);
  pt.end_resp = t + std::min(t_resp, tm.resp_window);
  return pt;
}

// Data of every dose with patients at time t. count_tox and count_resp
// receive the patients by category (index j * n_cat + k) with the
// categories observed by time t, a patient still under assessment being
// counted in the first category; rt and re the totals of the scores as
// given; pt and pe the mean normalized scores, in which a patient still
// under assessment counts by the fraction of the window followed so far.
// Doses without patients keep their values.
inline void tget_rates(const std::vector<TgetPatient>& pts, double t,
                       const EtTiming& tm, const GetDesign& des,
                       const std::vector<int>& nn,
                       std::vector<int>& count_tox,
                       std::vector<int>& count_resp,
                       std::vector<double>& rt, std::vector<double>& re,
                       std::vector<double>& pt, std::vector<double>& pe) {
  const int n_doses = static_cast<int>(nn.size());
  const int n_cat_t = static_cast<int>(des.score_tox.size());
  const int n_cat_e = static_cast<int>(des.score_resp.size());
  for (int j = 0; j < n_doses; ++j) {
    if (nn[j] == 0) continue;
    int* ct = &count_tox[j * n_cat_t];
    int* ce = &count_resp[j * n_cat_e];
    std::fill(ct, ct + n_cat_t, 0);
    std::fill(ce, ce + n_cat_e, 0);
    int done_t = 0, done_e = 0;
    double s_t = 0.0, s_e = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
      const TgetPatient& p = pts[i];
      if (p.dose != j) continue;
      if (p.end_tox <= t) {
        ++ct[p.cat_tox];
        ++done_t;
      } else {
        ++ct[0];
        s_t += t - p.enter;
      }
      if (p.end_resp <= t) {
        ++ce[p.cat_resp];
        ++done_e;
      } else {
        ++ce[0];
        s_e += t - p.enter;
      }
    }
    const double n_t = static_cast<double>(done_t) + s_t / tm.window;
    const double n_e = static_cast<double>(done_e) + s_e / tm.resp_window;
    rt[j] = raw_score_total(ct, des.score_tox);
    re[j] = raw_score_total(ce, des.score_resp);
    pt[j] = tget_mean(rt[j], des.max_tox, n_t);
    pe[j] = tget_mean(re[j], des.max_resp, n_e);
  }
}

// Simulate a single TITE-gBOIN-ET trial.
//
// des supplies the categories, the scores, the cohorts, the stopping sample
// size, the boundaries, the elimination rule and compat. nn receives the
// patients treated at every dose, and count_tox and count_resp the patients
// by dose and category (index j * n_cat + k) as observed at the last
// decision, a patient still under assessment being counted in the first
// category. duration is the time of the last decision.
//
// For every cohort, rng() is called for the times between the arrivals of
// its patients, then for the time to the next decision (except after the
// last cohort), then for each patient in arrival order as in
// tget_new_patient(), and once more for a tie under rule 4.
template <class Rng>
inline void simulate_tget_one(const GetDesign& des, const EtTiming& tm,
                              Rng& rng, std::vector<int>& nn,
                              std::vector<int>& count_tox,
                              std::vector<int>& count_resp,
                              int& cohorts_used, int& stop_code,
                              double& duration) {
  const int n_cat_t = static_cast<int>(des.score_tox.size());
  const int n_doses = static_cast<int>(des.cum_tox.size()) / n_cat_t;
  const int n_cohorts = static_cast<int>(des.cohort_size.size());
  std::fill(nn.begin(), nn.end(), 0);
  std::fill(count_tox.begin(), count_tox.end(), 0);
  std::fill(count_resp.begin(), count_resp.end(), 0);
  std::vector<double> rt(n_doses, 0.0), re(n_doses, 0.0);
  std::vector<double> pt(n_doses, 0.0), pe(n_doses, 0.0);
  std::vector<char> t_ok(n_doses), e_ok(n_doses), adm(n_doses);
  const char t0 = get_tox_ok(0.0, 0, des.p_eli, des.cut_tox) ? 1 : 0;
  const char e0 = get_resp_ok(0.0, 0, des.p_eli_resp, des.cut_resp) ? 1 : 0;
  std::fill(t_ok.begin(), t_ok.end(), t0);
  std::fill(e_ok.begin(), e_ok.end(), e0);
  std::vector<TgetPatient> pts;

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
      pts.push_back(tget_new_patient(d, enter[k], des, tm, rng));
    }
    nn[d] += m;

    // Wait until at least a share min_share of the patients at the current
    // dose have finished both assessments.
    for (;;) {
      int done = 0;
      double next = std::numeric_limits<double>::infinity();
      for (std::size_t q = 0; q < pts.size(); ++q) {
        const TgetPatient& p = pts[q];
        if (p.dose != d) continue;
        if (p.end_tox <= t_dec && p.end_resp <= t_dec) ++done;
        if (p.end_tox > t_dec && p.end_tox < next) next = p.end_tox;
        if (p.end_resp > t_dec && p.end_resp < next) next = p.end_resp;
      }
      if (static_cast<double>(done) / nn[d] >= tm.min_share) break;
      t_dec = next;
    }

    tget_rates(pts, t_dec, tm, des, nn, count_tox, count_resp, rt, re, pt,
               pe);
    const int nxt = et_next_dose_rates(cur, nn, pt, pe, des.lambda1,
                                       des.lambda2, des.eta1, rng);

    for (int j = 0; j < n_doses; ++j) {
      if (nn[j] == 0) continue;
      t_ok[j] = get_tox_ok(rt[j] / des.max_tox, nn[j], des.p_eli,
                           des.cut_tox) ? 1 : 0;
      e_ok[j] = get_resp_ok(re[j] / des.max_resp, nn[j], des.p_eli_resp,
                            des.cut_resp) ? 1 : 0;
    }
    const bool any_adm = get_admissible_set(t_ok, e_ok, des.compat != 0, adm);
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

#endif  // SIMFASTBOIN_TITE_GBOINET_CORE_H
