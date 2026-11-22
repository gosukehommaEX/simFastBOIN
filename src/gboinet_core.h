#ifndef SIMFASTBOIN_GBOINET_CORE_H
#define SIMFASTBOIN_GBOINET_CORE_H

// Core algorithms of the gBOIN-ET design (Takeda, Morita and Taguri, 2022),
// in which the toxicity and the efficacy outcomes have several categories
// that each carry a score between 0 and 1 (quasi-Bernoulli outcomes), written
// in plain C++ so that they can be compiled inside the package and inside a
// standalone verification harness.
//
// The rules are those of the BOIN-ET core (boinet_core.h) with the total
// scores in the place of the numbers of DLTs and of responses. The totals are
// not whole numbers, so the decisions and the posterior probabilities are
// evaluated here and not read from tables. With two categories of scores 0
// and 1 for both outcomes the trials are those of simulate_et_one().
//
// The only external dependency is beta_cdf(), which the caller supplies.

#include <vector>
#include <cmath>
#include <algorithm>
#include "boinet_core.h"
#include "gboin_core.h"

namespace simfastboin {

// Distribution function of the Beta(shape1, shape2) distribution at x.
// Defined by the caller (in the package, in comb.cpp).
double beta_cdf(double x, double shape1, double shape2);

// The settings of a gBOIN-ET trial that do not change from trial to trial.
struct GetDesign {
  // Probability of category k or higher at dose j, at index j * n_cat + k,
  // the categories in the order of the scores given; the first category is
  // the one without an event.
  std::vector<double> cum_tox;
  std::vector<double> cum_resp;
  std::vector<double> score_tox;   // normalized scores, one per category
  std::vector<double> score_resp;
  double resp_cor;                 // copula correlation of the two events
  std::vector<int> cohort_size;    // one entry per cohort
  int start_dose;                  // 1-based
  int n_stop;                      // stopping sample size at any dose
  double lambda1;
  double lambda2;
  double eta1;
  double p_eli;                    // toxicity score of the elimination rule
  double p_eli_resp;               // efficacy score of the elimination rule
  double cut_tox;                  // 1 - cutoff_eli
  double cut_resp;                 // 1 - cutoff_eli_resp
  int compat;                      // 1 for the rules of boinet 1.6.0
};

// Whether the toxicity data of a dose, a total score x among n patients, do
// not exclude it: Pr(mean score < p_eli | data) is at least cut under a
// uniform prior.
inline bool get_tox_ok(double x, int n, double p_eli, double cut) {
  return !(beta_cdf(p_eli, 1.0 + x, 1.0 + (n - x)) < cut);
}

// Whether the efficacy data of a dose do not exclude it:
// Pr(mean score > p_eli_resp | data) is at least cut.
inline bool get_resp_ok(double x, int n, double p_eli_resp, double cut) {
  return !((1.0 - beta_cdf(p_eli_resp, 1.0 + x, 1.0 + (n - x))) < cut);
}

// Admissibility of every dose from the flags of its toxicity and efficacy
// data. Unless compat is set, a dose excluded by its toxicity data also
// excludes every higher dose. Returns whether any dose is admissible.
inline bool get_admissible_set(const std::vector<char>& t_ok,
                               const std::vector<char>& e_ok, bool compat,
                               std::vector<char>& adm) {
  const int n_doses = static_cast<int>(t_ok.size());
  bool tox_out = false;
  bool any_adm = false;
  for (int j = 0; j < n_doses; ++j) {
    if (!t_ok[j] && !compat) tox_out = true;
    const bool ok = t_ok[j] && !tox_out && e_ok[j];
    adm[j] = ok ? 1 : 0;
    if (ok) any_adm = true;
  }
  return any_adm;
}

// Category of a patient at a dose from the variate of the event and the
// variate of the category. cum points at the probabilities of category k or
// higher at that dose. The patient has an event when event_variate is below
// cum[1]; given an event, the category is the highest k with
// cat_variate * cum[1] below cum[k].
inline int get_category(const double* cum, int n_cat, double event_variate,
                        double cat_variate) {
  if (!(event_variate < cum[1])) return 0;
  const double g = cat_variate * cum[1];
  for (int c = n_cat - 1; c >= 2; --c) {
    if (g < cum[c]) return c;
  }
  return 1;
}

// Next dose (1-based, from 0 to n_doses + 1) chosen by the decision rule
// from the total scores xt and xe of the patients nn.
//
// Under rule 4 the dose with the highest observed mean efficacy score among
// the current dose and its neighbors is chosen; a dose not yet tried counts
// as zero. Ties are broken with one uniform variate.
template <class Rng>
inline int get_next_dose(int cur, const std::vector<int>& nn,
                         const std::vector<double>& xt,
                         const std::vector<double>& xe,
                         const GetDesign& des, Rng& rng) {
  const int n_doses = static_cast<int>(nn.size());
  const int c = cur - 1;
  const double pt = xt[c] / static_cast<double>(nn[c]);
  const double pe_cur = xe[c] / static_cast<double>(nn[c]);
  if (pt <= des.lambda1 && pe_cur <= des.eta1) return cur + 1;
  if (pt < des.lambda2 && pe_cur > des.eta1) return cur;
  if (pt >= des.lambda2) return cur - 1;

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
    pe[k] = (nn[d] > 0) ? xe[d] / static_cast<double>(nn[d]) : 0.0;
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

// Simulate a single gBOIN-ET trial.
//
// rng() must return one uniform variate on (0, 1) per call. Each patient
// draws, in enrollment order, one variate for the toxicity event, one for
// the efficacy event, then one for the toxicity category when the toxicity
// outcome has more than two categories and one for the efficacy category
// when the efficacy outcome has more than two; the category variates are
// drawn whether or not the event occurs. The two events are linked by a
// Gaussian copula as in simulate_et_one(), and the category of a patient
// with an event is independent of everything else. A tie under rule 4 draws
// one more variate.
//
// count_tox and count_resp receive the numbers of patients by dose and
// category (index j * n_cat + k), nn the patients by dose.
template <class Rng>
inline void simulate_get_one(const GetDesign& des, Rng& rng,
                             std::vector<int>& nn,
                             std::vector<int>& count_tox,
                             std::vector<int>& count_resp,
                             int& cohorts_used, int& stop_code) {
  const int n_cat_t = static_cast<int>(des.score_tox.size());
  const int n_cat_e = static_cast<int>(des.score_resp.size());
  const int n_doses = static_cast<int>(des.cum_tox.size()) / n_cat_t;
  std::fill(nn.begin(), nn.end(), 0);
  std::fill(count_tox.begin(), count_tox.end(), 0);
  std::fill(count_resp.begin(), count_resp.end(), 0);

  std::vector<double> xt(n_doses, 0.0), xe(n_doses, 0.0);
  std::vector<char> t_ok(n_doses), e_ok(n_doses), adm(n_doses);
  const char t0 = get_tox_ok(0.0, 0, des.p_eli, des.cut_tox) ? 1 : 0;
  const char e0 = get_resp_ok(0.0, 0, des.p_eli_resp, des.cut_resp) ? 1 : 0;
  std::fill(t_ok.begin(), t_ok.end(), t0);
  std::fill(e_ok.begin(), e_ok.end(), e0);

  int cur = des.start_dose;
  cohorts_used = 0;
  stop_code = ET_STOP_MAX_COHORTS;

  for (std::size_t i = 0; i < des.cohort_size.size(); ++i) {
    ++cohorts_used;
    const int d = cur - 1;
    for (int k = 0; k < des.cohort_size[i]; ++k) {
      const double u = rng();
      const double v = rng();
      const double g = (n_cat_t > 2) ? rng() : 0.0;
      const double h = (n_cat_e > 2) ? rng() : 0.0;
      const double w = response_variate(u, v, des.resp_cor);
      const int ct = get_category(&des.cum_tox[d * n_cat_t], n_cat_t, u, g);
      const int ce = get_category(&des.cum_resp[d * n_cat_e], n_cat_e, w, h);
      nn[d] += 1;
      ++count_tox[d * n_cat_t + ct];
      ++count_resp[d * n_cat_e + ce];
    }
    xt[d] = score_total(&count_tox[d * n_cat_t], des.score_tox);
    xe[d] = score_total(&count_resp[d * n_cat_e], des.score_resp);
    t_ok[d] = get_tox_ok(xt[d], nn[d], des.p_eli, des.cut_tox) ? 1 : 0;
    e_ok[d] = get_resp_ok(xe[d], nn[d], des.p_eli_resp, des.cut_resp) ? 1 : 0;

    const int nxt = get_next_dose(cur, nn, xt, xe, des, rng);

    const bool any_adm = get_admissible_set(t_ok, e_ok, des.compat != 0, adm);
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

// Weighted least squares fit of y on the columns of x (n_col columns of
// length k, column c at x[c * k + i]) with weights w. The columns are
// orthogonalized one after the other, and a column that the earlier ones
// already span is dropped, so that the fitted values are the projection of
// y whatever the rank. Returns the weighted residual sum of squares and
// writes the fitted values.
inline double weighted_ls_fit(const std::vector<double>& x, int n_col, int k,
                              const std::vector<double>& y,
                              const std::vector<double>& w,
                              std::vector<double>& fit) {
  std::vector<double> q;
  std::vector<double> qq;
  q.reserve(static_cast<std::size_t>(n_col) * k);
  std::fill(fit.begin(), fit.end(), 0.0);
  int n_q = 0;
  for (int c = 0; c < n_col; ++c) {
    std::vector<double> v(x.begin() + c * k, x.begin() + (c + 1) * k);
    double norm0 = 0.0;
    for (int i = 0; i < k; ++i) norm0 += w[i] * v[i] * v[i];
    for (int pass = 0; pass < 2; ++pass) {
      for (int b = 0; b < n_q; ++b) {
        double dot = 0.0;
        for (int i = 0; i < k; ++i) dot += w[i] * v[i] * q[b * k + i];
        const double coef = dot / qq[b];
        for (int i = 0; i < k; ++i) v[i] -= coef * q[b * k + i];
      }
    }
    double norm = 0.0;
    for (int i = 0; i < k; ++i) norm += w[i] * v[i] * v[i];
    if (!(norm > 1e-10 * norm0)) continue;
    q.insert(q.end(), v.begin(), v.end());
    qq.push_back(norm);
    ++n_q;
    double dot = 0.0;
    for (int i = 0; i < k; ++i) dot += w[i] * y[i] * v[i];
    const double coef = dot / norm;
    for (int i = 0; i < k; ++i) fit[i] += coef * v[i];
  }
  double sse = 0.0;
  for (int i = 0; i < k; ++i) sse += w[i] * (y[i] - fit[i]) * (y[i] - fit[i]);
  return sse;
}

// Column of a fractional polynomial: dose^p, or log(dose) for p = 0.
inline double fp_power(double dose, double p) {
  return (p == 0.0) ? std::log(dose) : std::pow(dose, p);
}

// Mean efficacy scores estimated by the linear fractional polynomial of
// degree 2 of Takeda et al. (2022): mean = b0 + b1 dose^p1 + b2 dose^p2 with
// the powers among -2, -1, -0.5, 0 (the logarithm), 0.5, 1, 2 and 3, a
// repeated power p giving dose^p and dose^p log(dose). The 36 models are
// fitted by least squares to the scores of the patients, which is the
// weighted fit to the mean scores ybar of the doses with their numbers of
// patients w as weights, and the model with the smallest residual sum of
// squares is taken, the first of the list when several are within 1e-12 of
// the smallest. dose holds the dose levels (1-based) of the doses tried.
// With three doses or fewer every model passes through the observed means,
// which are returned. The fitted values are rounded to ten decimals, so
// that equal estimates compare equal.
inline void fp_linear_estimate(const std::vector<double>& dose,
                               const std::vector<double>& ybar,
                               const std::vector<double>& w,
                               std::vector<double>& est) {
  const int k = static_cast<int>(dose.size());
  est.assign(ybar.begin(), ybar.end());
  if (k <= 3) return;
  static const double powers[8] = {-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0};
  std::vector<double> x(3 * static_cast<std::size_t>(k));
  std::vector<double> fit(k), best_fit(k);
  double best = 0.0;
  bool first = true;
  for (int a = 0; a < 8; ++a) {
    for (int b = a; b < 8; ++b) {
      for (int i = 0; i < k; ++i) {
        const double x1 = fp_power(dose[i], powers[a]);
        x[i] = 1.0;
        x[k + i] = x1;
        x[2 * k + i] = (b == a) ? x1 * std::log(dose[i])
                                : fp_power(dose[i], powers[b]);
      }
      const double sse = weighted_ls_fit(x, 3, k, ybar, w, fit);
      if (first || sse < best - 1e-12) {
        best = sse;
        best_fit = fit;
        first = false;
      }
    }
  }
  for (int i = 0; i < k; ++i) {
    est[i] = std::floor(best_fit[i] * 1e10 + 0.5) / 1e10;
  }
}

// Index, among the evaluated doses, of the OBD chosen by rule. cand marks
// the admissible doses, from which the MTD is found, and pick the doses that
// can be selected, a subset of cand. With pick equal to cand this is
// et_obd_choose(). Returns -1 when no dose can be selected.
inline int get_obd_choose(const std::vector<double>& probt,
                          const std::vector<double>& probe,
                          const std::vector<char>& cand,
                          const std::vector<char>& pick,
                          const EtObdRule& rule) {
  const int k_doses = static_cast<int>(probt.size());
  bool any_pick = false;
  bool same = true;
  for (int k = 0; k < k_doses; ++k) {
    if (pick[k]) any_pick = true;
    if (pick[k] != cand[k]) same = false;
  }
  if (!any_pick) return -1;
  if (same) return et_obd_choose(probt, probe, cand, rule);
  if (rule.method != ET_OBD_MAX_RESP) {
    return et_obd_choose(probt, probe, pick, rule);
  }

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
    if (!pick[k]) continue;
    if (first || probe[k] > meff) {
      meff = probe[k];
      first = false;
    }
  }
  if (first) return -1;
  for (int k = 0; k <= mtd; ++k) {
    if (pick[k] && probe[k] == meff) return k;
  }
  return -1;
}

// Efficacy estimators of select_obd_get_one().
enum GetEstimate {
  GET_EST_GIVEN = 0,      // probe_all holds the estimates
  GET_EST_FP_LINEAR = 1   // linear fractional polynomial, computed here
};

// Select the OBD from the final data of a single trial: the patients n and
// the total scores xt and xe of every dose.
//
// With estimate GET_EST_GIVEN, probe_all holds the efficacy estimate of
// every dose and only the doses with at least one patient are read; with
// GET_EST_FP_LINEAR the estimates are those of fp_linear_estimate() and are
// written to probe_all. Only doses with at least min_obd_sample patients
// can be selected; the MTD is found among all admissible doses. Unless
// compat is set, a dose excluded by its toxicity data also excludes every
// higher dose. Returns the OBD (1-based) or 0 when no dose is selected, and
// sets reason to an EtObdReason.
inline int select_obd_get_one(const int* n, const double* xt, const double* xe,
                              double* probe_all, int n_doses, bool early_stop,
                              double p_eli, double cut_tox, double p_eli_resp,
                              double cut_resp, bool compat, int min_obd_sample,
                              int estimate, const EtObdRule& rule,
                              int& reason) {
  if (early_stop) {
    reason = ET_OBD_EARLY_STOP;
    return 0;
  }
  std::vector<int> eva;
  for (int j = 0; j < n_doses; ++j) {
    if (n[j] > 0) eva.push_back(j);
  }
  const int k_doses = static_cast<int>(eva.size());
  std::vector<char> t_ok(n_doses), e_ok(n_doses), adm_all(n_doses);
  for (int j = 0; j < n_doses; ++j) {
    t_ok[j] = get_tox_ok(xt[j], n[j], p_eli, cut_tox) ? 1 : 0;
    e_ok[j] = get_resp_ok(xe[j], n[j], p_eli_resp, cut_resp) ? 1 : 0;
  }
  get_admissible_set(t_ok, e_ok, compat, adm_all);
  std::vector<char> cand(k_doses), pick(k_doses);
  bool any_cand = false;
  for (int k = 0; k < k_doses; ++k) {
    cand[k] = adm_all[eva[k]];
    pick[k] = (cand[k] && n[eva[k]] >= min_obd_sample) ? 1 : 0;
    if (cand[k]) any_cand = true;
  }
  if (!any_cand) {
    reason = ET_OBD_NO_ADMISSIBLE;
    return 0;
  }
  if (k_doses == 1) {
    if (!pick[0]) {
      reason = ET_OBD_NO_ADMISSIBLE;
      return 0;
    }
    reason = ET_OBD_SELECTED;
    return eva[0] + 1;
  }

  std::vector<double> probt(k_doses), probe(k_doses);
  for (int k = 0; k < k_doses; ++k) {
    const int j = eva[k];
    probt[k] = xt[j] / static_cast<double>(n[j]);
    probe[k] = probe_all[j];
  }
  pava_equal(probt);
  if (estimate == GET_EST_FP_LINEAR) {
    std::vector<double> dose(k_doses), ybar(k_doses), w(k_doses);
    for (int k = 0; k < k_doses; ++k) {
      const int j = eva[k];
      dose[k] = j + 1.0;
      ybar[k] = xe[j] / static_cast<double>(n[j]);
      w[k] = n[j];
    }
    fp_linear_estimate(dose, ybar, w, probe);
    for (int k = 0; k < k_doses; ++k) probe_all[eva[k]] = probe[k];
  }

  const int k = get_obd_choose(probt, probe, cand, pick, rule);
  if (k < 0) {
    reason = ET_OBD_NO_ADMISSIBLE;
    return 0;
  }
  reason = ET_OBD_SELECTED;
  return eva[k] + 1;
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_GBOINET_CORE_H
