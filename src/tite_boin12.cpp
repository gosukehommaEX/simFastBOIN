// Rcpp glue for the core algorithm of the TITE-BOIN12 design. The algorithm
// itself lives in tite_boin12_core.h so that it can be compiled and verified
// outside of R as well.

#include <Rcpp.h>
#include "tite_boin12_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnifTb12 {
  double operator()() { return unif_rand(); }
};

// Lower or upper tail of the beta distribution.
struct RPbetaTb12 {
  double operator()(double q, double a, double b, bool lower) const {
    return R::pbeta(q, a, b, lower ? 1 : 0, 0);
  }
};

simfastboin::Tb12Design make_tb12_design(
    NumericVector p_true, NumericVector p_resp, NumericVector q0,
    NumericVector q1, int cor_type, double resp_cor, IntegerVector cohort_size,
    int start_dose, int n_earlystop, bool early_stop_simple, int n_star,
    int n_explore, int tie_break, int no_candidate, int explore_from,
    double lambda_e, double lambda_d, double psi00, double psi11,
    double alpha, double beta, double u_benchmark, double target,
    double target_resp, double cutoff_eli, double cutoff_eli_resp,
    int pending) {
  simfastboin::Tb12Design des;
  simfastboin::B12Design& b = des.b12;
  b.p_true.assign(p_true.begin(), p_true.end());
  b.p_resp.assign(p_resp.begin(), p_resp.end());
  b.q0.assign(q0.begin(), q0.end());
  b.q1.assign(q1.begin(), q1.end());
  b.cor_type = cor_type;
  b.resp_cor = resp_cor;
  b.cohort_size.assign(cohort_size.begin(), cohort_size.end());
  b.start_dose = start_dose;
  b.n_earlystop = n_earlystop;
  b.early_stop_simple = early_stop_simple ? 1 : 0;
  b.n_star = n_star;
  b.n_explore = n_explore;
  b.tie_random = tie_break == 1 ? 1 : 0;
  b.stop_no_candidate = no_candidate == 1 ? 1 : 0;
  b.explore_admissible_only = explore_from == 1 ? 1 : 0;
  b.lambda_e = lambda_e;
  b.lambda_d = lambda_d;
  b.psi00 = psi00;
  b.psi11 = psi11;
  b.alpha = alpha;
  b.beta = beta;
  b.u_benchmark = u_benchmark;
  b.max_n = 0;
  des.target = target;
  des.target_resp = target_resp;
  des.cutoff_eli = cutoff_eli;
  des.cutoff_eli_resp = cutoff_eli_resp;
  des.pending = pending;
  return des;
}

}  // namespace

// [[Rcpp::export]]
List tite_boin12_simulate_cpp(int n_trials, NumericVector p_true,
                              NumericVector p_resp, NumericVector q0,
                              NumericVector q1, int cor_type,
                              double resp_cor, IntegerVector cohort_size,
                              int start_dose, int n_earlystop,
                              bool early_stop_simple, int n_star,
                              int n_explore, int tie_break,
                              int no_candidate, int explore_from,
                              double lambda_e, double lambda_d, double psi00,
                              double psi11, double alpha, double beta,
                              double u_benchmark,
                              double target, double target_resp,
                              double cutoff_eli, double cutoff_eli_resp,
                              int pending, double window,
                              double resp_window, double accrual_rate,
                              int arrival, int event_time,
                              NumericVector shape_t, NumericVector scale_t,
                              NumericVector shape_e, NumericVector scale_e,
                              double min_completed_share) {

  const simfastboin::Tb12Design des = make_tb12_design(
    p_true, p_resp, q0, q1, cor_type, resp_cor, cohort_size, start_dose,
    n_earlystop, early_stop_simple, n_star, n_explore, tie_break,
    no_candidate, explore_from, lambda_e, lambda_d, psi00, psi11, alpha, beta,
    u_benchmark, target, target_resp, cutoff_eli, cutoff_eli_resp, pending);

  simfastboin::EtTiming tm;
  tm.window = window;
  tm.resp_window = resp_window;
  tm.accrual_rate = accrual_rate;
  tm.arrival = arrival;
  tm.event_time = event_time;
  tm.shape_t.assign(shape_t.begin(), shape_t.end());
  tm.scale_t.assign(scale_t.begin(), scale_t.end());
  tm.shape_e.assign(shape_e.begin(), shape_e.end());
  tm.scale_e.assign(scale_e.begin(), scale_e.end());
  tm.min_share = min_completed_share;

  const int n_doses = p_true.size();
  IntegerMatrix n_pts(n_trials, n_doses);
  IntegerMatrix n_tox(n_trials, n_doses);
  IntegerMatrix n_resp(n_trials, n_doses);
  IntegerMatrix n_tox_resp(n_trials, n_doses);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);
  NumericVector duration(n_trials);

  std::vector<int> nn(n_doses), yy(n_doses), xx(n_doses), zz(n_doses);
  std::vector<simfastboin::Tb12Patient> pts;
  RUnifTb12 rng;
  RPbetaTb12 pb;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    double dur = 0.0;
    simfastboin::tb12_simulate_one(des, tm, rng, pb, pts, nn, yy, xx, zz,
                                   used, code, dur);
    for (int j = 0; j < n_doses; ++j) {
      n_pts(t, j) = nn[j];
      n_tox(t, j) = yy[j];
      n_resp(t, j) = xx[j];
      n_tox_resp(t, j) = zz[j];
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
    duration[t] = dur;
  }

  return List::create(
    _["n_pts"] = n_pts,
    _["n_tox"] = n_tox,
    _["n_resp"] = n_resp,
    _["n_tox_resp"] = n_tox_resp,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code,
    _["duration"] = duration
  );
}

// Interim quantities of every dose at time `time` from patient data, as used
// by the decisions of tite_boin12_simulate(). A patient is described by the
// dose (1-based), the arrival time, the times at which the DLT and the
// response outcomes are observed (the time of the event, or the end of its
// window) and the two outcomes. Patients arriving after `time` are ignored.
// Not exported; used by the tests and by inst/reproduce/reproduce_Zhou2022.R.
// [[Rcpp::export]]
List tite_boin12_interim_cpp(IntegerVector dose, NumericVector enter,
                             NumericVector end_tox, NumericVector end_resp,
                             IntegerVector dlt, IntegerVector resp,
                             double time, int n_doses, double window,
                             double resp_window, double psi00, double psi11,
                             double alpha, double beta, double u_benchmark,
                             double target, double target_resp,
                             double cutoff_eli, double cutoff_eli_resp,
                             int pending) {

  NumericVector empty(n_doses);
  IntegerVector one(1, 1);
  const simfastboin::Tb12Design des = make_tb12_design(
    empty, empty, empty, empty, 0, 0.0, one, 1, 0, false, 6, 9, 0, 0, 0,
    0.0, 1.0, psi00, psi11, alpha, beta, u_benchmark, target, target_resp,
    cutoff_eli, cutoff_eli_resp, pending);

  simfastboin::EtTiming tm;
  tm.window = window;
  tm.resp_window = resp_window;

  std::vector<simfastboin::Tb12Patient> pts;
  for (int i = 0; i < dose.size(); ++i) {
    if (enter[i] > time) continue;
    simfastboin::Tb12Patient p;
    p.dose = dose[i] - 1;
    p.enter = enter[i];
    p.end_tox = end_tox[i];
    p.end_resp = end_resp[i];
    p.dlt = dlt[i] != 0 ? 1 : 0;
    p.resp = resp[i] != 0 ? 1 : 0;
    pts.push_back(p);
  }

  std::vector<simfastboin::Tb12Interim> est(n_doses);
  for (int j = 0; j < n_doses; ++j) {
    simfastboin::tb12_interim(pts, j, time, tm, des.b12, pending, est[j]);
  }
  RPbetaTb12 pb;
  std::vector<char> adm(n_doses);
  simfastboin::tb12_admissible_set(est, des, pb, adm);

  IntegerVector n(n_doses);
  NumericVector pi_tox(n_doses), pi_resp(n_doses), quasi(n_doses),
    desir(n_doses);
  LogicalVector admissible(n_doses);
  for (int j = 0; j < n_doses; ++j) {
    n[j] = est[j].n;
    pi_tox[j] = est[j].pi_t;
    pi_resp[j] = est[j].pi_e;
    quasi[j] = est[j].q;
    desir[j] = simfastboin::tb12_desirability(est[j], des.b12, pb);
    admissible[j] = adm[j] != 0;
  }

  return List::create(
    _["n"] = n,
    _["pi_tox"] = pi_tox,
    _["pi_resp"] = pi_resp,
    _["quasi_events"] = quasi,
    _["desirability"] = desir,
    _["admissible"] = admissible
  );
}
