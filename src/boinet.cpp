// Rcpp glue for the core algorithms of the BOIN-ET family (BOIN-ET,
// TITE-BOIN-ET and BF-BOIN-ET). The algorithms themselves live in
// boinet_core.h, tite_boinet_core.h and bf_boinet_core.h so that they can be
// compiled and verified outside of R as well.

#include <Rcpp.h>
#include "boinet_core.h"
#include "tite_boinet_core.h"
#include "bf_boinet_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnifEt {
  double operator()() { return unif_rand(); }
};

}  // namespace

// [[Rcpp::export]]
List boin_et_simulate_cpp(int n_trials,
                          NumericVector p_true,
                          NumericVector p_resp,
                          NumericVector q0,
                          NumericVector q1,
                          int cor_type,
                          double resp_cor,
                          IntegerVector cohort_size,
                          int start_dose,
                          int n_earlystop,
                          IntegerVector dec,
                          IntegerVector tox_ok,
                          IntegerVector resp_ok,
                          int max_n,
                          int boinet_compatible) {

  simfastboin::EtDesign des;
  des.p_true.assign(p_true.begin(), p_true.end());
  des.p_resp.assign(p_resp.begin(), p_resp.end());
  des.q0.assign(q0.begin(), q0.end());
  des.q1.assign(q1.begin(), q1.end());
  des.cor_type = cor_type;
  des.resp_cor = resp_cor;
  des.cohort_size.assign(cohort_size.begin(), cohort_size.end());
  des.start_dose = start_dose;
  des.n_stop = n_earlystop;
  des.max_n = max_n;
  des.dec.assign(dec.begin(), dec.end());
  des.tox_ok.assign(tox_ok.begin(), tox_ok.end());
  des.resp_ok.assign(resp_ok.begin(), resp_ok.end());
  des.compat = boinet_compatible;

  const int n_doses = p_true.size();

  IntegerMatrix n_pts(n_trials, n_doses);
  IntegerMatrix n_tox(n_trials, n_doses);
  IntegerMatrix n_resp(n_trials, n_doses);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nn(n_doses), yy(n_doses), xx(n_doses);
  RUnifEt rng;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    simfastboin::simulate_et_one(des, rng, nn, yy, xx, used, code);
    for (int j = 0; j < n_doses; ++j) {
      n_pts(t, j) = nn[j];
      n_tox(t, j) = yy[j];
      n_resp(t, j) = xx[j];
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
  }

  return List::create(
    _["n_pts"] = n_pts,
    _["n_tox"] = n_tox,
    _["n_resp"] = n_resp,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code
  );
}

// [[Rcpp::export]]
List boin_et_select_obd_cpp(IntegerMatrix n_pts,
                            IntegerMatrix n_tox,
                            IntegerMatrix n_resp,
                            NumericMatrix probe,
                            LogicalVector early_stop,
                            IntegerVector tox_ok,
                            IntegerVector resp_ok,
                            int max_n,
                            int method,
                            double target,
                            double tox_upper,
                            double w1,
                            double w2,
                            double t_low,
                            double t_upp,
                            double e_low,
                            double e_upp,
                            double psi00,
                            double psi11,
                            int boinet_compatible) {

  const int n_trials = n_pts.nrow();
  const int n_doses = n_pts.ncol();

  const std::vector<int> tok(tox_ok.begin(), tox_ok.end());
  const std::vector<int> rok(resp_ok.begin(), resp_ok.end());

  simfastboin::EtObdRule rule;
  rule.method = method;
  rule.target = target;
  rule.tox_upper = tox_upper;
  rule.w1 = w1;
  rule.w2 = w2;
  rule.t_low = t_low;
  rule.t_upp = t_upp;
  rule.e_low = e_low;
  rule.e_upp = e_upp;
  rule.psi00 = psi00;
  rule.psi11 = psi11;

  IntegerVector obd(n_trials);
  IntegerVector reason(n_trials);
  std::vector<int> nn(n_doses), yy(n_doses), xx(n_doses);
  std::vector<double> pe(n_doses);

  for (int t = 0; t < n_trials; ++t) {
    for (int j = 0; j < n_doses; ++j) {
      nn[j] = n_pts(t, j);
      yy[j] = n_tox(t, j);
      xx[j] = n_resp(t, j);
      pe[j] = probe(t, j);
    }
    int why = 0;
    obd[t] = simfastboin::select_obd_et_one(&nn[0], &yy[0], &xx[0], &pe[0],
                                            n_doses, early_stop[t] == TRUE,
                                            tok, rok, max_n,
                                            boinet_compatible != 0, rule, why);
    reason[t] = why;
  }

  return List::create(
    _["obd"] = obd,
    _["reason"] = reason
  );
}

namespace {

// Event-time settings shared by the time-to-event designs.
simfastboin::EtTiming make_timing(double window, double resp_window,
                                  double accrual_rate, int arrival,
                                  int event_time, NumericVector shape_t,
                                  NumericVector scale_t, NumericVector shape_e,
                                  NumericVector scale_e, double min_share) {
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
  tm.min_share = min_share;
  return tm;
}

}  // namespace

// [[Rcpp::export]]
List tite_boin_et_simulate_cpp(int n_trials,
                               NumericVector p_true,
                               NumericVector p_resp,
                               double resp_cor,
                               IntegerVector cohort_size,
                               int start_dose,
                               int n_earlystop,
                               IntegerVector tox_ok,
                               IntegerVector resp_ok,
                               int max_n,
                               int boinet_compatible,
                               double lambda1,
                               double lambda2,
                               double eta1,
                               double window,
                               double resp_window,
                               double accrual_rate,
                               int arrival,
                               int event_time,
                               NumericVector shape_t,
                               NumericVector scale_t,
                               NumericVector shape_e,
                               NumericVector scale_e,
                               double min_completed_share) {

  simfastboin::EtDesign des;
  des.p_true.assign(p_true.begin(), p_true.end());
  des.p_resp.assign(p_resp.begin(), p_resp.end());
  des.cor_type = simfastboin::ET_COR_COPULA;
  des.resp_cor = resp_cor;
  des.cohort_size.assign(cohort_size.begin(), cohort_size.end());
  des.start_dose = start_dose;
  des.n_stop = n_earlystop;
  des.max_n = max_n;
  des.tox_ok.assign(tox_ok.begin(), tox_ok.end());
  des.resp_ok.assign(resp_ok.begin(), resp_ok.end());
  des.compat = boinet_compatible;
  const simfastboin::EtTiming tm = make_timing(
    window, resp_window, accrual_rate, arrival, event_time, shape_t, scale_t,
    shape_e, scale_e, min_completed_share);

  const int n_doses = p_true.size();
  IntegerMatrix n_pts(n_trials, n_doses);
  IntegerMatrix n_tox(n_trials, n_doses);
  IntegerMatrix n_resp(n_trials, n_doses);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);
  NumericVector duration(n_trials);

  std::vector<int> nn(n_doses), yy(n_doses), xx(n_doses);
  RUnifEt rng;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    double dur = 0.0;
    simfastboin::simulate_tite_et_one(des, tm, lambda1, lambda2, eta1, rng, nn,
                                      yy, xx, used, code, dur);
    for (int j = 0; j < n_doses; ++j) {
      n_pts(t, j) = nn[j];
      n_tox(t, j) = yy[j];
      n_resp(t, j) = xx[j];
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
    duration[t] = dur;
  }

  return List::create(
    _["n_pts"] = n_pts,
    _["n_tox"] = n_tox,
    _["n_resp"] = n_resp,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code,
    _["duration"] = duration
  );
}

// [[Rcpp::export]]
List bf_boin_et_simulate_cpp(int n_trials,
                             NumericVector p_true,
                             NumericVector p_resp,
                             double resp_cor,
                             IntegerVector cohort_size,
                             int start_dose,
                             double lambda1,
                             double lambda2,
                             double eta1,
                             int n_elimination,
                             int n_earlystop,
                             double suspend,
                             double extension,
                             IntegerVector tox_ok,
                             IntegerVector resp_ok,
                             int max_n,
                             double target,
                             double psi00,
                             double psi01,
                             double psi10,
                             double psi11,
                             double window,
                             double resp_window,
                             double accrual_rate,
                             int arrival,
                             int event_time,
                             NumericVector shape_t,
                             NumericVector scale_t,
                             NumericVector shape_e,
                             NumericVector scale_e) {

  simfastboin::BfEtDesign des;
  des.p_true.assign(p_true.begin(), p_true.end());
  des.p_resp.assign(p_resp.begin(), p_resp.end());
  des.resp_cor = resp_cor;
  des.cohort_size.assign(cohort_size.begin(), cohort_size.end());
  des.start_dose = start_dose;
  des.lambda1 = lambda1;
  des.lambda2 = lambda2;
  des.eta1 = eta1;
  des.n_elim = n_elimination;
  des.n_stop = n_earlystop;
  des.suspend = suspend;
  des.extension = extension;
  des.tox_ok.assign(tox_ok.begin(), tox_ok.end());
  des.resp_ok.assign(resp_ok.begin(), resp_ok.end());
  des.max_n = max_n;
  des.target = target;
  des.psi00 = psi00;
  des.psi01 = psi01;
  des.psi10 = psi10;
  des.psi11 = psi11;
  // The share of finished assessments is a rule of TITE-BOIN-ET only.
  const simfastboin::EtTiming tm = make_timing(
    window, resp_window, accrual_rate, arrival, event_time, shape_t, scale_t,
    shape_e, scale_e, 0.0);

  const int n_doses = p_true.size();
  IntegerMatrix n_pts(n_trials, n_doses);
  IntegerMatrix n_bf(n_trials, n_doses);
  IntegerMatrix n_tox(n_trials, n_doses);
  IntegerMatrix n_resp(n_trials, n_doses);
  IntegerMatrix n_tox_bf(n_trials, n_doses);
  IntegerMatrix n_resp_bf(n_trials, n_doses);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);
  LogicalVector early_stop(n_trials);
  IntegerVector obd(n_trials);
  NumericVector duration(n_trials);

  std::vector<int> nn(n_doses), nb(n_doses), yy(n_doses), xx(n_doses),
    yb(n_doses), xb(n_doses);
  RUnifEt rng;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    bool early = false;
    int sel = 0;
    double dur = 0.0;
    simfastboin::simulate_bf_et_one(des, tm, rng, nn, nb, yy, xx, yb, xb, used,
                                    code, early, sel, dur);
    for (int j = 0; j < n_doses; ++j) {
      n_pts(t, j) = nn[j];
      n_bf(t, j) = nb[j];
      n_tox(t, j) = yy[j];
      n_resp(t, j) = xx[j];
      n_tox_bf(t, j) = yb[j];
      n_resp_bf(t, j) = xb[j];
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
    early_stop[t] = early;
    obd[t] = sel;
    duration[t] = dur;
  }

  return List::create(
    _["n_pts"] = n_pts,
    _["n_bf"] = n_bf,
    _["n_tox"] = n_tox,
    _["n_resp"] = n_resp,
    _["n_tox_bf"] = n_tox_bf,
    _["n_resp_bf"] = n_resp_bf,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code,
    _["early_stop"] = early_stop,
    _["obd"] = obd,
    _["duration"] = duration
  );
}
