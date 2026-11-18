// Rcpp glue for the core algorithm of the TITE-Comb-BOIN12 design. The
// algorithm itself lives in tite_comb_boin12_core.h so that it can be
// compiled and verified outside of R as well.
//
// Trial data are three-dimensional arrays with dimensions (trial, level of
// drug A, level of drug B), stored as R stores them.

#include <Rcpp.h>
#include <cmath>
#include "tite_comb_boin12_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnifTcb12 {
  double operator()() { return unif_rand(); }
};

// Lower or upper tail of the beta distribution.
struct RPbetaTcb12 {
  double operator()(double q, double a, double b, bool lower) const {
    return R::pbeta(q, a, b, lower ? 1 : 0, 0);
  }
};

}  // namespace

// [[Rcpp::export]]
List tite_comb_boin12_simulate_cpp(int n_trials, NumericVector p_true,
                                   NumericVector p_resp, NumericVector q0,
                                   NumericVector q1, int cor_type,
                                   double resp_cor, int n_rows, int n_cols,
                                   IntegerVector cohort_size, int start_row,
                                   int start_col, int n_earlystop,
                                   bool early_stop_simple, int n_star,
                                   int tie_break, int no_candidate,
                                   double lambda_e, double lambda_d,
                                   double psi00, double psi11, double alpha,
                                   double beta, double u_benchmark,
                                   double target, double target_resp,
                                   double cutoff_eli, double cutoff_lowest,
                                   double cutoff_eli_resp, int pending,
                                   int pending_use,
                                   double window, double resp_window,
                                   double accrual_rate, int arrival,
                                   int event_time, NumericVector shape_t,
                                   NumericVector scale_t,
                                   NumericVector shape_e,
                                   NumericVector scale_e,
                                   double min_completed_share) {

  simfastboin::Tcb12Design des;
  simfastboin::Cb12Design& cb = des.cb;
  simfastboin::B12Design& b = cb.b12;
  b.p_true.assign(p_true.begin(), p_true.end());
  b.p_resp.assign(p_resp.begin(), p_resp.end());
  b.q0.assign(q0.begin(), q0.end());
  b.q1.assign(q1.begin(), q1.end());
  b.cor_type = cor_type;
  b.resp_cor = resp_cor;
  b.cohort_size.assign(cohort_size.begin(), cohort_size.end());
  b.start_dose = 1;
  b.n_earlystop = n_earlystop;
  b.early_stop_simple = early_stop_simple ? 1 : 0;
  b.n_star = n_star;
  b.n_explore = 0;
  b.tie_random = tie_break == 1 ? 1 : 0;
  b.stop_no_candidate = no_candidate == 1 ? 1 : 0;
  b.explore_admissible_only = 0;
  b.lambda_e = lambda_e;
  b.lambda_d = lambda_d;
  b.psi00 = psi00;
  b.psi11 = psi11;
  b.alpha = alpha;
  b.beta = beta;
  b.u_benchmark = u_benchmark;
  b.max_n = 0;
  cb.n_rows = n_rows;
  cb.n_cols = n_cols;
  cb.start_row = start_row - 1;
  cb.start_col = start_col - 1;
  cb.both_weight = psi00 + psi11 - 100.0;
  if (std::fabs(cb.both_weight) < 1e-9) cb.both_weight = 0.0;
  des.target = target;
  des.target_resp = target_resp;
  des.cutoff_eli = cutoff_eli;
  des.cutoff_lowest = cutoff_lowest;
  des.cutoff_eli_resp = cutoff_eli_resp;
  des.pending = pending;
  des.pending_use = pending_use;

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

  const int n_cells = n_rows * n_cols;
  IntegerVector n_pts(n_trials * n_cells);
  IntegerVector n_tox(n_trials * n_cells);
  IntegerVector n_resp(n_trials * n_cells);
  IntegerVector n_tox_resp(n_trials * n_cells);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);
  NumericVector duration(n_trials);

  std::vector<int> nn(n_cells), yy(n_cells), xx(n_cells), zz(n_cells);
  std::vector<simfastboin::Tb12Patient> pts;
  RUnifTcb12 rng;
  RPbetaTcb12 pb;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    double dur = 0.0;
    simfastboin::tcb12_simulate_one(des, tm, rng, pb, pts, nn, yy, xx, zz,
                                    used, code, dur);
    for (int c = 0; c < n_cells; ++c) {
      n_pts[t + n_trials * c] = nn[c];
      n_tox[t + n_trials * c] = yy[c];
      n_resp[t + n_trials * c] = xx[c];
      n_tox_resp[t + n_trials * c] = zz[c];
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
    duration[t] = dur;
  }

  const IntegerVector dims = IntegerVector::create(n_trials, n_rows, n_cols);
  n_pts.attr("dim") = dims;
  n_tox.attr("dim") = dims;
  n_resp.attr("dim") = dims;
  n_tox_resp.attr("dim") = dims;

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
