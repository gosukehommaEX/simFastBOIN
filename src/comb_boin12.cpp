// Rcpp glue for the core algorithms of the Comb-BOIN12 design. The
// algorithms themselves live in comb_boin12_core.h so that they can be
// compiled and verified outside of R as well.
//
// Trial data are three-dimensional arrays with dimensions (trial, level of
// drug A, level of drug B), stored as R stores them.

#include <Rcpp.h>
#include <cmath>
#include "comb_boin12_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnifCb12 {
  double operator()() { return unif_rand(); }
};

// Upper tail Pr(Beta(a, b) > q) of the beta distribution.
struct RPbetaUpperCb12 {
  double operator()(double q, double a, double b) const {
    return R::pbeta(q, a, b, 0, 0);
  }
};

IntegerVector cb12_dims(int n_trials, int n_rows, int n_cols) {
  return IntegerVector::create(n_trials, n_rows, n_cols);
}

simfastboin::Cb12Design make_cb12_design(
    NumericVector p_true, NumericVector p_resp, NumericVector q0,
    NumericVector q1, int cor_type, double resp_cor, int n_rows, int n_cols,
    IntegerVector cohort_size, int start_row, int start_col, int n_earlystop,
    bool early_stop_simple, int n_star, int tie_break, int no_candidate,
    double lambda_e, double lambda_d, double psi00, double psi11,
    double alpha, double beta, double u_benchmark, int max_n,
    IntegerVector tox_ok, IntegerVector tox_ok_lowest,
    IntegerVector resp_ok) {
  simfastboin::Cb12Design des;
  simfastboin::B12Design& b = des.b12;
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
  b.max_n = max_n;
  b.tox_ok.assign(tox_ok.begin(), tox_ok.end());
  b.resp_ok.assign(resp_ok.begin(), resp_ok.end());
  des.n_rows = n_rows;
  des.n_cols = n_cols;
  des.start_row = start_row - 1;
  des.start_col = start_col - 1;
  des.tox_ok_lowest.assign(tox_ok_lowest.begin(), tox_ok_lowest.end());
  des.both_weight = psi00 + psi11 - 100.0;
  if (std::fabs(des.both_weight) < 1e-9) des.both_weight = 0.0;
  return des;
}

simfastboin::Cb12Logistic make_cb12_logistic(NumericVector fa,
                                             NumericVector fb,
                                             double scale0, double scale1,
                                             NumericMatrix z,
                                             NumericVector w, double df,
                                             int rounds) {
  simfastboin::Cb12Logistic lg;
  lg.fa.assign(fa.begin(), fa.end());
  lg.fb.assign(fb.begin(), fb.end());
  lg.scale0 = scale0;
  lg.scale1 = scale1;
  const int n_draws = z.nrow();
  lg.z.resize(n_draws * 5);
  for (int d = 0; d < n_draws; ++d) {
    for (int k = 0; k < 5; ++k) lg.z[d * 5 + k] = z(d, k);
  }
  lg.w.assign(w.begin(), w.end());
  lg.df = df;
  lg.rounds = rounds;
  return lg;
}

}  // namespace

// [[Rcpp::export]]
List comb_boin12_simulate_cpp(int n_trials, NumericVector p_true,
                              NumericVector p_resp, NumericVector q0,
                              NumericVector q1, int cor_type,
                              double resp_cor, int n_rows, int n_cols,
                              IntegerVector cohort_size, int start_row,
                              int start_col, int n_earlystop,
                              bool early_stop_simple, int n_star,
                              int tie_break, int no_candidate,
                              double lambda_e, double lambda_d, double psi00,
                              double psi11, double alpha, double beta,
                              double u_benchmark, int max_n,
                              IntegerVector tox_ok,
                              IntegerVector tox_ok_lowest,
                              IntegerVector resp_ok) {

  const simfastboin::Cb12Design des = make_cb12_design(
    p_true, p_resp, q0, q1, cor_type, resp_cor, n_rows, n_cols, cohort_size,
    start_row, start_col, n_earlystop, early_stop_simple, n_star, tie_break,
    no_candidate, lambda_e, lambda_d, psi00, psi11, alpha, beta, u_benchmark,
    max_n, tox_ok, tox_ok_lowest, resp_ok);
  const int n_cells = n_rows * n_cols;

  IntegerVector n_pts(n_trials * n_cells);
  IntegerVector n_tox(n_trials * n_cells);
  IntegerVector n_resp(n_trials * n_cells);
  IntegerVector n_tox_resp(n_trials * n_cells);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nn(n_cells), yy(n_cells), xx(n_cells), zz(n_cells);
  RUnifCb12 rng;
  RPbetaUpperCb12 pb;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    simfastboin::cb12_simulate_one(des, rng, pb, nn, yy, xx, zz, used, code);
    for (int c = 0; c < n_cells; ++c) {
      n_pts[t + n_trials * c] = nn[c];
      n_tox[t + n_trials * c] = yy[c];
      n_resp[t + n_trials * c] = xx[c];
      n_tox_resp[t + n_trials * c] = zz[c];
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
  }

  n_pts.attr("dim") = cb12_dims(n_trials, n_rows, n_cols);
  n_tox.attr("dim") = cb12_dims(n_trials, n_rows, n_cols);
  n_resp.attr("dim") = cb12_dims(n_trials, n_rows, n_cols);
  n_tox_resp.attr("dim") = cb12_dims(n_trials, n_rows, n_cols);

  return List::create(
    _["n_pts"] = n_pts,
    _["n_tox"] = n_tox,
    _["n_resp"] = n_resp,
    _["n_tox_resp"] = n_tox_resp,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code
  );
}

// [[Rcpp::export]]
NumericMatrix comb_boin12_logistic_cpp(IntegerMatrix n_pts,
                                       IntegerMatrix n_resp,
                                       NumericVector fa, NumericVector fb,
                                       double scale0, double scale1,
                                       NumericMatrix z, NumericVector w,
                                       double df, int rounds) {
  const int n_rows = n_pts.nrow();
  const int n_cols = n_pts.ncol();
  const int n_cells = n_rows * n_cols;
  const simfastboin::Cb12Logistic lg = make_cb12_logistic(
    fa, fb, scale0, scale1, z, w, df, rounds);
  std::vector<int> nn(n_pts.begin(), n_pts.end());
  std::vector<int> xx(n_resp.begin(), n_resp.end());
  std::vector<double> pe(n_cells);
  simfastboin::cb12_logistic_mean(&nn[0], &xx[0], n_cells, lg, &pe[0]);
  NumericMatrix out(n_rows, n_cols);
  for (int c = 0; c < n_cells; ++c) out[c] = pe[c];
  return out;
}

// [[Rcpp::export]]
List comb_boin12_select_cpp(IntegerVector n_pts, IntegerVector n_tox,
                            IntegerVector n_resp,
                            LogicalVector early_stop, double target,
                            double psi00, double psi11, int max_n,
                            IntegerVector tox_ok, IntegerVector tox_ok_lowest,
                            IntegerVector resp_ok, int est, int bound,
                            NumericVector fa, NumericVector fb,
                            double scale0, double scale1, NumericMatrix z,
                            NumericVector w, double df, int rounds) {

  const IntegerVector dims = n_pts.attr("dim");
  const int n_trials = dims[0];
  const int n_rows = dims[1];
  const int n_cols = dims[2];
  const int n_cells = n_rows * n_cols;

  NumericVector empty(n_cells);
  IntegerVector no_cohort(1, 1);
  const simfastboin::Cb12Design des = make_cb12_design(
    empty, empty, empty, empty, 0, 0.0, n_rows, n_cols, no_cohort, 1, 1, 0,
    false, 6, 0, 0, 0.0, 1.0, psi00, psi11, 1.0, 1.0, 0.5, max_n, tox_ok,
    tox_ok_lowest, resp_ok);
  const simfastboin::Cb12Logistic lg = make_cb12_logistic(
    fa, fb, scale0, scale1, z, w, df, rounds);

  IntegerVector obd_a(n_trials), obd_b(n_trials);
  IntegerVector mtd_a(n_trials), mtd_b(n_trials);
  IntegerVector reason(n_trials);
  NumericVector utility(n_trials * n_cells, NA_REAL);
  NumericVector resp_est(n_trials * n_cells, NA_REAL);
  NumericVector tox_est(n_trials * n_cells);
  int not_converged = 0;

  std::vector<int> nn(n_cells), yy(n_cells), xx(n_cells);
  std::vector<double> ii(n_cells), pe(n_cells), util(n_cells);
  std::vector<char> adm(n_cells);

  for (int t = 0; t < n_trials; ++t) {
    for (int c = 0; c < n_cells; ++c) {
      nn[c] = n_pts[t + n_trials * c];
      yy[c] = n_tox[t + n_trials * c];
      xx[c] = n_resp[t + n_trials * c];
      util[c] = NA_REAL;
      pe[c] = NA_REAL;
    }
    if (!simfastboin::cb12_biviso(&nn[0], &yy[0], n_rows, n_cols, &ii[0])) {
      ++not_converged;
    }
    const bool stopped = early_stop[t] == TRUE;
    if (!stopped &&
        simfastboin::cb12_admissible_set(&nn[0], &yy[0], &xx[0], des,
                                         &adm[0])) {
      if (est == 0) {
        simfastboin::cb12_logistic_mean(&nn[0], &xx[0], n_cells, lg, &pe[0]);
      } else {
        for (int c = 0; c < n_cells; ++c) {
          pe[c] = nn[c] > 0 ? static_cast<double>(xx[c]) / nn[c] : NA_REAL;
        }
      }
    }
    int why = 0;
    int mtdc = -1;
    const int pick = simfastboin::cb12_select_obdc_one(
      &nn[0], &yy[0], &xx[0], &ii[0], &pe[0], stopped, des, target, bound,
      why, mtdc, &util[0]);
    obd_a[t] = pick < 0 ? NA_INTEGER : pick % n_rows + 1;
    obd_b[t] = pick < 0 ? NA_INTEGER : pick / n_rows + 1;
    mtd_a[t] = mtdc < 0 ? NA_INTEGER : mtdc % n_rows + 1;
    mtd_b[t] = mtdc < 0 ? NA_INTEGER : mtdc / n_rows + 1;
    reason[t] = why;
    for (int c = 0; c < n_cells; ++c) {
      tox_est[t + n_trials * c] = ii[c];
      if (nn[c] > 0) {
        utility[t + n_trials * c] = util[c];
        resp_est[t + n_trials * c] = pe[c];
      }
    }
  }

  utility.attr("dim") = cb12_dims(n_trials, n_rows, n_cols);
  resp_est.attr("dim") = cb12_dims(n_trials, n_rows, n_cols);
  tox_est.attr("dim") = cb12_dims(n_trials, n_rows, n_cols);

  return List::create(
    _["obd_a"] = obd_a,
    _["obd_b"] = obd_b,
    _["mtd_a"] = mtd_a,
    _["mtd_b"] = mtd_b,
    _["reason"] = reason,
    _["utility"] = utility,
    _["resp_est"] = resp_est,
    _["tox_est"] = tox_est,
    _["not_converged"] = not_converged
  );
}
