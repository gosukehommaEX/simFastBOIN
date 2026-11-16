// Rcpp glue for the core algorithms of the BOIN12 design. The algorithms
// themselves live in boin12_core.h so that they can be compiled and verified
// outside of R as well.

#include <Rcpp.h>
#include "boin12_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnif12 {
  double operator()() { return unif_rand(); }
};

// Upper tail Pr(Beta(a, b) > q) of the beta distribution.
struct RPbetaUpper {
  double operator()(double q, double a, double b) const {
    return R::pbeta(q, a, b, 0, 0);
  }
};

simfastboin::B12Design make_b12_design(
    NumericVector p_true, NumericVector p_resp, NumericVector q0,
    NumericVector q1, int cor_type, double resp_cor, IntegerVector cohort_size,
    int start_dose, int n_earlystop, bool early_stop_simple, int n_star,
    int n_explore, double lambda_e, double lambda_d, double psi00,
    double psi11, double alpha, double beta, double u_benchmark, int max_n,
    IntegerVector tox_ok, IntegerVector resp_ok) {
  simfastboin::B12Design des;
  des.p_true.assign(p_true.begin(), p_true.end());
  des.p_resp.assign(p_resp.begin(), p_resp.end());
  des.q0.assign(q0.begin(), q0.end());
  des.q1.assign(q1.begin(), q1.end());
  des.cor_type = cor_type;
  des.resp_cor = resp_cor;
  des.cohort_size.assign(cohort_size.begin(), cohort_size.end());
  des.start_dose = start_dose;
  des.n_earlystop = n_earlystop;
  des.early_stop_simple = early_stop_simple ? 1 : 0;
  des.n_star = n_star;
  des.n_explore = n_explore;
  des.lambda_e = lambda_e;
  des.lambda_d = lambda_d;
  des.psi00 = psi00;
  des.psi11 = psi11;
  des.alpha = alpha;
  des.beta = beta;
  des.u_benchmark = u_benchmark;
  des.max_n = max_n;
  des.tox_ok.assign(tox_ok.begin(), tox_ok.end());
  des.resp_ok.assign(resp_ok.begin(), resp_ok.end());
  return des;
}

}  // namespace

// [[Rcpp::export]]
List boin12_simulate_cpp(int n_trials, NumericVector p_true,
                         NumericVector p_resp, NumericVector q0,
                         NumericVector q1, int cor_type, double resp_cor,
                         IntegerVector cohort_size, int start_dose,
                         int n_earlystop, bool early_stop_simple, int n_star,
                         int n_explore, double lambda_e, double lambda_d,
                         double psi00, double psi11, double alpha,
                         double beta, double u_benchmark, int max_n,
                         IntegerVector tox_ok, IntegerVector resp_ok) {

  const simfastboin::B12Design des = make_b12_design(
    p_true, p_resp, q0, q1, cor_type, resp_cor, cohort_size, start_dose,
    n_earlystop, early_stop_simple, n_star, n_explore, lambda_e, lambda_d,
    psi00, psi11, alpha, beta, u_benchmark, max_n, tox_ok, resp_ok);
  const int n_doses = p_true.size();

  IntegerMatrix n_pts(n_trials, n_doses);
  IntegerMatrix n_tox(n_trials, n_doses);
  IntegerMatrix n_resp(n_trials, n_doses);
  IntegerMatrix n_tox_resp(n_trials, n_doses);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nn(n_doses), yy(n_doses), xx(n_doses), zz(n_doses);
  RUnif12 rng;
  RPbetaUpper pb;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    simfastboin::b12_simulate_one(des, rng, pb, nn, yy, xx, zz, used, code);
    for (int j = 0; j < n_doses; ++j) {
      n_pts(t, j) = nn[j];
      n_tox(t, j) = yy[j];
      n_resp(t, j) = xx[j];
      n_tox_resp(t, j) = zz[j];
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
  }

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
List boin12_select_obd_cpp(IntegerMatrix n_pts, IntegerMatrix n_tox,
                           IntegerMatrix n_resp, IntegerMatrix n_tox_resp,
                           LogicalVector early_stop, double target,
                           double psi00, double psi11, double alpha,
                           double beta, double u_benchmark, int est,
                           int max_n, IntegerVector tox_ok,
                           IntegerVector resp_ok) {

  const int n_trials = n_pts.nrow();
  const int n_doses = n_pts.ncol();
  NumericVector empty(n_doses);
  IntegerVector no_cohort(1, 1);
  const simfastboin::B12Design des = make_b12_design(
    empty, empty, empty, empty, 0, 0.0, no_cohort, 1, 0, false, 6, 9, 0.0,
    1.0, psi00, psi11, alpha, beta, u_benchmark, max_n, tox_ok, resp_ok);
  RPbetaUpper pb;

  IntegerVector obd(n_trials);
  IntegerVector reason(n_trials);
  NumericMatrix utility(n_trials, n_doses);
  std::fill(utility.begin(), utility.end(), NA_REAL);

  std::vector<int> nn(n_doses), yy(n_doses), xx(n_doses), zz(n_doses);
  std::vector<double> util(n_doses);

  for (int t = 0; t < n_trials; ++t) {
    for (int j = 0; j < n_doses; ++j) {
      nn[j] = n_pts(t, j);
      yy[j] = n_tox(t, j);
      xx[j] = n_resp(t, j);
      zz[j] = n_tox_resp(t, j);
      util[j] = NA_REAL;
    }
    int why = 0;
    const int pick = simfastboin::b12_select_obd_one(
      &nn[0], &yy[0], &xx[0], &zz[0], n_doses, early_stop[t] == TRUE, des,
      target, est, pb, why, &util[0]);
    obd[t] = pick == 0 ? NA_INTEGER : pick;
    reason[t] = why;
    for (int j = 0; j < n_doses; ++j) utility(t, j) = util[j];
  }

  return List::create(
    _["obd"] = obd,
    _["reason"] = reason,
    _["utility"] = utility
  );
}
