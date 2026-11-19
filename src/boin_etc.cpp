// Rcpp glue for the core algorithms of the BOIN-ETC design. The algorithms
// themselves live in boin_etc_core.h so that they can be compiled and
// verified outside of R as well.
//
// Trial data are three-dimensional arrays with dimensions (trial, level of
// drug A, level of drug B), stored as R stores them.

#include <Rcpp.h>
#include <cmath>
#include "boin_etc_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnifEtc {
  double operator()() { return unif_rand(); }
};

IntegerVector etc_dims(int n_trials, int n_rows, int n_cols) {
  return IntegerVector::create(n_trials, n_rows, n_cols);
}

}  // namespace

// [[Rcpp::export]]
List boin_etc_simulate_cpp(int n_trials, NumericVector p_true,
                           NumericVector p_resp, NumericVector q0,
                           NumericVector q1, int cor_type, double resp_cor,
                           int n_rows, int n_cols, IntegerVector cohort_size,
                           int subtrials, int start_rule, int rule4,
                           int n_stop, int max_n, IntegerVector dec,
                           IntegerVector tox_ok) {

  simfastboin::EtcDesign des;
  des.n_rows = n_rows;
  des.n_cols = n_cols;
  des.p_true.assign(p_true.begin(), p_true.end());
  des.p_resp.assign(p_resp.begin(), p_resp.end());
  des.q0.assign(q0.begin(), q0.end());
  des.q1.assign(q1.begin(), q1.end());
  des.cor_type = cor_type;
  des.resp_cor = resp_cor;
  des.cohort_size.assign(cohort_size.begin(), cohort_size.end());
  des.subtrials = subtrials;
  des.start_rule = start_rule;
  des.rule4 = rule4;
  des.n_stop = n_stop;
  des.max_n = max_n;
  des.dec.assign(dec.begin(), dec.end());
  des.tox_ok.assign(tox_ok.begin(), tox_ok.end());
  const int n_cells = n_rows * n_cols;

  IntegerVector n_pts(n_trials * n_cells);
  IntegerVector n_tox(n_trials * n_cells);
  IntegerVector n_resp(n_trials * n_cells);
  IntegerVector n_tox_resp(n_trials * n_cells);
  IntegerVector excluded(n_trials * n_cells);
  IntegerVector cohorts_used(n_trials);
  IntegerVector n_subtrials(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nn(n_cells), yy(n_cells), xx(n_cells), zz(n_cells);
  std::vector<char> ex(n_cells);
  RUnifEtc rng;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int subs = 0;
    int code = 0;
    simfastboin::etc_simulate_one(des, rng, nn, yy, xx, zz, ex, used, subs,
                                  code);
    for (int c = 0; c < n_cells; ++c) {
      n_pts[t + n_trials * c] = nn[c];
      n_tox[t + n_trials * c] = yy[c];
      n_resp[t + n_trials * c] = xx[c];
      n_tox_resp[t + n_trials * c] = zz[c];
      excluded[t + n_trials * c] = ex[c];
    }
    cohorts_used[t] = used;
    n_subtrials[t] = subs;
    stop_code[t] = code;
  }

  n_pts.attr("dim") = etc_dims(n_trials, n_rows, n_cols);
  n_tox.attr("dim") = etc_dims(n_trials, n_rows, n_cols);
  n_resp.attr("dim") = etc_dims(n_trials, n_rows, n_cols);
  n_tox_resp.attr("dim") = etc_dims(n_trials, n_rows, n_cols);
  excluded.attr("dim") = etc_dims(n_trials, n_rows, n_cols);

  return List::create(
    _["n_pts"] = n_pts,
    _["n_tox"] = n_tox,
    _["n_resp"] = n_resp,
    _["n_tox_resp"] = n_tox_resp,
    _["excluded"] = excluded,
    _["cohorts_used"] = cohorts_used,
    _["n_subtrials"] = n_subtrials,
    _["stop_code"] = stop_code
  );
}

// [[Rcpp::export]]
List boin_etc_select_cpp(IntegerVector n_pts, IntegerVector n_tox,
                         IntegerVector n_resp, IntegerVector n_tox_resp,
                         IntegerVector excluded, LogicalVector early_stop,
                         double target, double psi00, double both_weight,
                         double alpha, double beta) {

  const IntegerVector dims = n_pts.attr("dim");
  const int n_trials = dims[0];
  const int n_rows = dims[1];
  const int n_cols = dims[2];
  const int n_cells = n_rows * n_cols;

  IntegerVector obd_a(n_trials), obd_b(n_trials), reason(n_trials);
  IntegerVector n_tied(n_trials);
  IntegerVector tied(n_trials * n_cells);
  IntegerVector mtd_col(n_trials * n_rows);
  NumericVector utility(n_trials * n_cells, NA_REAL);
  NumericVector tox_est(n_trials * n_cells);
  int not_converged = 0;

  std::vector<int> nn(n_cells), yy(n_cells), xx(n_cells), zz(n_cells);
  std::vector<char> ex(n_cells), tt(n_cells);
  std::vector<double> ii(n_cells), util(n_cells);
  std::vector<int> mc(n_rows);

  for (int t = 0; t < n_trials; ++t) {
    for (int c = 0; c < n_cells; ++c) {
      nn[c] = n_pts[t + n_trials * c];
      yy[c] = n_tox[t + n_trials * c];
      xx[c] = n_resp[t + n_trials * c];
      zz[c] = n_tox_resp[t + n_trials * c];
      ex[c] = excluded[t + n_trials * c] != 0;
      util[c] = NA_REAL;
    }
    if (!simfastboin::etc_biviso(&nn[0], &yy[0], n_rows, n_cols, &ii[0])) {
      ++not_converged;
    }
    int why = 0;
    const int pick = simfastboin::etc_select_one(
      &nn[0], &yy[0], &xx[0], &zz[0], &ex[0], &ii[0], early_stop[t] == TRUE,
      n_rows, n_cols, target, psi00, both_weight, alpha, beta, why, &mc[0],
      &tt[0], &util[0]);
    obd_a[t] = pick < 0 ? NA_INTEGER : pick % n_rows + 1;
    obd_b[t] = pick < 0 ? NA_INTEGER : pick / n_rows + 1;
    reason[t] = why;
    int count = 0;
    for (int c = 0; c < n_cells; ++c) {
      tox_est[t + n_trials * c] = ii[c];
      tied[t + n_trials * c] = tt[c];
      count += tt[c];
      if (nn[c] > 0) utility[t + n_trials * c] = util[c];
    }
    n_tied[t] = count;
    for (int a = 0; a < n_rows; ++a) {
      mtd_col[t + n_trials * a] = mc[a] < 0 ? NA_INTEGER : mc[a] + 1;
    }
  }

  tied.attr("dim") = etc_dims(n_trials, n_rows, n_cols);
  utility.attr("dim") = etc_dims(n_trials, n_rows, n_cols);
  tox_est.attr("dim") = etc_dims(n_trials, n_rows, n_cols);
  mtd_col.attr("dim") = IntegerVector::create(n_trials, n_rows);

  return List::create(
    _["obd_a"] = obd_a,
    _["obd_b"] = obd_b,
    _["reason"] = reason,
    _["n_tied"] = n_tied,
    _["tied"] = tied,
    _["mtd_col"] = mtd_col,
    _["utility"] = utility,
    _["tox_est"] = tox_est,
    _["not_converged"] = not_converged
  );
}
