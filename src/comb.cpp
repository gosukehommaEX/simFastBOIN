// Rcpp glue for the BOIN drug-combination design. The algorithms themselves
// live in comb_core.h so that they can be compiled and verified outside of R
// as well. post_prob_above() is defined in boin.cpp.
//
// Trial data are three-dimensional arrays with dimensions (trial, level of
// drug A, level of drug B), stored as R stores them.

#include <Rcpp.h>
#include "comb_core.h"

using namespace Rcpp;

namespace simfastboin {

// Distribution function of the beta distribution, evaluated with R's
// incomplete beta function.
double beta_cdf(double x, double shape1, double shape2) {
  return R::pbeta(x, shape1, shape2, 1, 0);
}

}  // namespace simfastboin

namespace {

// Draws one uniform variate from R's random number stream.
struct CombUnif {
  double operator()() { return unif_rand(); }
};

IntegerVector trial_dims(int n_trials, int n_rows, int n_cols) {
  return IntegerVector::create(n_trials, n_rows, n_cols);
}

}  // namespace

// [[Rcpp::export]]
List comb_boin_simulate_cpp(int n_trials,
                            NumericVector p_true,
                            int n_rows,
                            int n_cols,
                            IntegerVector cohort_size,
                            int start_row,
                            int start_col,
                            int n_earlystop,
                            bool titration,
                            bool extrasafe,
                            double target,
                            double cutoff_eli,
                            double offset,
                            double lambda_e,
                            double lambda_d,
                            IntegerVector b_esc,
                            IntegerVector b_deesc,
                            IntegerVector b_elim) {

  const int n_cells = n_rows * n_cols;

  const std::vector<double> p(p_true.begin(), p_true.end());
  const std::vector<int> cs(cohort_size.begin(), cohort_size.end());
  const std::vector<int> be(b_esc.begin(), b_esc.end());
  const std::vector<int> bd(b_deesc.begin(), b_deesc.end());
  const std::vector<int> bl(b_elim.begin(), b_elim.end());

  IntegerVector n_pts(n_trials * n_cells);
  IntegerVector n_tox(n_trials * n_cells);
  LogicalVector eliminated(n_trials * n_cells);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nn(n_cells), yy(n_cells);
  std::vector<char> el(n_cells);
  CombUnif rng;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    simfastboin::comb_simulate_one(p, n_rows, n_cols, cs, start_row - 1,
                                   start_col - 1, n_earlystop, titration,
                                   extrasafe, target, cutoff_eli, offset,
                                   lambda_e, lambda_d, be, bd, bl, rng,
                                   nn, yy, el, used, code);
    for (int c = 0; c < n_cells; ++c) {
      n_pts[t + n_trials * c] = nn[c];
      n_tox[t + n_trials * c] = yy[c];
      eliminated[t + n_trials * c] = (el[c] == 1);
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
  }

  n_pts.attr("dim") = trial_dims(n_trials, n_rows, n_cols);
  n_tox.attr("dim") = trial_dims(n_trials, n_rows, n_cols);
  eliminated.attr("dim") = trial_dims(n_trials, n_rows, n_cols);

  return List::create(
    _["n_pts"] = n_pts,
    _["n_tox"] = n_tox,
    _["eliminated"] = eliminated,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code
  );
}

// [[Rcpp::export]]
LogicalVector comb_boin_eliminate_cpp(IntegerVector n_pts,
                                      IntegerVector n_tox,
                                      double target,
                                      double cutoff_eli,
                                      bool extrasafe,
                                      double offset) {

  const IntegerVector dims = n_pts.attr("dim");
  const int n_trials = dims[0];
  const int n_rows = dims[1];
  const int n_cols = dims[2];
  const int n_cells = n_rows * n_cols;

  LogicalVector out(n_trials * n_cells);
  std::vector<int> nn(n_cells), yy(n_cells);
  std::vector<char> el(n_cells);

  for (int t = 0; t < n_trials; ++t) {
    for (int c = 0; c < n_cells; ++c) {
      nn[c] = n_pts[t + n_trials * c];
      yy[c] = n_tox[t + n_trials * c];
    }
    simfastboin::comb_selection_elimination(&nn[0], &yy[0], n_rows, n_cols,
                                            target, cutoff_eli, extrasafe,
                                            offset, &el[0]);
    for (int c = 0; c < n_cells; ++c) out[t + n_trials * c] = (el[c] == 1);
  }

  out.attr("dim") = trial_dims(n_trials, n_rows, n_cols);
  return out;
}

// [[Rcpp::export]]
List comb_boin_choose_cpp(NumericVector iso,
                          IntegerVector n_pts,
                          LogicalVector eliminated,
                          double target,
                          bool bound_mtd,
                          double lambda_d) {

  const IntegerVector dims = n_pts.attr("dim");
  const int n_trials = dims[0];
  const int n_rows = dims[1];
  const int n_cols = dims[2];
  const int n_cells = n_rows * n_cols;

  IntegerVector mtd_a(n_trials);
  IntegerVector mtd_b(n_trials);
  IntegerVector reason(n_trials);

  std::vector<double> ii(n_cells);
  std::vector<int> nn(n_cells);
  std::vector<char> el(n_cells);

  for (int t = 0; t < n_trials; ++t) {
    for (int c = 0; c < n_cells; ++c) {
      ii[c] = iso[t + n_trials * c];
      nn[c] = n_pts[t + n_trials * c];
      el[c] = (eliminated[t + n_trials * c] == TRUE) ? 1 : 0;
    }
    int row = 0;
    int col = 0;
    const int code = simfastboin::comb_choose_mtd(&ii[0], &nn[0], &el[0],
                                                  n_rows, n_cols, target,
                                                  bound_mtd, lambda_d,
                                                  row, col);
    mtd_a[t] = (row == 0) ? NA_INTEGER : row;
    mtd_b[t] = (col == 0) ? NA_INTEGER : col;
    reason[t] = code;
  }

  return List::create(_["mtd_a"] = mtd_a, _["mtd_b"] = mtd_b,
                      _["reason"] = reason);
}
