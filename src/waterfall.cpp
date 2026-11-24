// Rcpp glue for the waterfall design. The algorithms themselves live in
// waterfall_core.h so that they can be compiled and verified outside of R as
// well. post_prob_above() is defined in boin.cpp and beta_cdf() in comb.cpp.
//
// Trial data are three-dimensional arrays with dimensions (trial, level of
// drug A, level of drug B), stored as R stores them.

#include <Rcpp.h>
#include "waterfall_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct WaterfallUnif {
  double operator()() { return unif_rand(); }
};

IntegerVector waterfall_dims(int n_trials, int n_rows, int n_cols) {
  return IntegerVector::create(n_trials, n_rows, n_cols);
}

}  // namespace

// [[Rcpp::export]]
List waterfall_simulate_cpp(int n_trials,
                            NumericVector p_true,
                            int n_rows,
                            int n_cols,
                            IntegerVector n_cohort,
                            int cohort_size,
                            int n_earlystop,
                            bool titration,
                            bool extrasafe,
                            bool bound_mtd,
                            bool boin_compatible,
                            double target,
                            double cutoff_eli,
                            double offset,
                            double lambda_d,
                            IntegerVector b_esc,
                            IntegerVector b_deesc,
                            IntegerVector b_elim) {

  const int n_cells = n_rows * n_cols;

  simfastboin::WaterfallDesign des;
  des.n_rows = n_rows;
  des.n_cols = n_cols;
  des.p_true.assign(p_true.begin(), p_true.end());
  des.n_cohort.assign(n_cohort.begin(), n_cohort.end());
  des.cohort_size = cohort_size;
  des.n_earlystop = n_earlystop;
  des.titration = titration;
  des.extrasafe = extrasafe;
  des.bound_mtd = bound_mtd;
  des.boin_compatible = boin_compatible;
  des.target = target;
  des.cutoff_eli = cutoff_eli;
  des.offset = offset;
  des.lambda_d = lambda_d;
  des.b_esc.assign(b_esc.begin(), b_esc.end());
  des.b_deesc.assign(b_deesc.begin(), b_deesc.end());
  des.b_elim.assign(b_elim.begin(), b_elim.end());

  IntegerVector n_pts(n_trials * n_cells);
  IntegerVector n_tox(n_trials * n_cells);
  LogicalVector eliminated(n_trials * n_cells);
  IntegerVector n_subtrials(n_trials);
  IntegerVector n_enrolled(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nn(n_cells), yy(n_cells);
  std::vector<char> el(n_cells);
  WaterfallUnif rng;

  for (int t = 0; t < n_trials; ++t) {
    int subtrials = 0;
    int enrolled = 0;
    int code = 0;
    simfastboin::waterfall_simulate_one(des, rng, nn, yy, el, subtrials,
                                        enrolled, code);
    for (int c = 0; c < n_cells; ++c) {
      n_pts[t + n_trials * c] = nn[c];
      n_tox[t + n_trials * c] = yy[c];
      eliminated[t + n_trials * c] = (el[c] == 1);
    }
    n_subtrials[t] = subtrials;
    n_enrolled[t] = enrolled;
    stop_code[t] = code;
  }

  n_pts.attr("dim") = waterfall_dims(n_trials, n_rows, n_cols);
  n_tox.attr("dim") = waterfall_dims(n_trials, n_rows, n_cols);
  eliminated.attr("dim") = waterfall_dims(n_trials, n_rows, n_cols);

  return List::create(
    _["n_pts"] = n_pts,
    _["n_tox"] = n_tox,
    _["eliminated"] = eliminated,
    _["n_subtrials"] = n_subtrials,
    _["n_enrolled"] = n_enrolled,
    _["stop_code"] = stop_code
  );
}

// [[Rcpp::export]]
IntegerMatrix waterfall_choose_cpp(NumericVector iso,
                                   IntegerVector n_pts,
                                   LogicalVector eliminated,
                                   double target) {

  const IntegerVector dims = n_pts.attr("dim");
  const int n_trials = dims[0];
  const int n_rows = dims[1];
  const int n_cols = dims[2];
  const int n_cells = n_rows * n_cols;

  IntegerMatrix mtd(n_trials, n_rows);

  std::vector<double> ii(n_cells);
  std::vector<int> nn(n_cells), col(n_rows);
  std::vector<char> el(n_cells);

  for (int t = 0; t < n_trials; ++t) {
    bool has_iso = true;
    for (int c = 0; c < n_cells; ++c) {
      ii[c] = iso[t + n_trials * c];
      if (ISNAN(ii[c])) has_iso = false;
      nn[c] = n_pts[t + n_trials * c];
      el[c] = (eliminated[t + n_trials * c] == TRUE) ? 1 : 0;
    }
    if (has_iso) {
      simfastboin::waterfall_choose_mtd(&ii[0], &nn[0], &el[0], n_rows,
                                        n_cols, target, &col[0]);
    } else {
      std::fill(col.begin(), col.end(), 0);
    }
    for (int a = 0; a < n_rows; ++a) {
      mtd(t, a) = (col[a] == 0) ? NA_INTEGER : col[a];
    }
  }

  return mtd;
}
