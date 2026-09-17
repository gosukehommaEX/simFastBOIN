// Rcpp glue for the BOIN core algorithms. The algorithms themselves live in
// boin_core.h so that they can be compiled and verified outside of R as well.

#include <Rcpp.h>
#include "boin_core.h"

using namespace Rcpp;

namespace simfastboin {

// Posterior probability Pr(p > target | y DLTs out of n patients) under a
// uniform Beta(1, 1) prior, evaluated with R's incomplete beta function.
double post_prob_above(double target, int y, int n) {
  return 1.0 - R::pbeta(target, y + 1.0, n - y + 1.0, 1, 0);
}

}  // namespace simfastboin

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnif {
  double operator()() { return unif_rand(); }
};

}  // namespace

// [[Rcpp::export]]
List boin_simulate_cpp(int n_trials,
                       NumericVector p_true,
                       IntegerVector cohort_size,
                       int start_dose,
                       int n_earlystop,
                       bool early_stop_simple,
                       bool titration,
                       bool extrasafe,
                       double target,
                       double cutoff_eli,
                       double offset,
                       IntegerVector b_esc,
                       IntegerVector b_deesc,
                       IntegerVector b_elim,
                       int max_total_pts) {

  const int n_doses = p_true.size();

  const std::vector<double> p(p_true.begin(), p_true.end());
  const std::vector<int> cs(cohort_size.begin(), cohort_size.end());
  const std::vector<int> be(b_esc.begin(), b_esc.end());
  const std::vector<int> bd(b_deesc.begin(), b_deesc.end());
  const std::vector<int> bl(b_elim.begin(), b_elim.end());

  IntegerMatrix n_pts(n_trials, n_doses);
  IntegerMatrix n_tox(n_trials, n_doses);
  LogicalMatrix eliminated(n_trials, n_doses);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nn(n_doses), yy(n_doses);
  std::vector<char> el(n_doses);
  RUnif rng;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    simfastboin::simulate_one(p, cs, start_dose, n_earlystop, early_stop_simple,
                              titration, extrasafe, target, cutoff_eli, offset,
                              be, bd, bl, max_total_pts, rng,
                              nn, yy, el, used, code);
    for (int j = 0; j < n_doses; ++j) {
      n_pts(t, j) = nn[j];
      n_tox(t, j) = yy[j];
      eliminated(t, j) = (el[j] == 1);
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
  }

  return List::create(
    _["n_pts"] = n_pts,
    _["n_tox"] = n_tox,
    _["eliminated"] = eliminated,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code
  );
}

// [[Rcpp::export]]
List boin_select_mtd_cpp(IntegerMatrix n_pts,
                         IntegerMatrix n_tox,
                         double target,
                         double cutoff_eli,
                         bool extrasafe,
                         double offset,
                         bool bound_mtd,
                         double lambda_d,
                         int min_mtd_sample) {

  const int n_trials = n_pts.nrow();
  const int n_doses = n_pts.ncol();

  IntegerVector mtd(n_trials);
  IntegerVector reason(n_trials);

  std::vector<int> nn(n_doses), yy(n_doses);

  for (int t = 0; t < n_trials; ++t) {
    for (int j = 0; j < n_doses; ++j) {
      nn[j] = n_pts(t, j);
      yy[j] = n_tox(t, j);
    }
    int code = 0;
    const int d = simfastboin::select_mtd_one(&nn[0], &yy[0], n_doses, target,
                                              cutoff_eli, extrasafe, offset,
                                              bound_mtd, lambda_d,
                                              min_mtd_sample, code);
    mtd[t] = (d == 0) ? NA_INTEGER : d;
    reason[t] = code;
  }

  return List::create(_["mtd"] = mtd, _["reason"] = reason);
}

// [[Rcpp::export]]
NumericMatrix boin_isotonic_cpp(IntegerMatrix n_pts,
                                IntegerMatrix n_tox,
                                LogicalMatrix admissible) {

  const int n_trials = n_pts.nrow();
  const int n_doses = n_pts.ncol();

  NumericMatrix out(n_trials, n_doses);
  std::fill(out.begin(), out.end(), NA_REAL);

  for (int t = 0; t < n_trials; ++t) {
    std::vector<int> idx;
    std::vector<double> phat, weight;
    for (int j = 0; j < n_doses; ++j) {
      if (admissible(t, j) == 1 && n_pts(t, j) > 0) {
        double p, w;
        simfastboin::phat_and_weight(n_tox(t, j), n_pts(t, j), p, w);
        idx.push_back(j);
        phat.push_back(p);
        weight.push_back(w);
      }
    }
    if (idx.empty()) continue;

    simfastboin::pava(phat, weight);
    for (std::size_t m = 0; m < idx.size(); ++m) out(t, idx[m]) = phat[m];
  }

  return out;
}
