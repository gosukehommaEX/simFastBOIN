// Rcpp glue for the iBOIN core algorithm, which lives in iboin_core.h so that
// it can be compiled and verified outside of R as well.

#include <Rcpp.h>
#include "iboin_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnifIboin {
  double operator()() { return unif_rand(); }
};

}  // namespace

// [[Rcpp::export]]
List iboin_simulate_cpp(int n_trials,
                        NumericVector p_true,
                        IntegerVector cohort_size,
                        int start_dose,
                        int n_earlystop,
                        bool early_stop_simple,
                        bool extrasafe,
                        double target,
                        double cutoff_eli,
                        double offset,
                        IntegerMatrix b_esc,
                        IntegerMatrix b_deesc,
                        IntegerVector b_elim,
                        int max_total_pts) {

  const int n_doses = p_true.size();
  const int max_n = b_elim.size();

  const std::vector<double> p(p_true.begin(), p_true.end());
  const std::vector<int> cs(cohort_size.begin(), cohort_size.end());
  const std::vector<int> bl(b_elim.begin(), b_elim.end());

  // One block of max_n boundaries per dose.
  std::vector<int> be(n_doses * max_n), bd(n_doses * max_n);
  for (int j = 0; j < n_doses; ++j) {
    for (int k = 0; k < max_n; ++k) {
      be[j * max_n + k] = b_esc(j, k);
      bd[j * max_n + k] = b_deesc(j, k);
    }
  }

  IntegerMatrix n_pts(n_trials, n_doses);
  IntegerMatrix n_tox(n_trials, n_doses);
  LogicalMatrix eliminated(n_trials, n_doses);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nn(n_doses), yy(n_doses);
  std::vector<char> el(n_doses);
  RUnifIboin rng;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    simfastboin::simulate_one_by_dose(p, cs, start_dose, n_earlystop,
                                      early_stop_simple, extrasafe, target,
                                      cutoff_eli, offset, be, bd, bl,
                                      max_total_pts, rng, nn, yy, el, used,
                                      code);
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
