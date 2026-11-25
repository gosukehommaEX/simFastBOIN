// Rcpp glue for the U-BOIN design. The algorithms themselves live in
// uboin_core.h so that they can be compiled and verified outside of R as
// well.
//
// Trial data are arrays with dimensions (trial, dose, category), stored as R
// stores them.

#include <Rcpp.h>
#include "uboin_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct UbUnif {
  double operator()() { return unif_rand(); }
};

simfastboin::UbDesign ub_design(NumericMatrix p_outcome, NumericVector utility,
                                double prior_a, int max_n,
                                IntegerVector elim, IntegerVector toxic,
                                IntegerVector futile) {
  simfastboin::UbDesign des;
  des.n_doses = p_outcome.nrow();
  des.p_outcome.assign(4 * des.n_doses, 0.0);
  for (int j = 0; j < des.n_doses; ++j) {
    for (int k = 0; k < 4; ++k) des.p_outcome[4 * j + k] = p_outcome(j, k);
  }
  des.utility.assign(utility.begin(), utility.end());
  des.prior_a = prior_a;
  des.max_n = max_n;
  des.elim.assign(elim.begin(), elim.end());
  des.toxic.assign(toxic.begin(), toxic.end());
  des.futile.assign(futile.begin(), futile.end());
  des.n_cohort = 0;
  des.cohort_size = 0;
  des.start_dose = 0;
  des.n_stage1 = 0;
  des.n_stop = 0;
  return des;
}

}  // namespace

// [[Rcpp::export]]
List uboin_simulate_cpp(int n_trials,
                        NumericMatrix p_outcome,
                        int n_cohort,
                        int cohort_size,
                        int start_dose,
                        int n_stage1,
                        int n_stop,
                        NumericVector utility,
                        double prior_a,
                        int max_n,
                        IntegerVector elim,
                        IntegerVector toxic,
                        IntegerVector futile,
                        IntegerVector b_esc,
                        IntegerVector b_deesc) {

  simfastboin::UbDesign des = ub_design(p_outcome, utility, prior_a, max_n,
                                        elim, toxic, futile);
  des.n_cohort = n_cohort;
  des.cohort_size = cohort_size;
  des.start_dose = start_dose - 1;
  des.n_stage1 = n_stage1;
  des.n_stop = n_stop;
  des.b_esc.assign(b_esc.begin(), b_esc.end());
  des.b_deesc.assign(b_deesc.begin(), b_deesc.end());

  const int J = des.n_doses;
  IntegerVector n_outcome(n_trials * J * 4);
  LogicalVector eliminated(n_trials * J);
  IntegerVector cohorts_used(n_trials);
  LogicalVector stage2(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nk(4 * J);
  std::vector<char> el(J);
  UbUnif rng;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int second = 0;
    int code = 0;
    simfastboin::uboin_simulate_one(des, rng, nk, el, used, second, code);
    for (int j = 0; j < J; ++j) {
      for (int k = 0; k < 4; ++k) {
        n_outcome[t + n_trials * (j + J * k)] = nk[4 * j + k];
      }
      eliminated[t + n_trials * j] = (el[j] == 1);
    }
    cohorts_used[t] = used;
    stage2[t] = (second == 1);
    stop_code[t] = code;
  }

  n_outcome.attr("dim") = IntegerVector::create(n_trials, J, 4);
  eliminated.attr("dim") = IntegerVector::create(n_trials, J);

  return List::create(
    _["n_outcome"] = n_outcome,
    _["eliminated"] = eliminated,
    _["cohorts_used"] = cohorts_used,
    _["stage2"] = stage2,
    _["stop_code"] = stop_code
  );
}

// [[Rcpp::export]]
List uboin_select_cpp(IntegerVector n_outcome,
                      LogicalVector eliminated,
                      NumericVector utility,
                      double prior_a,
                      int max_n,
                      IntegerVector toxic,
                      IntegerVector futile) {

  const IntegerVector dims = n_outcome.attr("dim");
  const int n_trials = dims[0];
  const int J = dims[1];

  NumericMatrix dummy(J, 4);
  simfastboin::UbDesign des = ub_design(dummy, utility, prior_a, max_n,
                                        IntegerVector(0), toxic, futile);

  IntegerVector obd(n_trials);
  NumericMatrix util(n_trials, J);
  LogicalMatrix admissible(n_trials, J);
  std::vector<int> nk(4 * J);
  std::vector<char> el(J);

  for (int t = 0; t < n_trials; ++t) {
    for (int j = 0; j < J; ++j) {
      for (int k = 0; k < 4; ++k) {
        nk[4 * j + k] = n_outcome[t + n_trials * (j + J * k)];
      }
      el[j] = (eliminated[t + n_trials * j] == TRUE) ? 1 : 0;
    }
    for (int j = 0; j < J; ++j) {
      const int* c = &nk[4 * j];
      const int n = c[0] + c[1] + c[2] + c[3];
      util(t, j) = simfastboin::ub_utility(c, des);
      admissible(t, j) = n > 0 && el[j] == 0 &&
        simfastboin::ub_admissible(c, des);
    }
    const int best = simfastboin::ub_best_dose(&nk[0], &el[0], des);
    obd[t] = (best < 0) ? NA_INTEGER : best + 1;
  }

  return List::create(_["obd"] = obd, _["utility"] = util,
                      _["admissible"] = admissible);
}
