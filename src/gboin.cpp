// Rcpp glue for the gBOIN core algorithm, which lives in gboin_core.h so that
// it can be compiled and verified outside of R as well.

#include <Rcpp.h>
#include "gboin_core.h"

using namespace Rcpp;

namespace simfastboin {

// Posterior probability Pr(mu > target | score total x among n patients)
// under a uniform Beta(1, 1) prior, evaluated with R's incomplete beta
// function.
double post_prob_above_real(double target, double x, int n) {
  return 1.0 - R::pbeta(target, x + 1.0, n - x + 1.0, 1, 0);
}

}  // namespace simfastboin

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnifGboin {
  double operator()() { return unif_rand(); }
};

}  // namespace

// [[Rcpp::export]]
List gboin_simulate_cpp(int n_trials,
                        NumericMatrix cum_prob,
                        NumericVector score,
                        IntegerVector cohort_size,
                        int start_dose,
                        int n_earlystop,
                        bool early_stop_simple,
                        bool extrasafe,
                        double target,
                        double lambda_e,
                        double lambda_d,
                        double cutoff_eli,
                        double offset,
                        bool elim_by_count,
                        IntegerVector b_elim,
                        int max_total_pts) {

  // cum_prob has one row per category and one column per dose.
  const int n_cat = cum_prob.nrow();
  const int n_doses = cum_prob.ncol();

  std::vector<double> cp(n_doses * n_cat);
  for (int j = 0; j < n_doses; ++j) {
    for (int k = 0; k < n_cat; ++k) cp[j * n_cat + k] = cum_prob(k, j);
  }
  // The scores as given by the user; the core divides the totals by the
  // largest score.
  const std::vector<double> sc(score.begin(), score.end());
  const double largest = simfastboin::max_score(sc);
  const std::vector<int> cs(cohort_size.begin(), cohort_size.end());
  const std::vector<int> bl(b_elim.begin(), b_elim.end());

  IntegerMatrix n_pts(n_trials, n_doses);
  NumericMatrix score_sum(n_trials, n_doses);
  IntegerVector n_grade(n_trials * n_doses * n_cat);
  LogicalMatrix eliminated(n_trials, n_doses);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nn(n_doses), cnt(n_doses * n_cat);
  std::vector<char> el(n_doses);
  RUnifGboin rng;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    simfastboin::simulate_one_gboin(cp, sc, cs, n_doses, start_dose,
                                    n_earlystop, early_stop_simple,
                                    extrasafe, target, lambda_e, lambda_d,
                                    cutoff_eli, offset, elim_by_count, bl,
                                    max_total_pts, rng, nn, cnt, el, used,
                                    code);
    for (int j = 0; j < n_doses; ++j) {
      n_pts(t, j) = nn[j];
      score_sum(t, j) = simfastboin::score_total(&cnt[j * n_cat], sc,
                                                 largest);
      eliminated(t, j) = (el[j] == 1);
      // Array of trials by doses by categories, first index fastest.
      for (int k = 0; k < n_cat; ++k) {
        n_grade[t + n_trials * (j + n_doses * k)] = cnt[j * n_cat + k];
      }
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
  }

  n_grade.attr("dim") = IntegerVector::create(n_trials, n_doses, n_cat);

  return List::create(
    _["n_pts"] = n_pts,
    _["score_sum"] = score_sum,
    _["n_grade"] = n_grade,
    _["eliminated"] = eliminated,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code
  );
}
