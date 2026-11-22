// Rcpp glue for the core algorithms of the gBOIN-ET design, which live in
// gboinet_core.h so that they can be compiled and verified outside of R as
// well.

#include <Rcpp.h>
#include "gboinet_core.h"

using namespace Rcpp;

// beta_cdf(), which the core needs, is defined in comb.cpp with R's
// incomplete beta function.

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnifGet {
  double operator()() { return unif_rand(); }
};

// Copies a matrix with one row per category and one column per dose into
// the layout of the core, one block per dose.
std::vector<double> by_dose(NumericMatrix m) {
  const int n_cat = m.nrow();
  const int n_doses = m.ncol();
  std::vector<double> out(static_cast<std::size_t>(n_doses) * n_cat);
  for (int j = 0; j < n_doses; ++j) {
    for (int k = 0; k < n_cat; ++k) out[j * n_cat + k] = m(k, j);
  }
  return out;
}

}  // namespace

// [[Rcpp::export]]
List gboin_et_simulate_cpp(int n_trials,
                           NumericMatrix cum_tox,
                           NumericMatrix cum_resp,
                           NumericVector score_tox,
                           NumericVector score_resp,
                           double resp_cor,
                           IntegerVector cohort_size,
                           int start_dose,
                           int n_earlystop,
                           double lambda1,
                           double lambda2,
                           double eta1,
                           double p_eli,
                           double p_eli_resp,
                           double cutoff_eli,
                           double cutoff_eli_resp,
                           int boinet_compatible) {

  const int n_cat_t = cum_tox.nrow();
  const int n_cat_e = cum_resp.nrow();
  const int n_doses = cum_tox.ncol();

  simfastboin::GetDesign des;
  des.cum_tox = by_dose(cum_tox);
  des.cum_resp = by_dose(cum_resp);
  des.score_tox.assign(score_tox.begin(), score_tox.end());
  des.score_resp.assign(score_resp.begin(), score_resp.end());
  des.resp_cor = resp_cor;
  des.cohort_size.assign(cohort_size.begin(), cohort_size.end());
  des.start_dose = start_dose;
  des.n_stop = n_earlystop;
  des.lambda1 = lambda1;
  des.lambda2 = lambda2;
  des.eta1 = eta1;
  des.p_eli = p_eli;
  des.p_eli_resp = p_eli_resp;
  des.cut_tox = 1.0 - cutoff_eli;
  des.cut_resp = 1.0 - cutoff_eli_resp;
  des.compat = boinet_compatible;

  IntegerMatrix n_pts(n_trials, n_doses);
  NumericMatrix tox_sum(n_trials, n_doses);
  NumericMatrix resp_sum(n_trials, n_doses);
  IntegerVector n_tox_cat(n_trials * n_doses * n_cat_t);
  IntegerVector n_resp_cat(n_trials * n_doses * n_cat_e);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);

  std::vector<int> nn(n_doses), ct(n_doses * n_cat_t), ce(n_doses * n_cat_e);
  RUnifGet rng;

  for (int t = 0; t < n_trials; ++t) {
    int used = 0;
    int code = 0;
    simfastboin::simulate_get_one(des, rng, nn, ct, ce, used, code);
    for (int j = 0; j < n_doses; ++j) {
      n_pts(t, j) = nn[j];
      tox_sum(t, j) =
        simfastboin::score_total(&ct[j * n_cat_t], des.score_tox);
      resp_sum(t, j) =
        simfastboin::score_total(&ce[j * n_cat_e], des.score_resp);
      // Arrays of trials by doses by categories, first index fastest.
      for (int k = 0; k < n_cat_t; ++k) {
        n_tox_cat[t + n_trials * (j + n_doses * k)] = ct[j * n_cat_t + k];
      }
      for (int k = 0; k < n_cat_e; ++k) {
        n_resp_cat[t + n_trials * (j + n_doses * k)] = ce[j * n_cat_e + k];
      }
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
  }

  n_tox_cat.attr("dim") = IntegerVector::create(n_trials, n_doses, n_cat_t);
  n_resp_cat.attr("dim") = IntegerVector::create(n_trials, n_doses, n_cat_e);

  return List::create(
    _["n_pts"] = n_pts,
    _["tox_sum"] = tox_sum,
    _["resp_sum"] = resp_sum,
    _["n_tox_cat"] = n_tox_cat,
    _["n_resp_cat"] = n_resp_cat,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code
  );
}

// [[Rcpp::export]]
List gboin_et_select_obd_cpp(IntegerMatrix n_pts,
                             NumericMatrix tox_sum,
                             NumericMatrix resp_sum,
                             NumericMatrix probe,
                             LogicalVector early_stop,
                             double p_eli,
                             double p_eli_resp,
                             double cutoff_eli,
                             double cutoff_eli_resp,
                             int min_obd_sample,
                             int estimate,
                             int method,
                             double target,
                             double tox_upper,
                             double w1,
                             double w2,
                             double t_low,
                             double t_upp,
                             double e_low,
                             double e_upp,
                             double psi00,
                             double psi11,
                             int boinet_compatible) {

  const int n_trials = n_pts.nrow();
  const int n_doses = n_pts.ncol();

  simfastboin::EtObdRule rule;
  rule.method = method;
  rule.target = target;
  rule.tox_upper = tox_upper;
  rule.w1 = w1;
  rule.w2 = w2;
  rule.t_low = t_low;
  rule.t_upp = t_upp;
  rule.e_low = e_low;
  rule.e_upp = e_upp;
  rule.psi00 = psi00;
  rule.psi11 = psi11;

  IntegerVector obd(n_trials);
  IntegerVector reason(n_trials);
  NumericMatrix probe_out = clone(probe);
  std::vector<int> nn(n_doses);
  std::vector<double> xt(n_doses), xe(n_doses), pe(n_doses);

  for (int t = 0; t < n_trials; ++t) {
    for (int j = 0; j < n_doses; ++j) {
      nn[j] = n_pts(t, j);
      xt[j] = tox_sum(t, j);
      xe[j] = resp_sum(t, j);
      pe[j] = probe(t, j);
    }
    int why = 0;
    obd[t] = simfastboin::select_obd_get_one(
      &nn[0], &xt[0], &xe[0], &pe[0], n_doses, early_stop[t] == TRUE, p_eli,
      1.0 - cutoff_eli, p_eli_resp, 1.0 - cutoff_eli_resp,
      boinet_compatible != 0, min_obd_sample, estimate, rule, why);
    reason[t] = why;
    for (int j = 0; j < n_doses; ++j) probe_out(t, j) = pe[j];
  }

  return List::create(
    _["obd"] = obd,
    _["reason"] = reason,
    _["probe"] = probe_out
  );
}
