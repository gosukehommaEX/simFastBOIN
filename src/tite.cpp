// Rcpp glue for the TITE-BOIN core algorithms, which live in tite_core.h.
// post_prob_above() is defined in boin.cpp.

#include <Rcpp.h>
#include "tite_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnifTite {
  double operator()() { return unif_rand(); }
};

}  // namespace

// [[Rcpp::export]]
List tite_boin_simulate_cpp(int n_trials,
                            NumericVector p_true,
                            IntegerVector cohort_size,
                            int start_dose,
                            int n_earlystop,
                            bool early_stop_simple,
                            bool extrasafe,
                            double target,
                            double cutoff_eli,
                            double offset,
                            IntegerVector b_esc,
                            IntegerVector b_deesc,
                            IntegerVector b_elim,
                            int max_total_pts,
                            int method,
                            double lambda_e,
                            double lambda_d,
                            double max_pending_ratio,
                            int min_completed,
                            double min_follow_up,
                            double window,
                            int accrual,
                            double accrual_rate,
                            int dlt_time,
                            double late_fraction,
                            bool weighted,
                            NumericVector prior_weights,
                            double stream_seed) {

  const int n_doses = p_true.size();

  const std::vector<double> p(p_true.begin(), p_true.end());
  const std::vector<int> cs(cohort_size.begin(), cohort_size.end());
  const std::vector<int> be(b_esc.begin(), b_esc.end());
  const std::vector<int> bd(b_deesc.begin(), b_deesc.end());
  const std::vector<int> bl(b_elim.begin(), b_elim.end());
  const std::vector<double> pw(prior_weights.begin(), prior_weights.end());

  std::vector<simfastboin::DltTimeModel> models(n_doses);
  for (int j = 0; j < n_doses; ++j) {
    models[j] = simfastboin::make_dlt_time_model(dlt_time, p[j], window,
                                                 late_fraction);
  }

  IntegerMatrix n_pts(n_trials, n_doses);
  IntegerMatrix n_tox(n_trials, n_doses);
  LogicalMatrix eliminated(n_trials, n_doses);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);
  NumericVector duration(n_trials);
  IntegerVector n_suspensions(n_trials);
  NumericVector time_suspended(n_trials);

  std::vector<int> nn(n_doses), yy(n_doses);
  std::vector<char> el(n_doses);
  RUnifTite rng;
  simfastboin::Xoshiro256 rng2(
    static_cast<std::uint64_t>(static_cast<std::int64_t>(stream_seed)));

  for (int t = 0; t < n_trials; ++t) {
    int used = 0, code = 0, susp = 0;
    double dur = 0.0, susp_time = 0.0;
    simfastboin::simulate_tite_one(p, cs, start_dose, n_earlystop,
                                   early_stop_simple, extrasafe, target,
                                   cutoff_eli, offset, be, bd, bl,
                                   max_total_pts, method, lambda_e, lambda_d,
                                   max_pending_ratio, min_completed,
                                   min_follow_up, window,
                                   accrual, accrual_rate, models, weighted, pw,
                                   rng, rng2,
                                   nn, yy, el, used, code, dur, susp, susp_time);
    for (int j = 0; j < n_doses; ++j) {
      n_pts(t, j) = nn[j];
      n_tox(t, j) = yy[j];
      eliminated(t, j) = (el[j] == 1);
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
    duration[t] = dur;
    n_suspensions[t] = susp;
    time_suspended[t] = susp_time;
  }

  return List::create(
    _["n_pts"] = n_pts,
    _["n_tox"] = n_tox,
    _["eliminated"] = eliminated,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code,
    _["duration"] = duration,
    _["n_suspensions"] = n_suspensions,
    _["time_suspended"] = time_suspended
  );
}

// Decision of the engine at given states, for checking it against
// tite_boin_decision_table(). Codes: 0 escalate, 1 stay, 2 de-escalate,
// 3 suspend, 4 de-escalate and eliminate. mf is the shortest follow-up of the
// pending patients as a fraction of the window, compared with min_follow_up.
// [[Rcpp::export]]
IntegerVector tite_decision_cpp(IntegerVector n,
                                IntegerVector n_tox,
                                IntegerVector n_pending,
                                NumericVector stft,
                                NumericVector mf,
                                int method,
                                double target,
                                double lambda_e,
                                double lambda_d,
                                double max_pending_ratio,
                                int min_completed,
                                double min_follow_up,
                                IntegerVector b_esc,
                                IntegerVector b_deesc,
                                IntegerVector b_elim) {

  const int m = n.size();
  IntegerVector out(m);
  for (int k = 0; k < m; ++k) {
    const int idx = n[k] - 1;
    if (b_elim[idx] > 0 && n_tox[k] >= b_elim[idx]) {
      out[k] = 4;
    } else {
      out[k] = simfastboin::tite_decide(n[k], n_tox[k], n_pending[k], stft[k],
                                        method, target, lambda_e, lambda_d,
                                        max_pending_ratio, min_completed,
                                        mf[k] >= min_follow_up,
                                        b_esc[idx], b_deesc[idx]);
    }
  }
  return out;
}
