// Rcpp glue for the backfill engine, which lives in backfill_core.h.
// post_prob_above() is defined in boin.cpp.

#include <Rcpp.h>
#include "backfill_core.h"

using namespace Rcpp;

namespace {

// Draws one uniform variate from R's random number stream.
struct RUnifBackfill {
  double operator()() { return unif_rand(); }
};

}  // namespace

// [[Rcpp::export]]
List backfill_simulate_cpp(int n_trials,
                           int estimate,
                           NumericVector p_true,
                           NumericVector p_resp,
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
                           double resp_window,
                           double resp_late_fraction,
                           double resp_cor,
                           int n_cap,
                           int backfill_dose,
                           int conflict_dose,
                           int no_slot,
                           double stream_seed) {

  const int n_doses = p_true.size();

  simfastboin::BackfillDesign des;
  des.estimate = estimate;
  des.p_true.assign(p_true.begin(), p_true.end());
  des.p_resp.assign(p_resp.begin(), p_resp.end());
  des.cohort_size.assign(cohort_size.begin(), cohort_size.end());
  des.start_dose = start_dose;
  des.n_earlystop = n_earlystop;
  des.early_stop_simple = early_stop_simple;
  des.extrasafe = extrasafe;
  des.target = target;
  des.cutoff_eli = cutoff_eli;
  des.offset = offset;
  des.b_esc.assign(b_esc.begin(), b_esc.end());
  des.b_deesc.assign(b_deesc.begin(), b_deesc.end());
  des.b_elim.assign(b_elim.begin(), b_elim.end());
  des.max_total_pts = max_total_pts;
  des.lambda_e = lambda_e;
  des.lambda_d = lambda_d;
  des.max_pending_ratio = max_pending_ratio;
  des.min_completed = min_completed;
  des.min_follow_up = min_follow_up;
  des.window = window;
  des.accrual = accrual;
  des.accrual_rate = accrual_rate;
  des.time_model.resize(n_doses);
  des.resp_model.resize(n_doses);
  for (int j = 0; j < n_doses; ++j) {
    des.time_model[j] = simfastboin::make_dlt_time_model(
      dlt_time, des.p_true[j], window, late_fraction);
    des.resp_model[j] = simfastboin::make_dlt_time_model(
      dlt_time, des.p_resp[j], resp_window, resp_late_fraction);
  }
  des.resp_cor = resp_cor;
  des.n_cap = n_cap;
  des.backfill_dose = backfill_dose;
  des.conflict_dose = conflict_dose;
  des.no_slot = no_slot;

  IntegerMatrix n_pts(n_trials, n_doses);
  IntegerMatrix n_tox(n_trials, n_doses);
  IntegerMatrix n_bf(n_trials, n_doses);
  IntegerMatrix n_resp(n_trials, n_doses);
  LogicalMatrix eliminated(n_trials, n_doses);
  IntegerVector cohorts_used(n_trials);
  IntegerVector stop_code(n_trials);
  NumericVector duration(n_trials);
  IntegerVector n_suspensions(n_trials);
  NumericVector time_suspended(n_trials);
  IntegerVector n_turned_away(n_trials);

  std::vector<int> nn(n_doses), yy(n_doses), bb(n_doses), rr(n_doses);
  std::vector<char> el(n_doses);
  RUnifBackfill rng;
  const std::uint64_t seed64 =
    static_cast<std::uint64_t>(static_cast<std::int64_t>(stream_seed));
  simfastboin::Xoshiro256 rng2(seed64);
  simfastboin::Xoshiro256 rng3(seed64 ^ simfastboin::BACKFILL_STREAM_KEY);
  simfastboin::BackfillTrial trial(des);

  for (int t = 0; t < n_trials; ++t) {
    int used = 0, code = 0, susp = 0, away = 0;
    double dur = 0.0, susp_time = 0.0;
    trial.run(rng, rng2, rng3, nn, yy, bb, rr, el, used, code, dur, susp,
              susp_time, away);
    for (int j = 0; j < n_doses; ++j) {
      n_pts(t, j) = nn[j];
      n_tox(t, j) = yy[j];
      n_bf(t, j) = bb[j];
      n_resp(t, j) = rr[j];
      eliminated(t, j) = (el[j] == 1);
    }
    cohorts_used[t] = used;
    stop_code[t] = code;
    duration[t] = dur;
    n_suspensions[t] = susp;
    time_suspended[t] = susp_time;
    n_turned_away[t] = away;
  }

  return List::create(
    _["n_pts"] = n_pts,
    _["n_tox"] = n_tox,
    _["n_bf"] = n_bf,
    _["n_resp"] = n_resp,
    _["eliminated"] = eliminated,
    _["cohorts_used"] = cohorts_used,
    _["stop_code"] = stop_code,
    _["duration"] = duration,
    _["n_suspensions"] = n_suspensions,
    _["time_suspended"] = time_suspended,
    _["n_turned_away"] = n_turned_away
  );
}
