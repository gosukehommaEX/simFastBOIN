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

// Returns given values one after the other, for replaying a scripted trial.
struct ScriptedStream {
  const std::vector<double>& values;
  std::size_t next;
  const char* name;
  ScriptedStream(const std::vector<double>& v, const char* nm)
    : values(v), next(0), name(nm) {}
  double operator()() {
    if (next >= values.size()) {
      Rcpp::stop(std::string("the scripted stream '") + name + "' ran out");
    }
    return values[next++];
  }
};

// The design shared by every trial of a simulation.
simfastboin::BackfillDesign make_design(int estimate,
                                        const NumericVector& p_true,
                                        const NumericVector& p_resp,
                                        const IntegerVector& cohort_size,
                                        int start_dose, int n_earlystop,
                                        bool early_stop_simple, bool extrasafe,
                                        double target, double cutoff_eli,
                                        double offset,
                                        const IntegerVector& b_esc,
                                        const IntegerVector& b_deesc,
                                        const IntegerVector& b_elim,
                                        int max_total_pts, double lambda_e,
                                        double lambda_d,
                                        double max_pending_ratio,
                                        int min_completed, double min_follow_up,
                                        double window, int accrual,
                                        double accrual_rate, int dlt_time,
                                        double late_fraction,
                                        double resp_window,
                                        double resp_late_fraction,
                                        double resp_cor, int n_cap,
                                        int backfill_dose, int conflict_dose,
                                        int no_slot) {
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
  return des;
}

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

  const simfastboin::BackfillDesign des = make_design(
    estimate, p_true, p_resp, cohort_size, start_dose, n_earlystop,
    early_stop_simple, extrasafe, target, cutoff_eli, offset, b_esc, b_deesc,
    b_elim, max_total_pts, lambda_e, lambda_d, max_pending_ratio,
    min_completed, min_follow_up, window, accrual, accrual_rate, dlt_time,
    late_fraction, resp_window, resp_late_fraction, resp_cor, n_cap,
    backfill_dose, conflict_dose, no_slot);

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

// Replay one trial with scripted random numbers instead of random ones, for
// checking the engine against published trial examples. de_u holds the DLT
// variates of the escalation patients in cohort blocks, gap_u the variates of
// the times between arrivals, and extra_u the variates of the third stream
// (for each patient in enrollment order: the DLT variate if backfilled, then
// the response variate). Returns the patients in enrollment order.
// [[Rcpp::export]]
List backfill_replay_cpp(NumericVector de_u,
                         NumericVector gap_u,
                         NumericVector extra_u,
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
                         int no_slot) {

  const int n_doses = p_true.size();
  const simfastboin::BackfillDesign des = make_design(
    estimate, p_true, p_resp, cohort_size, start_dose, n_earlystop,
    early_stop_simple, extrasafe, target, cutoff_eli, offset, b_esc, b_deesc,
    b_elim, max_total_pts, lambda_e, lambda_d, max_pending_ratio,
    min_completed, min_follow_up, window, accrual, accrual_rate, dlt_time,
    late_fraction, resp_window, resp_late_fraction, resp_cor, n_cap,
    backfill_dose, conflict_dose, no_slot);

  const std::vector<double> v1(de_u.begin(), de_u.end());
  const std::vector<double> v2(gap_u.begin(), gap_u.end());
  const std::vector<double> v3(extra_u.begin(), extra_u.end());
  ScriptedStream rng(v1, "de_u");
  ScriptedStream rng2(v2, "gap_u");
  ScriptedStream rng3(v3, "extra_u");

  std::vector<int> nn(n_doses), yy(n_doses), bb(n_doses), rr(n_doses);
  std::vector<char> el(n_doses);
  int used = 0, code = 0, susp = 0, away = 0;
  double dur = 0.0, susp_time = 0.0;
  simfastboin::BackfillTrial trial(des);
  trial.run(rng, rng2, rng3, nn, yy, bb, rr, el, used, code, dur, susp,
            susp_time, away);

  const std::vector<int>& dose = trial.patient_dose();
  const int m = static_cast<int>(dose.size());
  IntegerVector out_dose(m);
  LogicalVector out_bf(m), out_dlt(m);
  NumericVector out_entry(m), out_done(m), out_resp(m);
  for (int k = 0; k < m; ++k) {
    out_dose[k] = dose[k] + 1;
    out_bf[k] = trial.patient_backfill()[k] == 1;
    out_dlt[k] = trial.patient_dlt()[k] == 1;
    out_entry[k] = trial.patient_entry()[k];
    out_done[k] = trial.patient_done()[k];
    out_resp[k] = trial.patient_response()[k];
  }

  return List::create(
    _["dose"] = out_dose,
    _["backfill"] = out_bf,
    _["dlt"] = out_dlt,
    _["entry"] = out_entry,
    _["done"] = out_done,
    _["response"] = out_resp,
    _["cohorts_used"] = used,
    _["stop_code"] = code,
    _["n_suspensions"] = susp,
    _["used"] = IntegerVector::create(static_cast<int>(rng.next),
                                      static_cast<int>(rng2.next),
                                      static_cast<int>(rng3.next))
  );
}
