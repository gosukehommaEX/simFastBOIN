#ifndef SIMFASTBOIN_BACKFILL_CORE_H
#define SIMFASTBOIN_BACKFILL_CORE_H

// Core algorithm of the BOIN designs with backfilling, written in plain C++
// like tite_core.h so that it can be verified outside of R.
//
// Two designs share the engine. BF-BOIN (Zhao et al., 2024) staggers the dose
// escalation (DE) cohorts completely and decides from the patients who have
// completed the DLT assessment only. BE-BOIN (Chen et al., 2026) takes the DE
// decisions of TITE-BOIN, with the pending patients imputed, and the same
// imputation enters every estimate used for backfilling.
//
// Random numbers come from three streams. The first supplies one uniform
// variate per DE patient, drawn in the cohort blocks of simulate_one() and
// simulate_tite_one(); the second, a generator of its own, the times between
// arrivals, one per arrival; the third, another generator of its own, the DLT
// variate of every backfill patient and the response variate of every patient.
// When no dose is ever opened for backfilling the first two streams are used
// exactly as in simulate_tite_one(), so that the trials agree with it.

#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include "tite_core.h"

namespace simfastboin {

// How the estimates treat the patients whose DLT assessment is pending.
enum BackfillEstimate {
  BACKFILL_OBSERVED = 0,    // BF-BOIN: completed patients only, staggered DE
  BACKFILL_IMPUTATION = 1   // BE-BOIN: TITE-BOIN single mean imputation
};

// What happens to a patient who arrives when no dose can take a patient.
enum NoSlot { NO_SLOT_WAIT = 0, NO_SLOT_LEAVE = 1 };

// Which conflicting backfill dose anchors the pooled estimate, and which of
// several open doses receives a backfill patient.
enum DoseChoice { CHOOSE_HIGHEST = 0, CHOOSE_LOWEST = 1 };

// Category of a decision, as in Table 2 of Zhao et al. (2024).
enum Category { CAT_ESCALATE = 0, CAT_STAY = 1, CAT_DEESCALATE = 2 };

// The third random number stream is seeded with the seed of the second, the
// one for the arrivals, combined with this key by exclusive or.
const std::uint64_t BACKFILL_STREAM_KEY = 0xD1B54A32D192ED03ULL;

// Standard normal distribution function.
inline double std_normal_cdf(double x) {
  return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// Standard normal quantile function, algorithm AS 241 (Wichura, 1988), which
// is also the algorithm of qnorm() in R.
inline double std_normal_quantile(double p) {
  if (p <= 0.0) return -std::numeric_limits<double>::infinity();
  if (p >= 1.0) return std::numeric_limits<double>::infinity();
  const double q = p - 0.5;
  double r, val;
  if (std::fabs(q) <= 0.425) {
    r = 0.180625 - q * q;
    val = q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r +
                    67265.770927008700853) * r + 45921.953931549871457) * r +
                  13731.693765509461125) * r + 1971.5909503065514427) * r +
                133.14166789178437745) * r + 3.387132872796366608) /
      (((((((r * 5226.495278852545925 + 28729.085735721942674) * r +
            39307.89580009271061) * r + 21213.794301586595867) * r +
          5394.1960214247511077) * r + 687.1870074920579083) * r +
        42.313330701600911252) * r + 1.0);
    return val;
  }
  r = (q < 0.0) ? p : 1.0 - p;
  r = std::sqrt(-std::log(r));
  if (r <= 5.0) {
    r -= 1.6;
    val = (((((((r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r +
                0.24178072517745061177) * r + 1.27045825245236838258) * r +
              3.64784832476320460504) * r + 5.7694972214606914055) * r +
            4.6303378461565452959) * r + 1.42343711074968357734) /
      (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r +
            0.0151986665636164571966) * r + 0.14810397642748007459) * r +
          0.68976733498510000455) * r + 1.6763848301838038494) * r +
        2.05319162663775882187) * r + 1.0);
  } else {
    r -= 5.0;
    val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r +
                0.0012426609473880784386) * r + 0.026532189526576123093) * r +
              0.29656057182850489123) * r + 1.7848265399172913358) * r +
            5.4637849111641143699) * r + 6.6579046435011037772) /
      (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r +
            1.8463183175100546818e-5) * r + 7.868691311456132591e-4) * r +
          0.0148753612908506148525) * r + 0.13692988092273580531) * r +
        0.59983220655588793769) * r + 1.0);
  }
  return (q < 0.0) ? -val : val;
}

// Response variate of a patient whose DLT variate is u, from an independent
// uniform v, under a Gaussian copula with correlation rho. With rho = 0 the
// response variate is v itself.
inline double response_variate(double u, double v, double rho) {
  if (rho == 0.0) return v;
  const double z = rho * std_normal_quantile(u) +
    std::sqrt(1.0 - rho * rho) * std_normal_quantile(v);
  return std_normal_cdf(z);
}

// The settings of a backfill trial that do not change from trial to trial.
struct BackfillDesign {
  int estimate;                    // BackfillEstimate
  std::vector<double> p_true;
  std::vector<double> p_resp;
  std::vector<int> cohort_size;
  int start_dose;                  // 1-based
  int n_earlystop;
  bool early_stop_simple;
  bool extrasafe;
  double target;
  double cutoff_eli;
  double offset;
  std::vector<int> b_esc;          // indexed by n - 1
  std::vector<int> b_deesc;
  std::vector<int> b_elim;         // 0 when there is no boundary
  int max_total_pts;               // patients in dose escalation
  double lambda_e;
  double lambda_d;
  double max_pending_ratio;        // imputation only
  int min_completed;               // imputation only
  double min_follow_up;            // imputation only
  double window;
  int accrual;
  double accrual_rate;
  std::vector<DltTimeModel> time_model;
  std::vector<DltTimeModel> resp_model;
  double resp_cor;
  int n_cap;
  int backfill_dose;               // DoseChoice
  int conflict_dose;               // DoseChoice
  int no_slot;                     // NoSlot
};

// Data at every dose at one moment.
struct DoseData {
  std::vector<int> n;          // treated
  std::vector<int> n_done;     // completed the DLT assessment
  std::vector<int> y;          // DLTs observed
  std::vector<int> pending;    // pending
  std::vector<double> stft;    // standardized total follow-up of the pending
  explicit DoseData(int k) : n(k), n_done(k), y(k), pending(k), stft(k) {}
};

// Boundary of a vector indexed by sample size, capped at its length.
inline int boundary_at(const std::vector<int>& b, int n) {
  int idx = n - 1;
  const int last = static_cast<int>(b.size()) - 1;
  if (idx > last) idx = last;
  return b[idx];
}

// The trial state and its rules.
class BackfillTrial {
 public:
  explicit BackfillTrial(const BackfillDesign& des)
    : des_(des), n_doses_(static_cast<int>(des.p_true.size())),
      data_(static_cast<int>(des.p_true.size())) {}

  template <class Rng, class Rng2, class Rng3>
  void run(Rng& rng, Rng2& rng2, Rng3& rng3,
           std::vector<int>& n_pts, std::vector<int>& n_tox,
           std::vector<int>& n_bf, std::vector<int>& n_resp,
           std::vector<char>& eliminated, int& cohorts_used, int& stop_code,
           double& duration, int& n_suspensions, double& time_suspended,
           int& n_turned_away);

 private:
  const BackfillDesign& des_;
  const int n_doses_;
  DoseData data_;

  std::vector<int> pt_dose_, pt_cohort_;
  std::vector<double> pt_entry_, pt_done_, pt_resp_;
  std::vector<char> pt_dlt_, pt_bf_;
  std::vector<int> bf_count_;

  bool imputation() const { return des_.estimate == BACKFILL_IMPUTATION; }

  // Data at time t.
  void collect(double t) {
    std::fill(data_.n.begin(), data_.n.end(), 0);
    std::fill(data_.n_done.begin(), data_.n_done.end(), 0);
    std::fill(data_.y.begin(), data_.y.end(), 0);
    std::fill(data_.pending.begin(), data_.pending.end(), 0);
    std::fill(data_.stft.begin(), data_.stft.end(), 0.0);
    const int m = static_cast<int>(pt_dose_.size());
    std::vector<double> follow_up(n_doses_, 0.0);
    for (int k = 0; k < m; ++k) {
      const int j = pt_dose_[k];
      ++data_.n[j];
      if (pt_done_[k] <= t) {
        ++data_.n_done[j];
        if (pt_dlt_[k]) ++data_.y[j];
      } else {
        ++data_.pending[j];
        follow_up[j] += t - pt_entry_[k];
      }
    }
    for (int j = 0; j < n_doses_; ++j) data_.stft[j] = follow_up[j] / des_.window;
  }

  // Numerator of the imputed estimate at dose j: the observed DLTs plus the
  // expected DLTs of the pending patients (Chen et al., 2026, equation 1).
  double imputed_dlts(int j) const {
    const double alpha = 0.5 * des_.target;
    const double p_post = (data_.y[j] + alpha) / (data_.n_done[j] + 1);
    return data_.y[j] + p_post / (1 - p_post) *
      (data_.pending[j] - data_.stft[j]);
  }

  // Category of y DLTs among n patients by the BOIN decision table.
  int category_counts(int y, int n) const {
    if (y <= boundary_at(des_.b_esc, n)) return CAT_ESCALATE;
    if (y >= boundary_at(des_.b_deesc, n)) return CAT_DEESCALATE;
    return CAT_STAY;
  }

  // Category of an estimated DLT rate.
  int category_rate(double q) const {
    if (q <= des_.lambda_e) return CAT_ESCALATE;
    if (q > des_.lambda_d) return CAT_DEESCALATE;
    return CAT_STAY;
  }

  // Category of the data pooled over doses lo to hi; returns -1 without data.
  int category_pooled(int lo, int hi) const {
    if (imputation()) {
      double num = 0.0;
      int den = 0;
      for (int j = lo; j <= hi; ++j) {
        if (data_.n[j] == 0) continue;
        num += imputed_dlts(j);
        den += data_.n[j];
      }
      if (den == 0) return -1;
      return category_rate(num / den);
    }
    int num = 0, den = 0;
    for (int j = lo; j <= hi; ++j) {
      num += data_.y[j];
      den += data_.n_done[j];
    }
    if (den == 0) return -1;
    return category_counts(num, den);
  }

  // Patients that count for elimination: every treated patient, the pending
  // ones as without DLT, for imputation; the completed ones otherwise.
  int n_for_elimination(int j) const {
    return imputation() ? data_.n[j] : data_.n_done[j];
  }

  // Dose for a backfill patient at time t, or -1 when no dose is open.
  int backfill_dose(int d, double t, const std::vector<char>& elim) const {
    if (d == 0) return -1;
    // A response at or below the dose is required.
    int lowest_response = n_doses_;
    const int m = static_cast<int>(pt_dose_.size());
    for (int k = 0; k < m; ++k) {
      if (pt_resp_[k] <= t && pt_dose_[k] < lowest_response) {
        lowest_response = pt_dose_[k];
      }
    }
    if (lowest_response >= d) return -1;
    // A dose is closed when its own estimate and the estimate pooled with the
    // next higher dose both exceed lambda_d; the doses above it close too.
    int closed_from = d;
    for (int b = 0; b < d; ++b) {
      if (data_.n[b] == 0) continue;
      const int own = imputation() ?
        category_rate(imputed_dlts(b) / data_.n[b]) :
        (data_.n_done[b] > 0 ? category_counts(data_.y[b], data_.n_done[b]) : -1);
      if (own != CAT_DEESCALATE) continue;
      if (category_pooled(b, b + 1) == CAT_DEESCALATE) {
        closed_from = b;
        break;
      }
    }
    const int upper = std::min(closed_from, d);
    if (des_.backfill_dose == CHOOSE_HIGHEST) {
      for (int b = upper - 1; b >= lowest_response; --b) {
        if (!elim[b] && data_.n[b] < des_.n_cap) return b;
      }
    } else {
      for (int b = lowest_response; b < upper; ++b) {
        if (!elim[b] && data_.n[b] < des_.n_cap) return b;
      }
    }
    return -1;
  }

  // Reconcile the decision at the current dose d with the backfilled doses
  // below it (Zhao et al., 2024, Table 2). Returns the new dose for a
  // de-escalation that the pooled data send below d - 1, or -1.
  int reconcile(int d, int& decision) const {
    const int cat_c = (decision == TITE_ESCALATE) ? CAT_ESCALATE :
      (decision == TITE_STAY) ? CAT_STAY : CAT_DEESCALATE;
    int b_star = -1;
    for (int b = 0; b < d; ++b) {
      if (bf_count_[b] == 0) continue;
      int cat_b;
      if (imputation()) {
        cat_b = category_rate(imputed_dlts(b) / data_.n[b]);
      } else {
        if (data_.n_done[b] == 0) continue;
        cat_b = category_counts(data_.y[b], data_.n_done[b]);
      }
      const bool conflict = (cat_b == CAT_DEESCALATE) ||
        (cat_b == CAT_STAY && cat_c == CAT_ESCALATE);
      if (!conflict) continue;
      if (des_.conflict_dose == CHOOSE_LOWEST) {
        if (b_star < 0) b_star = b;
      } else {
        b_star = b;
      }
    }
    if (b_star < 0) return -1;

    const int cat_q = category_pooled(b_star, d);
    if (cat_q == CAT_ESCALATE) { decision = TITE_ESCALATE; return -1; }
    if (cat_q == CAT_STAY) { decision = TITE_STAY; return -1; }
    decision = TITE_DEESCALATE;
    for (int k = d - 1; k >= b_star; --k) {
      if (category_pooled(b_star, k) != CAT_DEESCALATE) return k;
    }
    return (b_star > 0) ? b_star - 1 : 0;
  }

  template <class Rng3>
  void enroll(int dose, double t, double u, bool bf, int cohort, Rng3& rng3) {
    const double v = rng3();
    const double w = response_variate(u, v, des_.resp_cor);
    const bool dlt = u < des_.p_true[dose];
    pt_dose_.push_back(dose);
    pt_entry_.push_back(t);
    pt_done_.push_back(dlt ? t + dlt_time(des_.time_model[dose], u) :
                       t + des_.window);
    pt_dlt_.push_back(dlt ? 1 : 0);
    pt_resp_.push_back(w < des_.p_resp[dose] ?
                       t + dlt_time(des_.resp_model[dose], w) :
                       std::numeric_limits<double>::infinity());
    pt_bf_.push_back(bf ? 1 : 0);
    pt_cohort_.push_back(cohort);
    if (bf) ++bf_count_[dose];
  }
};

template <class Rng, class Rng2, class Rng3>
void BackfillTrial::run(Rng& rng, Rng2& rng2, Rng3& rng3,
                        std::vector<int>& n_pts, std::vector<int>& n_tox,
                        std::vector<int>& n_bf, std::vector<int>& n_resp,
                        std::vector<char>& eliminated, int& cohorts_used,
                        int& stop_code, double& duration, int& n_suspensions,
                        double& time_suspended, int& n_turned_away) {

  const double inf = std::numeric_limits<double>::infinity();
  const int n_cohort = static_cast<int>(des_.cohort_size.size());

  pt_dose_.clear(); pt_cohort_.clear(); pt_entry_.clear(); pt_done_.clear();
  pt_resp_.clear(); pt_dlt_.clear(); pt_bf_.clear();
  bf_count_.assign(n_doses_, 0);
  std::fill(eliminated.begin(), eliminated.end(), static_cast<char>(0));
  char* elim = &eliminated[0];

  int max_cs = 1;
  for (int i = 0; i < n_cohort; ++i) max_cs = std::max(max_cs, des_.cohort_size[i]);
  std::vector<double> cohort_u(max_cs, 0.0);

  int d = des_.start_dose - 1;
  int total = 0;
  double t = 0.0;
  stop_code = STOP_MAX_COHORTS;
  cohorts_used = 0;
  n_suspensions = 0;
  time_suspended = 0.0;
  n_turned_away = 0;

  for (int i = 0; i < n_cohort; ++i) {
    const int cs = des_.cohort_size[i];

    if (i > 0) {
      bool stop = false;
      bool suspended = false;

      for (;;) {
        collect(t);
        const int m = static_cast<int>(pt_dose_.size());

        // BF-BOIN waits for the whole of the last DE cohort before deciding
        // anything, as the BOIN design does.
        bool staggered_wait = false;
        if (!imputation()) {
          for (int k = 0; k < m; ++k) {
            if (pt_cohort_[k] == i - 1 && pt_done_[k] > t) {
              staggered_wait = true;
              break;
            }
          }
        }

        int decision = -1;
        int target = -1;
        if (staggered_wait) {
          decision = TITE_SUSPEND;
        } else {
          // Elimination at the current dose and the extra safety rule, nested
          // as in simulate_one().
          int nd = n_for_elimination(d);
          if (nd > 0) {
            const int b_el = boundary_at(des_.b_elim, nd);
            if (b_el > 0) {
              if (data_.y[d] >= b_el) {
                for (int j = d; j < n_doses_; ++j) elim[j] = 1;
                if (d == 0) { stop_code = STOP_LOWEST_ELIMINATED; stop = true; break; }
              }
              if (des_.extrasafe && d == 0 && nd >= 3) {
                if (post_prob_above(des_.target, data_.y[0], nd) >
                    des_.cutoff_eli - des_.offset) {
                  stop_code = STOP_LOWEST_TOO_TOXIC;
                  stop = true;
                  break;
                }
              }
            }
          }

          // Elimination at any dose.
          for (int j = 0; j < n_doses_; ++j) {
            const int nj = n_for_elimination(j);
            if (nj == 0) continue;
            const int b_el = boundary_at(des_.b_elim, nj);
            if (b_el > 0 && data_.y[j] >= b_el) {
              for (int k = j; k < n_doses_; ++k) elim[k] = 1;
              break;
            }
          }
          if (elim[0] == 1) { stop_code = STOP_LOWEST_ELIMINATED; stop = true; break; }

          if (elim[d] == 0) {
            if (imputation()) {
              double t_follow_up = -inf;
              for (int k = 0; k < m; ++k) {
                if (pt_dose_[k] == d && pt_done_[k] > t) {
                  const double t_k = pt_entry_[k] + des_.min_follow_up * des_.window;
                  if (t_k > t_follow_up) t_follow_up = t_k;
                }
              }
              decision = tite_decide(data_.n[d], data_.y[d], data_.pending[d],
                                     data_.stft[d], TITE_IMPUTATION, des_.target,
                                     des_.lambda_e, des_.lambda_d,
                                     des_.max_pending_ratio, des_.min_completed,
                                     t >= t_follow_up,
                                     boundary_at(des_.b_esc, data_.n[d]),
                                     boundary_at(des_.b_deesc, data_.n[d]));
            } else {
              const int cat = category_counts(data_.y[d], data_.n_done[d]);
              decision = (cat == CAT_ESCALATE) ? TITE_ESCALATE :
                (cat == CAT_STAY) ? TITE_STAY : TITE_DEESCALATE;
            }
            if (decision != TITE_SUSPEND) target = reconcile(d, decision);
          }

          // Early stopping on the patients at the current dose, backfilled
          // ones included.
          bool stop_now = false;
          if (data_.n[d] >= des_.n_earlystop) {
            if (des_.early_stop_simple) {
              stop_now = true;
            } else {
              const bool stay_here = (decision == TITE_STAY);
              const bool at_lowest = (d == 0) && (decision == TITE_DEESCALATE);
              const bool at_highest =
                ((d == n_doses_ - 1) || (elim[d + 1] == 1)) &&
                (decision == TITE_ESCALATE);
              stop_now = stay_here || at_lowest || at_highest;
            }
          }
          if (stop_now) { stop_code = STOP_N_EARLYSTOP; stop = true; break; }
        }

        if (decision == TITE_SUSPEND) {
          if (!suspended) { ++n_suspensions; suspended = true; }

          // The arriving patient is backfilled if a dose is open.
          const int b = backfill_dose(d, t, eliminated);
          if (b >= 0) {
            enroll(b, t, rng3(), true, -1, rng3);
            const double g = arrival_gap(des_.accrual, des_.accrual_rate, rng2);
            time_suspended += g;
            t += g;
            continue;
          }

          if (des_.no_slot == NO_SLOT_LEAVE) {
            ++n_turned_away;
            const double g = arrival_gap(des_.accrual, des_.accrual_rate, rng2);
            time_suspended += g;
            t += g;
            continue;
          }

          // Otherwise the patient waits for the next event that can change a
          // decision: a completed assessment at the current dose, the minimum
          // follow-up, a response, and, once a response has been observed,
          // a completed assessment at any dose.
          bool any_response = false;
          for (int k = 0; k < m; ++k) {
            if (pt_resp_[k] <= t) { any_response = true; break; }
          }
          double t_next = inf;
          double t_follow_up = -inf;
          for (int k = 0; k < m; ++k) {
            if (pt_done_[k] > t && (pt_dose_[k] == d || any_response) &&
                pt_done_[k] < t_next) {
              t_next = pt_done_[k];
            }
            if (pt_resp_[k] > t && pt_resp_[k] < t_next) t_next = pt_resp_[k];
            if (imputation() && pt_dose_[k] == d && pt_done_[k] > t) {
              const double t_k = pt_entry_[k] + des_.min_follow_up * des_.window;
              if (t_k > t_follow_up) t_follow_up = t_k;
            }
          }
          if (t_follow_up > t && t_follow_up < t_next) t_next = t_follow_up;
          time_suspended += t_next - t;
          t = t_next;
          continue;
        }

        // Dose transition.
        if (elim[d] == 1) {
          int first = 0;
          while (elim[first] == 0) ++first;
          d = first - 1;
        } else if (decision == TITE_ESCALATE) {
          if (d != n_doses_ - 1 && elim[d + 1] == 0) d = d + 1;
        } else if (decision == TITE_DEESCALATE) {
          if (target >= 0) {
            d = target;
          } else if (d != 0) {
            d = d - 1;
          }
        }
        break;
      }
      if (stop) break;
    }

    cohorts_used = i + 1;

    // One variate per patient for the whole cohort, drawn before the check on
    // the maximum sample size, as in simulate_one().
    for (int k = 0; k < cs; ++k) cohort_u[k] = rng();

    int n_enrol = cs;
    bool last = false;
    if (total + cs >= des_.max_total_pts) {
      n_enrol = des_.max_total_pts - total;
      if (n_enrol < 0) n_enrol = 0;
      if (n_enrol > cs) n_enrol = cs;
      stop_code = STOP_MAX_SAMPLE_SIZE;
      last = true;
    }

    for (int k = 0; k < n_enrol; ++k) {
      if (k > 0) t += arrival_gap(des_.accrual, des_.accrual_rate, rng2);
      enroll(d, t, cohort_u[k], false, i, rng3);
    }
    total += n_enrol;
    if (last) break;

    // Arrival of the next patient.
    t += arrival_gap(des_.accrual, des_.accrual_rate, rng2);
  }

  std::fill(n_pts.begin(), n_pts.end(), 0);
  std::fill(n_tox.begin(), n_tox.end(), 0);
  std::fill(n_bf.begin(), n_bf.end(), 0);
  std::fill(n_resp.begin(), n_resp.end(), 0);
  duration = 0.0;
  const int m = static_cast<int>(pt_dose_.size());
  for (int k = 0; k < m; ++k) {
    const int j = pt_dose_[k];
    ++n_pts[j];
    if (pt_dlt_[k]) ++n_tox[j];
    if (pt_bf_[k]) ++n_bf[j];
    if (pt_resp_[k] < inf) ++n_resp[j];
    if (pt_done_[k] > duration) duration = pt_done_[k];
  }
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_BACKFILL_CORE_H
