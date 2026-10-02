#ifndef SIMFASTBOIN_TITE_CORE_H
#define SIMFASTBOIN_TITE_CORE_H

// Core algorithms of the time-to-event BOIN (TITE-BOIN) design, written in
// plain C++ like boin_core.h so that they can be verified outside of R.
//
// Random numbers come from two streams. The first supplies exactly one uniform
// variate per patient, drawn in enrollment order and in the same cohort blocks
// as simulate_one() in boin_core.h. That variate decides whether the patient
// has a DLT (u < p) and, if so, when it occurs. The second stream, a generator
// of its own, supplies the times between arrivals. When no patient is ever
// pending at a decision, the first stream is consumed exactly as in
// simulate_one() and the two engines produce the same trials.

#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include "boin_core.h"

namespace simfastboin {

// Decision codes of tite_decide().
enum TiteDecision {
  TITE_ESCALATE = 0,
  TITE_STAY = 1,
  TITE_DEESCALATE = 2,
  TITE_SUSPEND = 3
};

// Estimation methods.
enum TiteMethod {
  TITE_IMPUTATION = 0,
  TITE_ESS = 1
};

// Distributions of the time to DLT and of the times between arrivals.
enum DltTime { DLT_TIME_WEIBULL = 0, DLT_TIME_UNIFORM = 1 };
enum Accrual { ACCRUAL_EXPONENTIAL = 0, ACCRUAL_UNIFORM = 1, ACCRUAL_FIXED = 2 };

// xoshiro256** 1.0 (Blackman and Vigna), seeded through splitmix64. Returns
// uniform variates on [0, 1) with 53 random bits.
class Xoshiro256 {
 public:
  explicit Xoshiro256(std::uint64_t seed) {
    std::uint64_t x = seed;
    for (int i = 0; i < 4; ++i) s_[i] = splitmix64(x);
  }

  double operator()() {
    return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0);
  }

 private:
  std::uint64_t s_[4];

  static std::uint64_t splitmix64(std::uint64_t& x) {
    x += 0x9e3779b97f4a7c15ULL;
    std::uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }

  static std::uint64_t rotl(std::uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
  }

  std::uint64_t next() {
    const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
  }
};

// Decision at the current dose when no elimination applies.
//
// n patients have been treated, n_tox DLTs have been observed and n_pending
// patients are pending, whose standardized total follow-up time is stft. The
// arithmetic follows tite_boin_bounds() in the R code expression by expression,
// so that the engine and tite_boin_decision_table() agree. b_esc_n and
// b_deesc_n are the BOIN integer boundaries at n patients, used when nobody is
// pending. follow_up_ok tells whether every pending patient has been followed
// for at least the minimum required before an escalation (rule 2 of Chen et
// al., 2025); when it is false, an escalation becomes a suspension.
inline int tite_decide(int n, int n_tox, int n_pending, double stft,
                       int method, double target,
                       double lambda_e, double lambda_d,
                       double max_pending_ratio, int min_completed,
                       bool follow_up_ok,
                       int b_esc_n, int b_deesc_n) {

  if (n_pending == 0) {
    if (n_tox <= b_esc_n) return TITE_ESCALATE;
    if (n_tox >= b_deesc_n) return TITE_DEESCALATE;
    return TITE_STAY;
  }

  const double inf = std::numeric_limits<double>::infinity();
  const int n_done = n - n_pending;
  double esc, deesc, lower, upper, stat;

  if (method == TITE_IMPUTATION) {
    const double alpha = 0.5 * target;
    const double p_post = (n_tox + alpha) / (n_done + 1);
    const double odds = (1 - p_post) / p_post;
    const bool at_or_above =
      static_cast<double>(n_tox) / n >= target - 1e-12;
    esc = at_or_above ? inf : n_pending - odds * (n * lambda_e - n_tox);
    deesc = at_or_above ? n_pending - odds * (n * lambda_d - n_tox) : -inf;
    lower = 0.0;
    upper = static_cast<double>(n_pending);
    stat = stft;
  } else {
    esc = n_tox / lambda_e;
    deesc = (n_tox > 0) ? n_tox / lambda_d : -inf;
    lower = static_cast<double>(n_done);
    upper = static_cast<double>(n);
    stat = n_done + stft;
  }

  const int esc_action = (n_done >= min_completed && follow_up_ok) ?
    TITE_ESCALATE : TITE_SUSPEND;

  if (deesc >= upper) return TITE_DEESCALATE;
  if (static_cast<double>(n_pending) / n > max_pending_ratio) return TITE_SUSPEND;
  if (esc <= lower) return esc_action;
  if (esc < upper && stat >= esc) return esc_action;
  if (deesc > lower && stat <= deesc) return TITE_DEESCALATE;
  return TITE_STAY;
}

// Parameters of the time to DLT at one dose. For the Weibull distribution,
// Pr(X <= window) = p and a fraction late_fraction of the DLTs fall in the
// second half of the window, so Pr(X <= window / 2) = (1 - late_fraction) p.
struct DltTimeModel {
  int type;
  double window;
  double shape;
  double rate;
  double p;
};

inline DltTimeModel make_dlt_time_model(int type, double p, double window,
                                        double late_fraction) {
  DltTimeModel m;
  m.type = type;
  m.window = window;
  m.p = p;
  m.shape = 1.0;
  m.rate = 0.0;
  if (type == DLT_TIME_WEIBULL && p > 0.0 && p < 1.0) {
    const double p_half = (1.0 - late_fraction) * p;
    m.shape = std::log(std::log(1.0 - p) / std::log(1.0 - p_half)) / std::log(2.0);
    m.rate = -std::log(1.0 - p) / std::pow(window, m.shape);
  }
  return m;
}

// Time to DLT of a patient whose variate u is below p. Under the Weibull model
// this is the inverse distribution function at u, which lies inside the window
// exactly when u < p. Under the uniform model u / p is uniform on (0, 1)
// given u < p.
inline double dlt_time(const DltTimeModel& m, double u) {
  if (m.type == DLT_TIME_WEIBULL) {
    return std::pow(-std::log(1.0 - u) / m.rate, 1.0 / m.shape);
  }
  return m.window * (u / m.p);
}

// Weighted follow-up of a pending patient who has completed a fraction u of the
// window, under a piecewise uniform prior for the time to DLT that puts
// probabilities w[0], w[1] and w[2] on the three thirds of the window (Yuan et
// al., 2018, Supplementary Appendix D; Lin and Yuan, 2020, Supplementary S1).
// Each third contributes the share of it already covered. With equal weights
// the result is u itself.
inline double weighted_follow_up(double u, const std::vector<double>& w) {
  double s = 0.0;
  for (int k = 0; k < 3; ++k) {
    double part = 3.0 * u - k;
    if (part < 0.0) part = 0.0;
    if (part > 1.0) part = 1.0;
    s += w[k] * part;
  }
  return s;
}

// Time from one arrival to the next.
template <class Rng2>
inline double arrival_gap(int accrual, double rate, Rng2& rng2) {
  if (accrual == ACCRUAL_FIXED) return 1.0 / rate;
  const double v = rng2();
  if (accrual == ACCRUAL_UNIFORM) return v * 2.0 / rate;
  return -std::log(1.0 - v) / rate;
}

// Simulate a single TITE-BOIN trial.
//
// The dose for a cohort is decided when its first patient arrives, from the
// data available at that moment; the other patients of the cohort arrive later
// and receive the same dose. At every decision the rules are applied in the
// order of simulate_one(): elimination at the current dose and the extra safety
// rule, then elimination at any dose whose data have changed through late DLTs,
// then early stopping, then the dose transition. A decision to suspend accrual
// waits until the next pending patient at the current dose completes the
// assessment, or, when a minimum follow-up is required, until the last of them
// to arrive has been followed for that long, whichever comes first, and the
// rules are applied again at that moment.
//
// The duration of a trial is the time from the first arrival until every
// enrolled patient has completed the assessment.
template <class Rng, class Rng2>
inline void simulate_tite_one(const std::vector<double>& p_true,
                              const std::vector<int>& cohort_size,
                              int start_dose,
                              int n_earlystop,
                              bool early_stop_simple,
                              bool extrasafe,
                              double target,
                              double cutoff_eli,
                              double offset,
                              const std::vector<int>& b_esc,
                              const std::vector<int>& b_deesc,
                              const std::vector<int>& b_elim,
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
                              const std::vector<DltTimeModel>& time_model,
                              bool weighted,
                              const std::vector<double>& prior_weights,
                              Rng& rng,
                              Rng2& rng2,
                              std::vector<int>& n_pts,
                              std::vector<int>& n_tox,
                              std::vector<char>& eliminated,
                              int& cohorts_used,
                              int& stop_code,
                              double& duration,
                              int& n_suspensions,
                              double& time_suspended) {

  const int n_doses = static_cast<int>(p_true.size());
  const int n_cohort = static_cast<int>(cohort_size.size());
  const int max_n = static_cast<int>(b_esc.size());

  std::fill(n_pts.begin(), n_pts.end(), 0);
  std::fill(n_tox.begin(), n_tox.end(), 0);
  std::fill(eliminated.begin(), eliminated.end(), static_cast<char>(0));
  char* elim = &eliminated[0];

  // Patients enrolled so far: dose, arrival time, time at which the assessment
  // is complete, and whether a DLT occurs.
  std::vector<int> pt_dose;
  std::vector<double> pt_entry, pt_done;
  std::vector<char> pt_dlt;
  pt_dose.reserve(max_total_pts);
  pt_entry.reserve(max_total_pts);
  pt_done.reserve(max_total_pts);
  pt_dlt.reserve(max_total_pts);

  // Data observed at the time of a decision.
  std::vector<int> obs_tox(n_doses, 0);

  int max_cs = 1;
  for (int i = 0; i < n_cohort; ++i) max_cs = std::max(max_cs, cohort_size[i]);
  std::vector<double> cohort_u(max_cs, 0.0);

  int d = start_dose - 1;
  int total = 0;
  double t = 0.0;
  stop_code = STOP_MAX_COHORTS;
  cohorts_used = 0;
  n_suspensions = 0;
  time_suspended = 0.0;

  for (int i = 0; i < n_cohort; ++i) {
    const int cs = cohort_size[i];

    if (i > 0) {
      bool stop = false;
      bool suspended = false;

      for (;;) {
        // Observed DLTs at every dose, and the pending patients at the current
        // dose, at time t.
        std::fill(obs_tox.begin(), obs_tox.end(), 0);
        int n_pending = 0;
        double follow_up = 0.0;
        // Time at which every pending patient at the current dose will have
        // been followed for the minimum required before an escalation.
        double t_follow_up = -std::numeric_limits<double>::infinity();
        const int n_enrolled = static_cast<int>(pt_dose.size());
        for (int k = 0; k < n_enrolled; ++k) {
          if (pt_done[k] <= t) {
            if (pt_dlt[k]) ++obs_tox[pt_dose[k]];
          } else if (pt_dose[k] == d) {
            ++n_pending;
            const double t_k = pt_entry[k] + min_follow_up * window;
            if (t_k > t_follow_up) t_follow_up = t_k;
            if (weighted) {
              follow_up += weighted_follow_up((t - pt_entry[k]) / window,
                                              prior_weights);
            } else {
              follow_up += t - pt_entry[k];
            }
          }
        }
        // With equal prior weights STFT is computed exactly as before, so that
        // the default results do not depend on the weighting code.
        const double stft = weighted ? follow_up : follow_up / window;

        int nd = n_pts[d];
        if (nd > max_n) nd = max_n;
        const int idx = nd - 1;

        // Elimination at the current dose and the extra safety rule, nested as
        // in simulate_one().
        if (b_elim[idx] > 0) {
          if (obs_tox[d] >= b_elim[idx]) {
            for (int j = d; j < n_doses; ++j) elim[j] = 1;
            if (d == 0) { stop_code = STOP_LOWEST_ELIMINATED; stop = true; break; }
          }
          if (extrasafe && d == 0 && n_pts[0] >= 3) {
            if (post_prob_above(target, obs_tox[0], n_pts[0]) > cutoff_eli - offset) {
              stop_code = STOP_LOWEST_TOO_TOXIC;
              stop = true;
              break;
            }
          }
        }

        // Late DLTs can make any dose meet the elimination rule. Without
        // pending patients this finds nothing new.
        for (int j = 0; j < n_doses; ++j) {
          if (n_pts[j] == 0) continue;
          int nj = n_pts[j];
          if (nj > max_n) nj = max_n;
          if (b_elim[nj - 1] > 0 && obs_tox[j] >= b_elim[nj - 1]) {
            for (int k = j; k < n_doses; ++k) elim[k] = 1;
            break;
          }
        }
        if (elim[0] == 1) { stop_code = STOP_LOWEST_ELIMINATED; stop = true; break; }

        int decision = -1;
        if (elim[d] == 0) {
          decision = tite_decide(n_pts[d], obs_tox[d], n_pending, stft, method,
                                 target, lambda_e, lambda_d, max_pending_ratio,
                                 min_completed, t >= t_follow_up,
                                 b_esc[idx], b_deesc[idx]);
        }

        // Early stopping, evaluated at the current dose before the transition.
        bool stop_now = false;
        if (n_pts[d] >= n_earlystop) {
          if (early_stop_simple) {
            stop_now = true;
          } else {
            const bool stay_here = (decision == TITE_STAY);
            const bool at_lowest = (d == 0) && (decision == TITE_DEESCALATE);
            const bool at_highest =
              ((d == n_doses - 1) || (elim[d + 1] == 1)) &&
              (decision == TITE_ESCALATE);
            stop_now = stay_here || at_lowest || at_highest;
          }
        }
        if (stop_now) { stop_code = STOP_N_EARLYSTOP; stop = true; break; }

        if (decision == TITE_SUSPEND) {
          // Wait for the next pending patient at the current dose to complete,
          // or for the minimum follow-up to be reached if that comes first.
          // Without a minimum follow-up t_follow_up never lies ahead of t.
          double t_next = std::numeric_limits<double>::infinity();
          for (int k = 0; k < n_enrolled; ++k) {
            if (pt_dose[k] == d && pt_done[k] > t && pt_done[k] < t_next) {
              t_next = pt_done[k];
            }
          }
          if (t_follow_up > t && t_follow_up < t_next) t_next = t_follow_up;
          if (!suspended) { ++n_suspensions; suspended = true; }
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
          if (d != n_doses - 1 && elim[d + 1] == 0) d = d + 1;
        } else if (decision == TITE_DEESCALATE) {
          if (d != 0) d = d - 1;
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
    if (total + cs >= max_total_pts) {
      n_enrol = max_total_pts - total;
      if (n_enrol < 0) n_enrol = 0;
      if (n_enrol > cs) n_enrol = cs;
      stop_code = STOP_MAX_SAMPLE_SIZE;
      last = true;
    }

    for (int k = 0; k < n_enrol; ++k) {
      if (k > 0) t += arrival_gap(accrual, accrual_rate, rng2);
      const double u = cohort_u[k];
      const bool dlt = u < p_true[d];
      pt_dose.push_back(d);
      pt_entry.push_back(t);
      pt_done.push_back(dlt ? t + dlt_time(time_model[d], u) : t + window);
      pt_dlt.push_back(dlt ? 1 : 0);
      ++n_pts[d];
      if (dlt) ++n_tox[d];
    }
    total += n_enrol;
    if (last) break;

    // Arrival of the first patient of the next cohort.
    t += arrival_gap(accrual, accrual_rate, rng2);
  }

  duration = 0.0;
  for (std::size_t k = 0; k < pt_done.size(); ++k) {
    if (pt_done[k] > duration) duration = pt_done[k];
  }
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_TITE_CORE_H
