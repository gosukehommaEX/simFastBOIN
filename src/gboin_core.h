#ifndef SIMFASTBOIN_GBOIN_CORE_H
#define SIMFASTBOIN_GBOIN_CORE_H

// Core algorithm of the gBOIN design of Mu et al. (2019) for a toxicity
// outcome with several grades, each carrying a score that, divided by the
// largest score, lies between 0 and 1 (a quasi-Bernoulli outcome), written
// in plain C++ so that it can be compiled inside the package and inside a
// standalone verification harness.
//
// Floating point. The total normalized score of a dose is the total of the
// scores as given by the user, divided once by the largest score. With
// scores that are binary fractions, such as 0, 0.5, 1 and 1.5, the total of
// the scores as given is exact, so that the same data give the same total
// on every platform, also where the compiler fuses a multiplication and an
// addition into one operation (FMA). The gBOIN-ET core (gboinet_core.h)
// forms its totals with the same functions.
//
// The only external dependency is post_prob_above_real(), which the caller
// supplies.

#include <vector>
#include <algorithm>
#include <cmath>
#include "boin_core.h"

namespace simfastboin {

// Pr(mu > target | score total x among n patients) under a uniform
// Beta(1, 1) prior, that is one minus the Beta(x + 1, n - x + 1)
// distribution function at target. Defined by the caller.
double post_prob_above_real(double target, double x, int n);

// Total of the scores, as given by the user, of the patients of one dose,
// from the numbers of patients in each category. The sum runs over the
// categories in order, and every product is stored before it is added, so
// that it is not fused with the addition.
inline double raw_score_total(const int* count,
                              const std::vector<double>& score) {
  double x = 0.0;
  for (std::size_t k = 0; k < score.size(); ++k) {
    volatile double term = count[k] * score[k];
    x += term;
  }
  return x;
}

// Largest score, by which the scores are normalized.
inline double max_score(const std::vector<double>& score) {
  double m = score[0];
  for (std::size_t k = 1; k < score.size(); ++k) {
    if (score[k] > m) m = score[k];
  }
  return m;
}

// Total normalized score of the patients of one dose: the total of the
// scores as given, divided by the largest score.
inline double score_total(const int* count, const std::vector<double>& score,
                          double largest) {
  return raw_score_total(count, score) / largest;
}

// Simulate a single gBOIN trial.
//
// score holds the scores of the categories as given by the user, not
// normalized. cum_prob holds, for dose j and category k (0-based, the
// categories in the order of the scores given), the probability of category
// k or higher at index j * n_cat + k. rng() must return one uniform variate
// on (0, 1) per call. Exactly one variate is consumed per patient, drawn in
// enrollment order: the patient falls in the highest category k whose
// cum_prob exceeds the variate, and in category 0 when none does. A cohort
// cut short by the maximum sample size still draws a full cohort's worth.
// With two categories of scores 0 and 1 the trials are those of
// simulate_one() in boin_core.h without titration.
//
// The design escalates when the total normalized score x at the current
// dose is at most lambda_e * n and de-escalates when it is at least
// lambda_d * n. A dose with at least three patients and
// Pr(mu > target | x, n) above cutoff_eli is eliminated with every higher
// dose, and is always left for the dose below; the trial stops when the
// lowest dose is eliminated. With elim_by_count the dose is eliminated
// instead when x reaches b_elim, the elimination boundary of the number of
// DLTs of the BOIN design (indexed by the number of patients minus one, 0
// meaning that no boundary exists), as in the 'UnifiedDoseFinding' package.
//
// count receives the numbers of patients by dose and category (index
// j * n_cat + k), n_pts the patients by dose.
template <class Rng>
inline void simulate_one_gboin(const std::vector<double>& cum_prob,
                               const std::vector<double>& score,
                               const std::vector<int>& cohort_size,
                               int n_doses,
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
                               const std::vector<int>& b_elim,
                               int max_total_pts,
                               Rng& rng,
                               std::vector<int>& n_pts,
                               std::vector<int>& count,
                               std::vector<char>& eliminated,
                               int& cohorts_used,
                               int& stop_code) {

  const int n_cat = static_cast<int>(score.size());
  const int n_cohort = static_cast<int>(cohort_size.size());
  const double largest = max_score(score);

  std::fill(n_pts.begin(), n_pts.end(), 0);
  std::fill(count.begin(), count.end(), 0);
  std::fill(eliminated.begin(), eliminated.end(), static_cast<char>(0));

  int* n = &n_pts[0];
  char* elim = &eliminated[0];

  int d = start_dose - 1;
  int total = 0;
  stop_code = STOP_MAX_COHORTS;
  cohorts_used = 0;

  int max_cs = 1;
  for (int i = 0; i < n_cohort; ++i) max_cs = std::max(max_cs, cohort_size[i]);
  std::vector<int> cohort_cat(max_cs, 0);

  for (int i = 0; i < n_cohort; ++i) {
    cohorts_used = i + 1;
    const int cs = cohort_size[i];

    for (int k = 0; k < cs; ++k) {
      const double u = rng();
      int cat = 0;
      for (int c = n_cat - 1; c >= 1; --c) {
        if (u < cum_prob[d * n_cat + c]) { cat = c; break; }
      }
      cohort_cat[k] = cat;
    }

    bool last = false;
    int n_add = cs;
    if (total + cs >= max_total_pts) {
      n_add = max_total_pts - total;
      if (n_add < 0) n_add = 0;
      if (n_add > cs) n_add = cs;
      last = true;
    }
    for (int k = 0; k < n_add; ++k) ++count[d * n_cat + cohort_cat[k]];
    n[d] += n_add;
    total += n_add;
    if (last) { stop_code = STOP_MAX_SAMPLE_SIZE; break; }

    const double x = score_total(&count[d * n_cat], score, largest);

    // Whether some total score eliminates the dose at this sample size,
    // and whether the observed one does.
    bool has_boundary = false;
    bool gone = false;
    if (elim_by_count) {
      const int max_n = static_cast<int>(b_elim.size());
      const int idx = ((n[d] > max_n) ? max_n : n[d]) - 1;
      has_boundary = b_elim[idx] > 0;
      gone = has_boundary && x >= b_elim[idx];
    } else if (n[d] >= 3) {
      has_boundary =
        post_prob_above_real(target, static_cast<double>(n[d]), n[d]) >
          cutoff_eli;
      gone = post_prob_above_real(target, x, n[d]) > cutoff_eli;
    }
    const bool esc = (x <= lambda_e * n[d]) && !gone;
    const bool de = (x >= lambda_d * n[d]) || gone;

    // Dose elimination and the extra safety stopping rule. As in the BOIN
    // core, the safety rule is evaluated only at sample sizes at which some
    // total score eliminates the dose.
    if (has_boundary) {
      if (gone) {
        for (int j = d; j < n_doses; ++j) elim[j] = 1;
        if (d == 0) { stop_code = STOP_LOWEST_ELIMINATED; break; }
      }
      if (extrasafe && d == 0 && n[0] >= 3) {
        if (post_prob_above_real(target, x, n[0]) > cutoff_eli - offset) {
          stop_code = STOP_LOWEST_TOO_TOXIC;
          break;
        }
      }
    }

    // Early stopping on the number of patients at the current dose.
    bool stop_now = false;
    if (n[d] >= n_earlystop) {
      if (early_stop_simple) {
        stop_now = true;
      } else {
        const bool stay = !esc && !de;
        const bool at_lowest = (d == 0) && de;
        const bool at_highest =
          ((d == n_doses - 1) || (elim[d + 1] == 1)) && esc;
        stop_now = stay || at_lowest || at_highest;
      }
    }
    if (stop_now) { stop_code = STOP_N_EARLYSTOP; break; }

    // Dose transition.
    if (esc && d != n_doses - 1) {
      if (elim[d + 1] == 0) d = d + 1;
    } else if (de && d != 0) {
      d = d - 1;
    }
  }
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_GBOIN_CORE_H
