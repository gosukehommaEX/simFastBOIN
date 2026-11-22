#ifndef SIMFASTBOIN_SELECT_CORE_H
#define SIMFASTBOIN_SELECT_CORE_H

// Selection of the MTD from estimates that need not be whole numbers of
// DLTs among whole numbers of patients: the numbers used for the isotonic
// estimate are given apart from the patients treated. Used by the iBOIN
// design, whose estimate can add the prior information to the data, and by
// the gBOIN design, whose toxicity outcome is a score between 0 and 1.
//
// The arithmetic is that of phat_and_weight() and select_mtd_one() in
// boin_core.h, so that whole numbers give the same doses.

#include "boin_core.h"

namespace simfastboin {

// Posterior mean estimate and inverse-variance weight for y events, not
// necessarily a whole number, among n patients.
inline void phat_and_weight_real(double y, double n, double& phat,
                                 double& weight) {
  const double yy = y + 0.05;
  const double nn = n + 0.1;
  const double zz = n - y + 0.05;
  phat = yy / nn;
  weight = 1.0 / (yy * zz / (nn * nn * (nn + 1.0)));
}

// Select the MTD among the doses that are not eliminated and have been given
// to at least min_mtd_sample patients (and to at least one).
//
// elim marks the eliminated doses, n holds the patients treated, est_y and
// est_n the numbers of events and of patients that enter the isotonic
// estimate. Returns the selected dose (1-based) or 0 when no dose is
// selected.
inline int select_mtd_estimate(const char* elim, const int* n,
                               const double* est_y, const double* est_n,
                               int n_doses, double target,
                               int min_mtd_sample, int& reason) {

  if (elim[0] == 1) {
    reason = MTD_LOWEST_ELIMINATED;
    return 0;
  }

  std::vector<int> adm;
  for (int j = 0; j < n_doses; ++j) {
    if (elim[j] == 0 && n[j] > 0 && n[j] >= min_mtd_sample) adm.push_back(j);
  }
  if (adm.empty()) {
    reason = MTD_NO_ADMISSIBLE_DOSE;
    return 0;
  }

  const int k_adm = static_cast<int>(adm.size());
  std::vector<double> phat(k_adm), weight(k_adm);
  for (int k = 0; k < k_adm; ++k) {
    phat_and_weight_real(est_y[adm[k]], est_n[adm[k]], phat[k], weight[k]);
  }
  pava(phat, weight);
  for (int k = 0; k < k_adm; ++k) phat[k] += (k + 1) * 1e-10;

  int best = 0;
  double best_dist = std::fabs(phat[0] - target);
  for (int k = 1; k < k_adm; ++k) {
    const double dist = std::fabs(phat[k] - target);
    if (dist < best_dist) { best_dist = dist; best = k; }
  }

  reason = MTD_SELECTED;
  return adm[best] + 1;
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_SELECT_CORE_H
