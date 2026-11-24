#ifndef SIMFASTBOIN_WATERFALL_CORE_H
#define SIMFASTBOIN_WATERFALL_CORE_H

// Core algorithms of the waterfall design of Zhang and Yuan (2016), which
// finds the MTD contour of a drug combination by a sequence of
// one-dimensional BOIN subtrials. The rules are those of
// BOIN::get.oc.comb(mtd.contour = TRUE) in 'BOIN' 2.7.2. The code is plain
// C++ so that it can be compiled both inside the package (against R's
// numerical routines) and inside a standalone verification harness.
//
// The caller supplies post_prob_above() and beta_cdf().
//
// A dose matrix with J levels of drug A (rows) and K levels of drug B
// (columns) is stored column-major, as R stores a matrix: the combination in
// 0-based row a and column b is element a + J * b.
//
// Integer boundary vectors are indexed by the number of patients minus one.
// A value of 0 in b_elim means that no elimination boundary exists at that
// sample size, which corresponds to NA in the reference implementation.

#include <vector>
#include <algorithm>
#include <cmath>

namespace simfastboin {

// Pr(p > target | y DLTs out of n patients) under a uniform Beta(1, 1) prior.
// Defined by the caller.
double post_prob_above(double target, int y, int n);

// Distribution function of the beta distribution. Defined by the caller.
double beta_cdf(double x, double shape1, double shape2);

// Stop codes returned by waterfall_simulate_one().
enum WaterfallStopCode {
  WATERFALL_STOP_COMPLETED = 0,
  WATERFALL_STOP_FIRST_SUBTRIAL = 1,
  WATERFALL_STOP_LATER_SUBTRIAL = 2,
  WATERFALL_STOP_MAX_SAMPLE_SIZE = 3
};

// The design of a waterfall trial.
struct WaterfallDesign {
  int n_rows;                    // levels of drug A (J)
  int n_cols;                    // levels of drug B (K), with J <= K
  std::vector<double> p_true;    // true DLT probabilities, column-major
  std::vector<int> n_cohort;     // maximum number of cohorts of each subtrial
  int cohort_size;
  int n_earlystop;
  bool titration;
  bool extrasafe;
  bool bound_mtd;
  bool boin_compatible;
  double target;
  double cutoff_eli;
  double offset;
  double lambda_d;
  std::vector<int> b_esc;
  std::vector<int> b_deesc;
  std::vector<int> b_elim;
};

// Position in a boundary vector for n patients.
inline int waterfall_bound_index(const std::vector<int>& bound, int n) {
  const int max_n = static_cast<int>(bound.size());
  if (n > max_n) n = max_n;
  return n - 1;
}

// Isotonic regression by pooling the first adjacent violators, written as
// in the reference implementation so that the pooled values carry the same
// rounding: a pool takes the weighted mean of the two level sets it joins.
inline void waterfall_pava(std::vector<double>& x,
                           const std::vector<double>& wt) {
  const int m = static_cast<int>(x.size());
  if (m <= 1) return;
  std::vector<int> level(m);
  for (int i = 0; i < m; ++i) level[i] = i;
  for (;;) {
    int first = -1;
    for (int i = 0; i < m - 1; ++i) {
      if (x[i + 1] - x[i] < 0.0) {
        first = i;
        break;
      }
    }
    if (first < 0) break;
    const int level1 = level[first];
    const int level2 = level[first + 1];
    long double num = 0.0L;
    long double den = 0.0L;
    for (int i = 0; i < m; ++i) {
      if (level[i] == level1 || level[i] == level2) {
        const double product = x[i] * wt[i];
        num += product;
        den += wt[i];
      }
    }
    const double pooled = static_cast<double>(num) / static_cast<double>(den);
    for (int i = 0; i < m; ++i) {
      if (level[i] == level1 || level[i] == level2) {
        x[i] = pooled;
        level[i] = level1;
      }
    }
  }
}

// The candidate MTD of a completed subtrial from the patients n and the
// DLTs y of its doses, in increasing order. Returns the 0-based position of
// the candidate, or -1 when there is none. is_escalation is set when the
// data at the candidate would lead to an escalation.
inline int waterfall_candidate(const WaterfallDesign& des,
                               const std::vector<int>& n,
                               const std::vector<int>& y,
                               bool& is_escalation) {
  const int n_dose = static_cast<int>(n.size());
  is_escalation = false;

  // Elimination re-derived from the data of the subtrial: the first dose
  // with at least three patients that is too toxic, and every dose above.
  int first_eliminated = n_dose;
  for (int i = 0; i < n_dose; ++i) {
    if (n[i] >= 3 &&
        post_prob_above(des.target, y[i], n[i]) > des.cutoff_eli) {
      first_eliminated = i;
      break;
    }
  }
  if (des.extrasafe && n[0] >= 3 &&
      post_prob_above(des.target, y[0], n[0]) >
        des.cutoff_eli - des.offset) {
    first_eliminated = 0;
  }
  if (first_eliminated == 0) return -1;

  std::vector<int> adm;
  for (int i = 0; i < first_eliminated; ++i) {
    if (n[i] != 0) adm.push_back(i);
  }
  const int m = static_cast<int>(adm.size());
  if (m == 0) return -1;

  std::vector<double> phat(m), wt(m);
  for (int k = 0; k < m; ++k) {
    const double yi = y[adm[k]];
    const double ni = n[adm[k]];
    phat[k] = (yi + 0.05) / (ni + 0.1);
    const double var = (yi + 0.05) * (ni - yi + 0.05) /
      ((ni + 0.1) * (ni + 0.1) * (ni + 0.1 + 1.0));
    wt[k] = 1.0 / var;
  }
  waterfall_pava(phat, wt);
  for (int k = 0; k < m; ++k) {
    const double shift = (k + 1) * 1e-10;
    phat[k] = phat[k] + shift;
  }

  // With bound_mtd, only the estimates below the de-escalation boundary
  // remain; the position among them indexes the admissible doses, as in the
  // reference implementation.
  std::vector<double> kept;
  if (des.bound_mtd) {
    for (int k = 0; k < m; ++k) {
      if (phat[k] < des.lambda_d) kept.push_back(phat[k]);
    }
    if (kept.empty()) return -1;
  } else {
    kept = phat;
  }

  int best = 0;
  double best_dist = std::fabs(kept[0] - des.target);
  for (int k = 1; k < static_cast<int>(kept.size()); ++k) {
    const double dist = std::fabs(kept[k] - des.target);
    if (dist < best_dist) {
      best_dist = dist;
      best = k;
    }
  }
  const int sel = adm[best];
  const int ix = waterfall_bound_index(des.b_esc, n[sel]);
  is_escalation = (y[sel] <= des.b_esc[ix]);
  return sel;
}

// Run one subtrial on the doses whose cells are listed in dose, in
// increasing order, starting at position start, with at most n_coh cohorts.
//
// rng() must return one uniform variate on (0, 1) per call. With titration,
// one variate per dose of the subtrial is drawn first; afterwards one per
// patient, in enrollment order. n, y and elm receive the patients, the DLTs
// and the elimination flags of the doses. Returns the 0-based position of
// the candidate MTD, or -1 when the subtrial ends without one, in which case
// every dose is flagged when the lowest dose was found too toxic.
template <class Rng>
inline int waterfall_subtrial(const WaterfallDesign& des,
                              const std::vector<int>& dose,
                              int n_coh,
                              int start,
                              bool titration,
                              Rng& rng,
                              std::vector<int>& n,
                              std::vector<int>& y,
                              std::vector<char>& elm,
                              bool& is_escalation) {
  const int n_dose = static_cast<int>(dose.size());
  const int cs = des.cohort_size;
  n.assign(n_dose, 0);
  y.assign(n_dose, 0);
  elm.assign(n_dose, static_cast<char>(0));
  is_escalation = false;

  int d = start;
  bool early_stop = false;

  if (titration) {
    // One patient per dose, up to the first DLT or the highest dose. The
    // variates of every dose are drawn, as in the reference implementation.
    int first_dlt = -1;
    for (int i = 0; i < n_dose; ++i) {
      const bool dlt = rng() < des.p_true[dose[i]];
      if (dlt && first_dlt < 0) first_dlt = i;
    }
    if (first_dlt < 0) {
      d = n_dose - 1;
      for (int i = 0; i < n_dose; ++i) n[i] = 1;
    } else {
      d = first_dlt;
      for (int i = 0; i <= d; ++i) n[i] = 1;
      y[d] = 1;
    }
  }

  bool first_top_up = true;
  int total = 0;
  for (int i = 0; i < n_dose; ++i) total += n[i];

  for (int cohort = 0; cohort < n_coh; ++cohort) {
    int add = cs;
    if (titration && n[d] < cs && first_top_up) {
      // Fill the current dose up to the nominal cohort size, once.
      first_top_up = false;
      add = cs - 1;
    }
    int dlt = 0;
    const double p = des.p_true[dose[d]];
    for (int k = 0; k < add; ++k) {
      if (rng() < p) ++dlt;
    }
    y[d] += dlt;
    n[d] += add;
    total += add;

    // Elimination of the dose and the doses above, and the extra safety
    // rule, which is nested in the elimination branch and uses the
    // Beta(y + 0.05, n - y + 0.1) distribution, as in the reference
    // implementation.
    const int ix = waterfall_bound_index(des.b_elim, n[d]);
    if (des.b_elim[ix] > 0) {
      if (y[d] >= des.b_elim[ix]) {
        for (int i = d; i < n_dose; ++i) elm[i] = 1;
        if (d == 0) {
          early_stop = true;
          break;
        }
      }
      if (des.extrasafe && d == 0 && n[0] >= 3) {
        const double above =
          1.0 - beta_cdf(des.target, y[0] + 0.05, (n[0] - y[0]) + 0.1);
        if (above > des.cutoff_eli - des.offset) {
          early_stop = true;
          break;
        }
      }
    }

    if (y[d] <= des.b_esc[ix] && d != n_dose - 1) {
      if (elm[d + 1] == 0) d = d + 1;
    } else if (y[d] >= des.b_deesc[ix] && d != 0) {
      d = d - 1;
    }

    if (n[d] >= des.n_earlystop) break;
    if (total >= n_coh * cs) break;
  }

  if (early_stop) {
    for (int i = 0; i < n_dose; ++i) elm[i] = 1;
    return -1;
  }
  return waterfall_candidate(des, n, y, is_escalation);
}

// Simulate a single waterfall trial.
//
// The first subtrial runs up the first column (the lowest level of drug B)
// and then along the last row (the highest level of drug A). The following
// subtrials run along the rows below the row of the last candidate MTD, one
// row at a time, over the columns 2 to K, starting one column to the right
// of the last candidate. Three rules come from the reference implementation
// and not from the article: a candidate in the first column below the last
// row removes the rows above it and, when its data would lead to an
// escalation, is followed by a subtrial on the rest of its own row; the
// combinations to the right of a candidate in its row are excluded from the
// selection; and the trial ends when a subtrial ends without a candidate.
//
// When a subtrial other than the one added after the first ends without a
// candidate, the reference implementation drops its data. With
// boin_compatible the data are dropped in the same way; otherwise they are
// kept and every dose of that subtrial is flagged as eliminated.
// n_enrolled counts every patient treated, dropped or not.
template <class Rng>
inline void waterfall_simulate_one(const WaterfallDesign& des,
                                   Rng& rng,
                                   std::vector<int>& n_pts,
                                   std::vector<int>& n_tox,
                                   std::vector<char>& eliminated,
                                   int& n_subtrials,
                                   int& n_enrolled,
                                   int& stop_code) {
  const int J = des.n_rows;
  const int K = des.n_cols;
  const int cs = des.cohort_size;
  const int n_planned = static_cast<int>(des.n_cohort.size());

  std::fill(n_pts.begin(), n_pts.end(), 0);
  std::fill(n_tox.begin(), n_tox.end(), 0);
  std::fill(eliminated.begin(), eliminated.end(), static_cast<char>(0));
  n_subtrials = 0;
  n_enrolled = 0;
  stop_code = WATERFALL_STOP_MAX_SAMPLE_SIZE;

  int max_total = 0;
  for (int i = 0; i < n_planned; ++i) max_total += des.n_cohort[i] * cs;

  // The first subtrial: the first column below the last row, then the last
  // row.
  std::vector<int> dose;
  for (int a = 0; a < J - 1; ++a) dose.push_back(a);
  for (int b = 0; b < K; ++b) dose.push_back((J - 1) + J * b);

  std::vector<int> n, y;
  std::vector<char> elm;
  int start = 0;
  int sub = 0;
  int total = 0;

  while (total < max_total) {
    if (sub >= n_planned) break;
    bool is_escalation = false;
    const int sel = waterfall_subtrial(des, dose, des.n_cohort[sub], start,
                                       des.titration && sub == 0, rng, n, y,
                                       elm, is_escalation);
    ++n_subtrials;
    for (std::size_t i = 0; i < dose.size(); ++i) n_enrolled += n[i];

    if (sel < 0) {
      if (!des.boin_compatible) {
        for (std::size_t i = 0; i < dose.size(); ++i) {
          n_pts[dose[i]] = n[i];
          n_tox[dose[i]] = y[i];
          eliminated[dose[i]] = 1;
        }
      }
      stop_code = (sub == 0) ? WATERFALL_STOP_FIRST_SUBTRIAL :
        WATERFALL_STOP_LATER_SUBTRIAL;
      return;
    }

    for (std::size_t i = 0; i < dose.size(); ++i) {
      n_pts[dose[i]] = n[i];
      n_tox[dose[i]] = y[i];
      eliminated[dose[i]] = elm[i];
    }
    int di = dose[sel] % J;
    int dj = dose[sel] / J;

    if (sub == 0 && dj == 0 && di < J - 1) {
      // A candidate in the first column below the last row: the rows above
      // are removed, and the rest of its row is explored when its data
      // would lead to an escalation.
      for (int b = 0; b < K; ++b) {
        for (int a = di + 1; a < J; ++a) eliminated[a + J * b] = 1;
      }
      if (is_escalation && sub + 1 < n_planned) {
        ++sub;
        std::vector<int> dose_row;
        for (int b = 1; b < K; ++b) dose_row.push_back(di + J * b);
        bool row_escalation = false;
        const int sel_row = waterfall_subtrial(des, dose_row,
                                               des.n_cohort[sub], 0, false,
                                               rng, n, y, elm,
                                               row_escalation);
        ++n_subtrials;
        for (std::size_t i = 0; i < dose_row.size(); ++i) {
          n_enrolled += n[i];
          n_pts[dose_row[i]] = n[i];
          n_tox[dose_row[i]] = y[i];
          eliminated[dose_row[i]] = elm[i];
          if (sel_row < 0 && !des.boin_compatible) {
            eliminated[dose_row[i]] = 1;
          }
        }
        if (sel_row >= 0) {
          di = dose_row[sel_row] % J;
          dj = dose_row[sel_row] / J;
        }
      }
    }

    total = 0;
    for (int c = 0; c < J * K; ++c) total += n_pts[c];

    if (di == 0) {
      stop_code = WATERFALL_STOP_COMPLETED;
      return;
    }
    ++sub;
    for (int b = dj + 1; b < K; ++b) eliminated[di + J * b] = 1;

    // The next subtrial: the row below, columns 2 to K, starting one column
    // to the right of the candidate, or at the last column.
    dose.clear();
    for (int b = 1; b < K; ++b) dose.push_back((di - 1) + J * b);
    start = (dj == K - 1) ? dj - 1 : dj;
  }
}

// Select the MTD of every row of a single trial.
//
// iso holds the bivariate isotonic estimates, computed by the caller from
// (y + 0.05) / (n + 0.1), with 1.1 at the eliminated combinations, and the
// weights n + 0.1. 1e-5 * (row + column) (1-based) is added. From the
// highest level of drug A to the lowest, the treated combination of the row
// that is not eliminated and whose estimate is closest to the target is
// selected, the first among exact ties; a row whose first combination is
// eliminated, or without such a combination, selects none. A row whose MTD
// is not to the right of the MTD of the row above takes the column of that
// MTD. mtd_col receives the 1-based column of every row, or 0 for none.
inline void waterfall_choose_mtd(const double* iso, const int* n,
                                 const char* elim, int n_rows, int n_cols,
                                 double target, int* mtd_col) {
  const int J = n_rows;
  const int K = n_cols;
  for (int k = J - 1; k >= 0; --k) {
    int sel = 0;
    if (elim[k] != 1) {
      double best_dist = 0.0;
      for (int b = 0; b < K; ++b) {
        const int idx = k + J * b;
        if (n[idx] == 0 || elim[idx] == 1) continue;
        const double shift = 1e-5 * static_cast<double>((k + 1) + (b + 1));
        const double phat = iso[idx] + shift;
        const double dist = std::fabs(phat - target);
        if (sel == 0 || dist < best_dist) {
          best_dist = dist;
          sel = b + 1;
        }
      }
    }
    mtd_col[k] = sel;
    if (k < J - 1 && mtd_col[k] != 0 && mtd_col[k + 1] != 0 &&
        mtd_col[k] <= mtd_col[k + 1]) {
      mtd_col[k] = mtd_col[k + 1];
    }
  }
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_WATERFALL_CORE_H
