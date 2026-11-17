#ifndef SIMFASTBOIN_COMB_BOIN12_CORE_H
#define SIMFASTBOIN_COMB_BOIN12_CORE_H

// Core algorithms of the Comb-BOIN12 design (Lu, Zhang, Yuan and Lin, 2025)
// for two drugs given together, written in plain C++ like boin12_core.h so
// that they can be verified outside of R as well.
//
// A dose matrix with J levels of drug A (rows a) and K levels of drug B
// (columns b) is stored column-major: the combination (a, b), 0-based, is
// element c = a + J * b. The vectors of the B12Design (p_true, p_resp, q0,
// q1) are indexed in the same way, so that b12_patient() and
// b12_desirability() of boin12_core.h apply to a combination as to a dose.
//
// Rules for the next cohort at the current combination (j, k), with phat =
// y / n there (0 when untried):
//   - a combination is admissible when its DLT and response data both pass
//     (tables tox_ok and resp_ok; tox_ok_lowest at (0, 0), which holds the
//     tighter cutoff of extrasafe); a combination excluded by its DLT data
//     also excludes every combination (a, b) with a >= j and b >= k; with no
//     admissible combination the trial stops;
//   - phat >= lambda_d: the most desirable of the admissible combinations
//     (a, b) with a <= j and b <= k, other than (j, k), that have the
//     largest a + b; when there is none, the trial stops with
//     stop_no_candidate unless (j, k) = (0, 0), and otherwise stays at (j, k)
//     when it is admissible;
//   - lambda_e < phat < lambda_d and n >= n_star: the most desirable of
//     (j, k), (j - 1, k) and (j, k - 1); otherwise of (j, k) and its four
//     neighbors (j - 1, k), (j, k - 1), (j + 1, k), (j, k + 1); combinations
//     outside the matrix and inadmissible ones do not count;
//   - when none of these is admissible: the trial stops with
//     stop_no_candidate; otherwise the admissible combination nearest to
//     (j, k) in |a - j| + |b - k|, then with the lowest a + b, then the most
//     desirable;
//   - ties of desirability (exact equality): only the tied combinations with
//     the largest a + b remain, unless tie_random; when more than one
//     remains, one uniform variate u = rng() picks element floor(u m) of the
//     m remaining combinations sorted by a + b and then by a.
// With one row or one column, these are the rules of b12_next_dose() without
// exploration. The trial stops early when n_earlystop > 0, the current
// combination has n_earlystop patients and the next one is the same (or,
// with early_stop_simple, in any case).
//
// Selection at the end of a trial, from the final data and the bivariate
// isotonic estimates iso of the DLT rates (computed by the caller):
//   - the MTDC is the tried combination whose iso + 1e-5 (a + b + 2) is
//     closest to target, exact ties to the first in column-major order (the
//     rule of comb_choose_mtd(); below the target it favors the higher
//     combination);
//   - with the rectangle bound, the candidates are the combinations (a, b)
//     with a <= aT and b <= bT; with the contour bound, for every column b
//     the tried combination of that column chosen by the same rule bounds
//     the column, and the candidates are the combinations (a, b) with a at
//     most its row;
//   - the OBDC is the admissible tried candidate with the highest utility
//     estimate 100 pE (1 - pT) + psi00 (1 - pT)(1 - pE) + psi11 pT pE, with
//     pT = iso and pE the response estimate given by the caller; ties within
//     1e-12 go to the lowest a + b and then to the first in column-major
//     order.
//
// The response estimate of Lu et al. (2025) is the posterior mean under the
// logistic model logit(pE) = g0 + g1 a + g2 b + g3 a^2 + g4 b^2 with
// independent Cauchy(0, s) priors, computed by cb12_logistic_mean() with
// adaptive importance sampling (see there).

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include "boin12_core.h"

namespace simfastboin {

// Values of cb12_next() other than a combination.
const int CB12_NEXT_NONE = -1;   // no combination is admissible
const int CB12_NEXT_STOP = -2;   // stop_no_candidate stops the trial

// Candidate bounds of the OBDC selection.
enum Cb12Bound {
  CB12_BOUND_RECTANGLE = 0,
  CB12_BOUND_CONTOUR = 1
};

// The settings of a Comb-BOIN12 trial that do not change from trial to
// trial. b12.start_dose and b12.n_explore are not used. both_weight is
// psi00 + psi11 - 100, set to exactly 0 by the caller when the sum is 100
// up to rounding (see cb12_quasi_events()).
struct Cb12Design {
  B12Design b12;
  int n_rows;                      // J
  int n_cols;                      // K
  int start_row;                   // 0-based
  int start_col;                   // 0-based
  std::vector<int> tox_ok_lowest;  // tox_ok of the combination (0, 0)
  double both_weight;
};

// Quasi-binomial number of events of n patients with y DLTs, x responses and
// z patients with both: (100 x + psi00 (n - y - x) + (psi00 + psi11 - 100)
// z) / 100, the expression of b12_quasi_events() rearranged so that, when
// psi00 + psi11 is 100 (the utility of Lu et al., 2025), combinations with
// the same n, y and x have exactly the same value whatever z, and their tie
// is decided by the tie rule rather than by rounding errors.
inline double cb12_quasi_events(double n, double y, double x, double z,
                                const Cb12Design& des) {
  return (100.0 * x + des.b12.psi00 * (n - y - x) + des.both_weight * z) /
    100.0;
}

// Admissibility of every combination from the counts n, y (DLTs) and x
// (responses). Returns whether any combination is admissible.
inline bool cb12_admissible_set(const int* n, const int* y, const int* x,
                                const Cb12Design& des, char* adm) {
  const int J = des.n_rows;
  const int K = des.n_cols;
  const int w = des.b12.max_n + 1;
  bool any_adm = false;
  // blocked: the combination or one with both levels at most its own fails
  // the DLT rule.
  std::vector<char> blocked(J * K, 0);
  for (int b = 0; b < K; ++b) {
    for (int a = 0; a < J; ++a) {
      const int c = a + J * b;
      const std::vector<int>& tab = (c == 0) ? des.tox_ok_lowest :
        des.b12.tox_ok;
      bool blk = tab[n[c] * w + y[c]] != 1;
      if (a > 0 && blocked[c - 1]) blk = true;
      if (b > 0 && blocked[c - J]) blk = true;
      blocked[c] = blk ? 1 : 0;
      const bool ok = !blk && des.b12.resp_ok[n[c] * w + x[c]] == 1;
      adm[c] = ok ? 1 : 0;
      if (ok) any_adm = true;
    }
  }
  return any_adm;
}

// Desirability of combination c.
template <class Pb>
inline double cb12_desirability(int c, const int* n, const int* y,
                                const int* x, const int* z,
                                const Cb12Design& des, const Pb& pb) {
  const double q = cb12_quasi_events(n[c], y[c], x[c], z[c], des);
  return b12_desirability(n[c], q, des.b12, pb);
}

// The most desirable of the m combinations in cand, with the tie rule; rng()
// gives the variate of a tie that remains.
template <class Pb, class Rng>
inline int cb12_pick(const int* cand, int m, const int* n, const int* y,
                     const int* x, const int* z, const Cb12Design& des,
                     const Pb& pb, Rng& rng) {
  const int J = des.n_rows;
  std::vector<int> tied;
  double best = -1.0;
  for (int i = 0; i < m; ++i) {
    const double v = cb12_desirability(cand[i], n, y, x, z, des, pb);
    if (v > best) {
      best = v;
      tied.assign(1, cand[i]);
    } else if (v == best) {
      tied.push_back(cand[i]);
    }
  }
  if (!des.b12.tie_random) {
    int top = -1;
    for (std::size_t i = 0; i < tied.size(); ++i) {
      top = std::max(top, tied[i] % J + tied[i] / J);
    }
    std::vector<int> keep;
    for (std::size_t i = 0; i < tied.size(); ++i) {
      if (tied[i] % J + tied[i] / J == top) keep.push_back(tied[i]);
    }
    tied.swap(keep);
  }
  const int k = static_cast<int>(tied.size());
  if (k == 1) return tied[0];
  std::vector<int> order(tied);
  for (int i = 1; i < k; ++i) {
    const int c = order[i];
    const int s = c % J + c / J;
    int p = i - 1;
    while (p >= 0) {
      const int d = order[p];
      const int sd = d % J + d / J;
      if (sd > s || (sd == s && d % J > c % J)) {
        order[p + 1] = d;
        --p;
      } else {
        break;
      }
    }
    order[p + 1] = c;
  }
  int i = static_cast<int>(rng() * k);
  if (i >= k) i = k - 1;
  return order[i];
}

// Next combination (0-based index) from the current one cur, CB12_NEXT_NONE
// when no combination is admissible and CB12_NEXT_STOP when
// stop_no_candidate stops the trial. adm receives the admissibility of every
// combination.
template <class Pb, class Rng>
inline int cb12_next(int cur, const int* n, const int* y, const int* x,
                     const int* z, const Cb12Design& des, const Pb& pb,
                     Rng& rng, char* adm) {
  const B12Design& b12 = des.b12;
  const int J = des.n_rows;
  const int K = des.n_cols;
  if (!cb12_admissible_set(n, y, x, des, adm)) return CB12_NEXT_NONE;
  const int j = cur % J;
  const int k = cur / J;
  const double phat = n[cur] > 0 ? static_cast<double>(y[cur]) / n[cur] :
    0.0;
  std::vector<int> cand;

  if (phat >= b12.lambda_d) {
    int top = -1;
    for (int bb = 0; bb <= k; ++bb) {
      for (int aa = 0; aa <= j; ++aa) {
        const int c = aa + J * bb;
        if (c == cur || !adm[c]) continue;
        const int s = aa + bb;
        if (s > top) {
          top = s;
          cand.clear();
        }
        if (s == top) cand.push_back(c);
      }
    }
    if (!cand.empty()) {
      return cb12_pick(&cand[0], static_cast<int>(cand.size()), n, y, x, z,
                       des, pb, rng);
    }
    if (b12.stop_no_candidate && cur != 0) return CB12_NEXT_STOP;
    if (adm[cur]) return cur;
  } else {
    const bool pair = phat > b12.lambda_e && n[cur] >= b12.n_star;
    const int da[5] = {0, -1, 0, 1, 0};
    const int db[5] = {0, 0, -1, 0, 1};
    const int m = pair ? 3 : 5;
    for (int i = 0; i < m; ++i) {
      const int aa = j + da[i];
      const int bb = k + db[i];
      if (aa < 0 || aa >= J || bb < 0 || bb >= K) continue;
      const int c = aa + J * bb;
      if (adm[c]) cand.push_back(c);
    }
    if (!cand.empty()) {
      return cb12_pick(&cand[0], static_cast<int>(cand.size()), n, y, x, z,
                       des, pb, rng);
    }
  }

  if (b12.stop_no_candidate) return CB12_NEXT_STOP;
  int best_d = J + K;
  int best_s = J + K;
  cand.clear();
  for (int bb = 0; bb < K; ++bb) {
    for (int aa = 0; aa < J; ++aa) {
      const int c = aa + J * bb;
      if (!adm[c]) continue;
      const int d = std::abs(aa - j) + std::abs(bb - k);
      const int s = aa + bb;
      if (d < best_d || (d == best_d && s < best_s)) {
        best_d = d;
        best_s = s;
        cand.clear();
      }
      if (d == best_d && s == best_s) cand.push_back(c);
    }
  }
  return cb12_pick(&cand[0], static_cast<int>(cand.size()), n, y, x, z, des,
                   pb, rng);
}

// Simulate a single Comb-BOIN12 trial. nn, yy, xx and zz (J * K each)
// receive the patients, DLTs, responses and patients with both at every
// combination; stop_code is a B12Stop code.
//
// rng() must return one uniform variate on (0, 1) per call. Each patient
// draws two variates in enrollment order, the first for the DLT and the
// second for the response; a decision between tied combinations then draws
// one more variate when a tie remains.
template <class Rng, class Pb>
inline void cb12_simulate_one(const Cb12Design& des, Rng& rng, const Pb& pb,
                              std::vector<int>& nn, std::vector<int>& yy,
                              std::vector<int>& xx, std::vector<int>& zz,
                              int& cohorts_used, int& stop_code) {
  const B12Design& b12 = des.b12;
  std::fill(nn.begin(), nn.end(), 0);
  std::fill(yy.begin(), yy.end(), 0);
  std::fill(xx.begin(), xx.end(), 0);
  std::fill(zz.begin(), zz.end(), 0);
  std::vector<char> adm(nn.size());

  int cur = des.start_row + des.n_rows * des.start_col;
  cohorts_used = 0;
  stop_code = B12_STOP_MAX_COHORTS;

  for (std::size_t i = 0; i < b12.cohort_size.size(); ++i) {
    ++cohorts_used;
    for (int p = 0; p < b12.cohort_size[i]; ++p) {
      const double u = rng();
      const double v = rng();
      bool dlt, resp;
      b12_patient(b12, cur, u, v, dlt, resp);
      nn[cur] += 1;
      if (dlt) yy[cur] += 1;
      if (resp) xx[cur] += 1;
      if (dlt && resp) zz[cur] += 1;
    }

    const int nxt = cb12_next(cur, &nn[0], &yy[0], &xx[0], &zz[0], des, pb,
                              rng, &adm[0]);
    if (nxt == CB12_NEXT_NONE) {
      stop_code = B12_STOP_NO_ADMISSIBLE;
      return;
    }
    if (nxt == CB12_NEXT_STOP) {
      stop_code = B12_STOP_NO_CANDIDATE;
      return;
    }
    if (b12.n_earlystop > 0 && nn[cur] >= b12.n_earlystop &&
        (b12.early_stop_simple || nxt == cur)) {
      stop_code = B12_STOP_N_EARLYSTOP;
      return;
    }
    cur = nxt;
  }
}

// Index (0-based) of the tried combination, among those with
// in_line[c] != 0, whose iso + 1e-5 (a + b + 2) is closest to target, the
// first in column-major order among exact ties; -1 when there is none.
inline int cb12_closest(const double* iso, const int* n,
                        const std::vector<char>& in_line, int J, int K,
                        double target) {
  int pick = -1;
  double best = 0.0;
  for (int c = 0; c < J * K; ++c) {
    if (n[c] == 0 || !in_line[c]) continue;
    const double shift = 1e-5 * static_cast<double>(c % J + c / J + 2);
    const double dist = std::fabs(iso[c] + shift - target);
    if (pick < 0 || dist < best) {
      best = dist;
      pick = c;
    }
  }
  return pick;
}

// Select the OBDC of one trial. Returns its index (0-based) or -1, sets
// reason (EtObdReason) and mtdc (-1 when there is no selection), and writes
// the utility estimates (0 to 100) of the tried combinations to util (the
// others are left untouched).
inline int cb12_select_obdc_one(const int* n, const int* y, const int* x,
                                const double* iso, const double* pe,
                                bool early_stop, const Cb12Design& des,
                                double target, int bound, int& reason,
                                int& mtdc, double* util) {
  const int J = des.n_rows;
  const int K = des.n_cols;
  mtdc = -1;
  if (early_stop) {
    reason = ET_OBD_EARLY_STOP;
    return -1;
  }
  std::vector<char> adm(J * K);
  if (!cb12_admissible_set(n, y, x, des, &adm[0])) {
    reason = ET_OBD_NO_ADMISSIBLE;
    return -1;
  }
  std::vector<char> all(J * K, 1);
  mtdc = cb12_closest(iso, n, all, J, K, target);

  // Highest level of drug A allowed in every column.
  std::vector<int> top_row(K, -1);
  if (bound == CB12_BOUND_CONTOUR) {
    for (int b = 0; b < K; ++b) {
      std::vector<char> col(J * K, 0);
      for (int a = 0; a < J; ++a) col[a + J * b] = 1;
      const int c = cb12_closest(iso, n, col, J, K, target);
      if (c >= 0) top_row[b] = c % J;
    }
  } else {
    for (int b = 0; b <= mtdc / J; ++b) top_row[b] = mtdc % J;
  }

  int obd = -1;
  double best = 0.0;
  int best_s = 0;
  for (int c = 0; c < J * K; ++c) {
    if (n[c] == 0) continue;
    const double pt = iso[c];
    const double pr = pe[c];
    const double value = 100.0 * pr * (1.0 - pt) +
      des.b12.psi00 * (1.0 - pt) * (1.0 - pr) + des.b12.psi11 * pt * pr;
    util[c] = value;
    const int a = c % J;
    const int b = c / J;
    if (!adm[c] || a > top_row[b]) continue;
    const int s = a + b;
    if (obd < 0 || value > best + 1e-12 ||
        (value >= best - 1e-12 && s < best_s)) {
      obd = c;
      best = value;
      best_s = s;
    }
  }
  if (obd < 0) {
    reason = ET_OBD_NO_ADMISSIBLE;
    return -1;
  }
  reason = ET_OBD_SELECTED;
  return obd;
}

// The logistic model of the response probabilities and the points of the
// importance sampling. Element c of fa and fb is the standardized level of
// drug A and of drug B of combination c; z holds n_draws rows of five
// standard normal variates (row-major) and w n_draws chi-square variates with
// df degrees of freedom, both obtained by the caller from Halton points.
struct Cb12Logistic {
  std::vector<double> fa;
  std::vector<double> fb;
  double scale0;                   // Cauchy scale of the intercept
  double scale1;                   // Cauchy scale of the other coefficients
  std::vector<double> z;
  std::vector<double> w;
  double df;
  int rounds;                      // refits of the proposal
};

// log(1 + exp(e)) without overflow.
inline double cb12_log1pexp(double e) {
  return e > 0.0 ? e + std::log1p(std::exp(-e)) : std::log1p(std::exp(e));
}

// Cholesky factor (lower, row-major 5 x 5) of a symmetric matrix; returns
// false when it is not positive definite.
inline bool cb12_cholesky(const double* s, double* l) {
  for (int i = 0; i < 25; ++i) l[i] = 0.0;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = s[i * 5 + j];
      for (int k = 0; k < j; ++k) sum -= l[i * 5 + k] * l[j * 5 + k];
      if (i == j) {
        if (!(sum > 0.0)) return false;
        l[i * 5 + i] = std::sqrt(sum);
      } else {
        l[i * 5 + j] = sum / l[j * 5 + j];
      }
    }
  }
  return true;
}

// Solve (L L') v = r for v, L from cb12_cholesky().
inline void cb12_chol_solve(const double* l, const double* r, double* v) {
  double t[5];
  for (int i = 0; i < 5; ++i) {
    double sum = r[i];
    for (int k = 0; k < i; ++k) sum -= l[i * 5 + k] * t[k];
    t[i] = sum / l[i * 5 + i];
  }
  for (int i = 4; i >= 0; --i) {
    double sum = t[i];
    for (int k = i + 1; k < 5; ++k) sum -= l[k * 5 + i] * v[k];
    v[i] = sum / l[i * 5 + i];
  }
}

// Inverse of a symmetric positive definite matrix from its Cholesky factor.
inline void cb12_chol_inverse(const double* l, double* inv) {
  for (int j = 0; j < 5; ++j) {
    double e[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double col[5];
    e[j] = 1.0;
    cb12_chol_solve(l, e, col);
    for (int i = 0; i < 5; ++i) inv[i * 5 + j] = col[i];
  }
}

// Posterior mean of the response probability of every tried combination
// (pe, J * K; NaN for the untried ones) from the patients n and responses x
// of the tried combinations.
//
// The posterior mode is found by Newton steps with a positive definite
// curvature (the Fisher information of the data plus 2 / (s^2 + g^2) for
// each Cauchy prior) and step halving. The proposal is a multivariate t
// with df degrees of freedom centered at the mode, with the inverse of the
// negative Hessian there as scale matrix (the positive definite curvature
// when the Hessian is not). The points are mode + L z / sqrt(w / df). The
// proposal is then refitted rounds times to the weighted mean and covariance
// (times (df - 2) / df) of the weighted points, and the weights of the last
// proposal give the estimates. No random variate is drawn.
inline void cb12_logistic_mean(const int* n, const int* x, int n_cells,
                               const Cb12Logistic& lg, double* pe) {
  // Design rows of the tried combinations.
  std::vector<double> f;
  std::vector<double> nv, xv;
  for (int c = 0; c < n_cells; ++c) {
    if (n[c] == 0) continue;
    const double r[5] = {1.0, lg.fa[c], lg.fb[c], lg.fa[c] * lg.fa[c],
                         lg.fb[c] * lg.fb[c]};
    f.insert(f.end(), r, r + 5);
    nv.push_back(n[c]);
    xv.push_back(x[c]);
  }
  const int m = static_cast<int>(nv.size());
  const double s2[5] = {lg.scale0 * lg.scale0, lg.scale1 * lg.scale1,
                        lg.scale1 * lg.scale1, lg.scale1 * lg.scale1,
                        lg.scale1 * lg.scale1};

  // Negative log posterior (up to a constant).
  struct Obj {
    static double value(const double* th, const std::vector<double>& f,
                        const std::vector<double>& nv,
                        const std::vector<double>& xv, const double* s2,
                        int m) {
      double out = 0.0;
      for (int i = 0; i < m; ++i) {
        double eta = 0.0;
        for (int k = 0; k < 5; ++k) eta += f[i * 5 + k] * th[k];
        out -= xv[i] * eta - nv[i] * cb12_log1pexp(eta);
      }
      for (int k = 0; k < 5; ++k) out += std::log1p(th[k] * th[k] / s2[k]);
      return out;
    }
  };

  double th[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double cur = Obj::value(th, f, nv, xv, s2, m);
  for (int it = 0; it < 200; ++it) {
    double g[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double h[25];
    for (int i = 0; i < 25; ++i) h[i] = 0.0;
    for (int i = 0; i < m; ++i) {
      double eta = 0.0;
      for (int k = 0; k < 5; ++k) eta += f[i * 5 + k] * th[k];
      const double p = 1.0 / (1.0 + std::exp(-eta));
      const double r = nv[i] * p - xv[i];
      const double wgt = nv[i] * p * (1.0 - p);
      for (int a = 0; a < 5; ++a) {
        g[a] += r * f[i * 5 + a];
        for (int b = 0; b < 5; ++b) {
          h[a * 5 + b] += wgt * f[i * 5 + a] * f[i * 5 + b];
        }
      }
    }
    for (int k = 0; k < 5; ++k) {
      g[k] += 2.0 * th[k] / (s2[k] + th[k] * th[k]);
      h[k * 5 + k] += 2.0 / (s2[k] + th[k] * th[k]);
    }
    double l[25];
    if (!cb12_cholesky(h, l)) break;
    double step[5];
    cb12_chol_solve(l, g, step);
    double gnorm = 0.0;
    for (int k = 0; k < 5; ++k) gnorm += step[k] * g[k];
    if (gnorm < 1e-20) break;
    double t = 1.0;
    double nxt[5];
    double val = cur;
    bool moved = false;
    for (int half = 0; half < 60; ++half) {
      for (int k = 0; k < 5; ++k) nxt[k] = th[k] - t * step[k];
      val = Obj::value(nxt, f, nv, xv, s2, m);
      if (val <= cur) {
        moved = true;
        break;
      }
      t *= 0.5;
    }
    if (!moved) break;
    const double gain = cur - val;
    for (int k = 0; k < 5; ++k) th[k] = nxt[k];
    cur = val;
    if (gain < 1e-14 * (1.0 + std::fabs(cur))) break;
  }

  // Scale matrix: inverse of the negative Hessian at the mode, or of the
  // positive definite curvature when the Hessian is not positive definite.
  double hess[25], curv[25];
  for (int i = 0; i < 25; ++i) {
    hess[i] = 0.0;
    curv[i] = 0.0;
  }
  for (int i = 0; i < m; ++i) {
    double eta = 0.0;
    for (int k = 0; k < 5; ++k) eta += f[i * 5 + k] * th[k];
    const double p = 1.0 / (1.0 + std::exp(-eta));
    const double wgt = nv[i] * p * (1.0 - p);
    for (int a = 0; a < 5; ++a) {
      for (int b = 0; b < 5; ++b) {
        hess[a * 5 + b] += wgt * f[i * 5 + a] * f[i * 5 + b];
      }
    }
  }
  for (int i = 0; i < 25; ++i) curv[i] = hess[i];
  for (int k = 0; k < 5; ++k) {
    const double d = s2[k] + th[k] * th[k];
    hess[k * 5 + k] += 2.0 * (s2[k] - th[k] * th[k]) / (d * d);
    curv[k * 5 + k] += 2.0 / d;
  }
  double l[25], inv[25], scale[25];
  if (!cb12_cholesky(hess, l)) cb12_cholesky(curv, l);
  cb12_chol_inverse(l, inv);
  for (int i = 0; i < 25; ++i) scale[i] = inv[i];

  double mu[5];
  for (int k = 0; k < 5; ++k) mu[k] = th[k];
  const int n_draws = static_cast<int>(lg.w.size());
  std::vector<double> draws(n_draws * 5);
  std::vector<double> lw(n_draws);
  std::vector<double> prob(static_cast<std::size_t>(n_draws) * m);
  const double df = lg.df;

  for (int round = 0; round <= lg.rounds; ++round) {
    const bool last = round == lg.rounds;
    double lchol[25];
    if (!cb12_cholesky(scale, lchol)) {
      // Fall back to the curvature at the mode.
      double lc[25], ic[25];
      cb12_cholesky(curv, lc);
      cb12_chol_inverse(lc, ic);
      for (int i = 0; i < 25; ++i) scale[i] = ic[i];
      cb12_cholesky(scale, lchol);
    }
    double log_det = 0.0;
    for (int k = 0; k < 5; ++k) log_det += std::log(lchol[k * 5 + k]);
    double lw_max = -std::numeric_limits<double>::infinity();
    for (int d = 0; d < n_draws; ++d) {
      const double* zz = &lg.z[d * 5];
      const double r = std::sqrt(lg.w[d] / df);
      double* t = &draws[d * 5];
      double maha = 0.0;
      for (int a = 0; a < 5; ++a) {
        double sum = 0.0;
        for (int b = 0; b <= a; ++b) sum += lchol[a * 5 + b] * zz[b];
        t[a] = mu[a] + sum / r;
        maha += (zz[a] / r) * (zz[a] / r);
      }
      // Log posterior at t, and the response probabilities for the last
      // proposal: with e = exp(-|eta|), log(1 + exp(eta)) = max(eta, 0) +
      // log1p(e) and p = 1 / (1 + e) or e / (1 + e).
      double lp = 0.0;
      for (int i = 0; i < m; ++i) {
        double eta = 0.0;
        for (int k = 0; k < 5; ++k) eta += f[i * 5 + k] * t[k];
        const double e = std::exp(-std::fabs(eta));
        const double l1 = (eta > 0.0 ? eta : 0.0) + std::log1p(e);
        lp += xv[i] * eta - nv[i] * l1;
        if (last) {
          prob[static_cast<std::size_t>(d) * m + i] =
            eta >= 0.0 ? 1.0 / (1.0 + e) : e / (1.0 + e);
        }
      }
      for (int k = 0; k < 5; ++k) lp -= std::log1p(t[k] * t[k] / s2[k]);
      const double prop_lp = -0.5 * (df + 5.0) * std::log1p(maha / df) -
        log_det;
      lw[d] = lp - prop_lp;
      if (lw[d] > lw_max) lw_max = lw[d];
    }
    double wsum = 0.0;
    for (int d = 0; d < n_draws; ++d) {
      lw[d] = std::exp(lw[d] - lw_max);
      wsum += lw[d];
    }
    for (int d = 0; d < n_draws; ++d) lw[d] /= wsum;

    if (last) break;
    // Refit the proposal to the weighted points.
    for (int k = 0; k < 5; ++k) {
      double sum = 0.0;
      for (int d = 0; d < n_draws; ++d) sum += lw[d] * draws[d * 5 + k];
      mu[k] = sum;
    }
    for (int a = 0; a < 5; ++a) {
      for (int b = 0; b <= a; ++b) {
        double sum = 0.0;
        for (int d = 0; d < n_draws; ++d) {
          sum += lw[d] * (draws[d * 5 + a] - mu[a]) *
            (draws[d * 5 + b] - mu[b]);
        }
        const double v = sum * (df > 2.0 ? (df - 2.0) / df : 1.0);
        scale[a * 5 + b] = v;
        scale[b * 5 + a] = v;
      }
      scale[a * 5 + a] += 1e-9;
    }
  }

  // Posterior means of the tried combinations; the others are not needed
  // and are set to NaN.
  int i = 0;
  for (int c = 0; c < n_cells; ++c) {
    if (n[c] == 0) {
      pe[c] = std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    double sum = 0.0;
    for (int d = 0; d < n_draws; ++d) {
      sum += lw[d] * prob[static_cast<std::size_t>(d) * m + i];
    }
    pe[c] = sum;
    ++i;
  }
}

}  // namespace simfastboin

#endif  // SIMFASTBOIN_COMB_BOIN12_CORE_H
