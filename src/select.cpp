// Rcpp glue for the selection of the MTD from estimates that need not be
// whole numbers, see select_core.h.

#include <Rcpp.h>
#include "select_core.h"

using namespace Rcpp;

// [[Rcpp::export]]
List select_mtd_estimate_cpp(LogicalMatrix eliminated,
                             IntegerMatrix n_pts,
                             NumericMatrix est_y,
                             NumericMatrix est_n,
                             double target,
                             int min_mtd_sample) {

  const int n_trials = n_pts.nrow();
  const int n_doses = n_pts.ncol();

  IntegerVector mtd(n_trials);
  IntegerVector reason(n_trials);

  std::vector<char> el(n_doses);
  std::vector<int> nn(n_doses);
  std::vector<double> ey(n_doses), en(n_doses);

  for (int t = 0; t < n_trials; ++t) {
    for (int j = 0; j < n_doses; ++j) {
      el[j] = (eliminated(t, j) == 1) ? 1 : 0;
      nn[j] = n_pts(t, j);
      ey[j] = est_y(t, j);
      en[j] = est_n(t, j);
    }
    int code = 0;
    const int d = simfastboin::select_mtd_estimate(&el[0], &nn[0], &ey[0],
                                                   &en[0], n_doses, target,
                                                   min_mtd_sample, code);
    mtd[t] = (d == 0) ? NA_INTEGER : d;
    reason[t] = code;
  }

  return List::create(_["mtd"] = mtd, _["reason"] = reason);
}
