#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

//' Batch Stationary Point Solver
//'
//' @param h_array Hessian array with dimensions draws x p x p.
//' @param b_matrix Linear coefficient matrix with dimensions draws x p.
//' @param kappa_thresh Condition-number threshold for near-singularity.
//'
//' @return Numeric matrix of stationary points (draws x p), with NA rows for
//' near-singular or unsolved systems.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix stationary_points_batch(NumericVector h_array,
                                      NumericMatrix b_matrix,
                                      double kappa_thresh) {
  IntegerVector dims = h_array.attr("dim");
  if (dims.size() != 3) {
    stop("h_array must be a 3D array with dims (draws, p, p).");
  }

  int n_draws = dims[0];
  int p = dims[1];
  if (dims[2] != p) {
    stop("h_array must have square Hessian slices (p x p).");
  }
  if (b_matrix.nrow() != n_draws || b_matrix.ncol() != p) {
    stop("b_matrix must have dimensions (draws x p) matching h_array.");
  }

  NumericMatrix out(n_draws, p);
  std::fill(out.begin(), out.end(), NA_REAL);

  int n_rhs = 1;
  int lda = p;
  int ldb = p;
  std::vector<int> ipiv(p);
  std::vector<double> a(p * p);
  std::vector<double> b(p);

  for (int d = 0; d < n_draws; ++d) {
    for (int i = 0; i < p; ++i) {
      b[i] = -0.5 * b_matrix(d, i);
      for (int j = 0; j < p; ++j) {
        a[i + p * j] = h_array[d + n_draws * (i + p * j)];
      }
    }

    int info = 0;
    F77_CALL(dgesv)(&p, &n_rhs, a.data(), &lda, ipiv.data(), b.data(), &ldb, &info);

    if (info != 0) {
      continue;
    }

    // Use an LU-diagonal ratio as a low-cost near-singularity proxy.
    double min_diag = INFINITY;
    double max_diag = 0.0;
    for (int i = 0; i < p; ++i) {
      double di = std::fabs(a[i + p * i]);
      if (!std::isfinite(di) || di <= 0.0) {
        min_diag = 0.0;
        break;
      }
      min_diag = std::min(min_diag, di);
      max_diag = std::max(max_diag, di);
    }

    if (!std::isfinite(min_diag) || !std::isfinite(max_diag) || min_diag <= 0.0) {
      continue;
    }

    double kappa = max_diag / min_diag;
    if (!std::isfinite(kappa) || kappa > kappa_thresh) {
      continue;
    }

    for (int i = 0; i < p; ++i) {
      out(d, i) = std::isfinite(b[i]) ? b[i] : NA_REAL;
    }
  }

  return out;
}