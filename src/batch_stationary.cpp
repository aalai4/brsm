#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

enum StationaryStatus {
  STATUS_OK = 0,
  STATUS_LAPACK_FAIL = 1,
  STATUS_INVALID_LU_DIAG = 2,
  STATUS_KAPPA_EXCEEDED = 3
};

List stationary_points_batch_impl(NumericVector h_array,
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

  NumericMatrix x_star(n_draws, p);
  std::fill(x_star.begin(), x_star.end(), NA_REAL);

  IntegerVector status_code(n_draws, STATUS_LAPACK_FAIL);
  NumericVector kappa_proxy(n_draws, NA_REAL);

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
      status_code[d] = STATUS_LAPACK_FAIL;
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
      status_code[d] = STATUS_INVALID_LU_DIAG;
      continue;
    }

    double kappa = max_diag / min_diag;
    kappa_proxy[d] = kappa;
    if (!std::isfinite(kappa) || kappa > kappa_thresh) {
      status_code[d] = STATUS_KAPPA_EXCEEDED;
      continue;
    }

    for (int i = 0; i < p; ++i) {
      x_star(d, i) = std::isfinite(b[i]) ? b[i] : NA_REAL;
    }
    status_code[d] = STATUS_OK;
  }

  return List::create(
    _["x_star"] = x_star,
    _["status_code"] = status_code,
    _["kappa_proxy"] = kappa_proxy
  );
}

}  // namespace

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
  return as<NumericMatrix>(stationary_points_batch_impl(h_array, b_matrix, kappa_thresh)["x_star"]);
}

//' Detailed Batch Stationary Point Solver
//'
//' @param h_array Hessian array with dimensions draws x p x p.
//' @param b_matrix Linear coefficient matrix with dimensions draws x p.
//' @param kappa_thresh Condition-number threshold for near-singularity.
//'
//' @return A list with stationary points, per-draw status codes, and kappa proxies.
//' @keywords internal
// [[Rcpp::export]]
List stationary_points_batch_details(NumericVector h_array,
                                     NumericMatrix b_matrix,
                                     double kappa_thresh) {
  return stationary_points_batch_impl(h_array, b_matrix, kappa_thresh);
}