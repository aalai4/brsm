#include <Rcpp.h>
using namespace Rcpp;

//' Convert Prediction Matrix to Long-Format Data Frame
//'
//' @param pred_matrix A matrix of predictions (n_draws x n_points)
//' @param draw_ids Integer vector of draw identifiers
//' @param point_ids Integer vector of point identifiers
//'
//' @return A data frame with columns: draw, point_id, estimate
//' @keywords internal
// [[Rcpp::export]]
DataFrame matrix_to_long_format(NumericMatrix pred_matrix,
                                 IntegerVector draw_ids,
                                 IntegerVector point_ids) {
  int n_draws = pred_matrix.nrow();
  int n_points = pred_matrix.ncol();
  int n_total = n_draws * n_points;

  // Pre-allocate output vectors
  IntegerVector out_draw(n_total);
  IntegerVector out_point(n_total);
  NumericVector out_estimate(n_total);

  // Fill vectors efficiently by iterating row-by-row through the matrix
  // This maintains cache locality better than column-wise iteration
  int idx = 0;
  for (int i = 0; i < n_draws; ++i) {
    int draw_val = draw_ids[i];
    for (int j = 0; j < n_points; ++j) {
      out_draw[idx] = draw_val;
      out_point[idx] = point_ids[j];
      out_estimate[idx] = pred_matrix(i, j);
      ++idx;
    }
  }

  // Create and return the data frame
  return DataFrame::create(
    Named("draw") = out_draw,
    Named("point_id") = out_point,
    Named("estimate") = out_estimate
  );
}
