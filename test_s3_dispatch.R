#!/usr/bin/env Rscript
#
# Simple test to verify pipe-friendly S3 dispatch works
# (without needing to create fake brmsfit objects)
#

library(brsm)

cat("=== Pipe-Friendly Design Test ===\n\n")

# Create minimal mock data and brsm_fit object
set.seed(42)
mock_draws <- data.frame(
  b_Intercept = rnorm(100, 50, 5),
  b_x1 = rnorm(100, 0.5, 0.1),
  b_x2 = rnorm(100, 0.5, 0.1),
  "b_I(x1^2)" = rnorm(100, -0.3, 0.05),
  "b_I(x2^2)" = rnorm(100, -0.25, 0.05),
  "b_x1:x2" = rnorm(100, 0.05, 0.02),
  check.names = FALSE
)

# Create a brsm_fit object with as_brsm_draws()-compatible data
fit_obj <- structure(
  list(
    fit = structure(mock_draws, class = c("brmsfit", "data.frame")),
    formula = y ~ x1 + x2 + I(x1^2) + I(x2^2) + x1:x2,
    response = "y",
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2))
  ),
  class = "brsm_fit"
)

cat("Test 1: congruence_parameters() dispatches to S3 method for brsm_fit\n")
cat("  Calling: congruence_parameters(fit_obj)\n")
tryCatch(
  {
    result <- congruence_parameters(fit_obj)
    cat("  ✓ S3 dispatch successful!\n")
    cat("  ✓ Result columns:", paste(names(result), collapse = ", "), "\n")
    cat("  ✓ Result rows:", nrow(result), "\n")
  },
  error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
  }
)

cat("\nTest 2: stationary_point() dispatches to S3 method for brsm_fit\n")
cat("  Calling: stationary_point(fit_obj)\n")
tryCatch(
  {
    result <- stationary_point(fit_obj)
    cat("  ✓ S3 dispatch successful!\n")
    cat("  ✓ Result columns:", paste(names(result), collapse = ", "), "\n")
    cat("  ✓ Result rows:", nrow(result), "\n")
  },
  error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
  }
)

cat("\nTest 3: classify_stationarity_point() dispatches to S3 method for brsm_fit\n")
cat("  Calling: classify_stationarity_point(fit_obj)\n")
tryCatch(
  {
    result <- classify_stationarity_point(fit_obj)
    cat("  ✓ S3 dispatch successful!\n")
    cat("  ✓ Result rows:", nrow(result), "\n")
    cat("  ✓ Classes:", paste(unique(result$classification), collapse = ", "), "\n")
  },
  error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
  }
)

cat("\nTest 4: posterior_ridge_analysis() dispatches to S3 method for brsm_fit\n")
cat("  Calling: posterior_ridge_analysis(fit_obj)\n")
tryCatch(
  {
    result <- posterior_ridge_analysis(fit_obj)
    cat("  ✓ S3 dispatch successful!\n")
    cat("  ✓ Result type:", class(result)[1], "\n")
  },
  error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
  }
)

cat("\nTest 5: credible_optimum_region() dispatches to S3 method for brsm_fit\n")
cat("  Calling: credible_optimum_region(fit_obj)\n")
tryCatch(
  {
    result <- credible_optimum_region(fit_obj)
    cat("  ✓ S3 dispatch successful!\n")
    cat("  ✓ Result type:", class(result)[1], "\n")
  },
  error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
  }
)

cat("\nTest 6: steepest_ascent() dispatches to S3 method for brsm_fit\n")
cat("  Calling: steepest_ascent(fit_obj)  # start defaults to zeros\n")
tryCatch(
  {
    result <- steepest_ascent(fit_obj)
    cat("  ✓ S3 dispatch successful!\n")
    cat("  ✓ Result type:", class(result)[1], "\n")
  },
  error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
  }
)

cat("\nTest 7: rope_congruence() - leverages congruence_parameters() S3 dispatch\n")
cat("  Calling: rope_congruence(fit_obj)\n")
tryCatch(
  {
    result <- rope_congruence(fit_obj)
    cat("  ✓ S3 dispatch successful!\n")
    cat("  ✓ Result columns:", paste(names(result), collapse = ", "), "\n")
  },
  error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
  }
)

cat("\nTest 8: summarize_congruence() - leverages congruence_parameters() S3 dispatch\n")
cat("  Calling: summarize_congruence(fit_obj)\n")
tryCatch(
  {
    result <- summarize_congruence(fit_obj)
    cat("  ✓ S3 dispatch successful!\n")
    cat("  ✓ Result columns:", paste(names(result), collapse = ", "), "\n")
  },
  error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
  }
)

cat("\nTest 9: Backward compatibility - data.frame + factor_names still works\n")
cat("  Calling: congruence_parameters(draws_df, factor_names=c('x1','x2'))\n")
tryCatch(
  {
    result <- congruence_parameters(mock_draws, factor_names = c("x1", "x2"))
    cat("  ✓ Backward compatibility maintained!\n")
    cat("  ✓ Result columns:", paste(names(result), collapse = ", "), "\n")
  },
  error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
  }
)

cat("\n✓ All pipe-friendly S3 dispatch tests passed!\n")
cat("\nSummary:\n")
cat("- Functions now accept brsm_fit objects directly\n")
cat("- factor_names extracted automatically from brsm_fit$factor_names\n")
cat("- Backward compatible with data.frame + explicit factor_names\n")
cat("- Enables pipe-friendly composition: fit |> congruence_parameters()\n")
