#!/usr/bin/env Rscript
#
# Test pipe-friendly design with S3 methods for brsm_fit objects
#

library(brsm)

# Generate test data
set.seed(42)
dat <- data.frame(
  x1 = runif(30, -2, 2),
  x2 = runif(30, -2, 2)
)
dat$y <- 50 + 0.5 * dat$x1 + 0.5 * dat$x2 - 0.3 * dat$x1^2 - 0.25 * dat$x2^2 + rnorm(30, 0, 5)

cat("=== Test 1: congruence_parameters without explicit factor_names ===\n")
# Create mock draws for quick testing (avoid Bayesian fitting)
draws_raw <- data.frame(
  b_Intercept = rnorm(100, 50, 5),
  b_x1 = rnorm(100, 0.5, 0.1),
  b_x2 = rnorm(100, 0.5, 0.1),
  "b_I(x1^2)" = rnorm(100, -0.3, 0.05),
  "b_I(x2^2)" = rnorm(100, -0.25, 0.05),
  "b_x1:x2" = rnorm(100, 0.05, 0.02),
  check.names = FALSE
)

# Create a minimal brsm_fit object for testing
fit_obj <- list(
  fit = NULL,
  formula = y ~ x1 + x2 + I(x1^2) + I(x2^2) + x1:x2,
  response = "y",
  factor_names = c("x1", "x2"),
  ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
  call = NULL,
  draws_df = draws_raw
)
class(fit_obj) <- "brsm_fit"

# Store draws so as_brsm_draws can access them
fit_obj$fit <- list(draws_raw)
class(fit_obj$fit) <- "brmsfit"

cat("Testing S3 dispatch for congruence_parameters:\n")
tryCatch(
  {
    # This should work without explicit factor_names because S3 method extracts from fit_obj
    result <- congruence_parameters(fit_obj)
    cat("✓ congruence_parameters(fit_obj) works!\n")
    cat("  Columns:", names(result), "\n")
    cat("  Rows:", nrow(result), "\n")
  },
  error = function(e) {
    cat("✗ Error:", e$message, "\n")
  }
)

cat("\n=== Test 2: Chain multiple analyses with pipe operator ===\n")
cat("Testing pipe composition: fit_obj |> stationary_point()\n")
tryCatch(
  {
    sp <- stationary_point(fit_obj)
    cat("✓ stationary_point(fit_obj) works!\n")
    cat("  Columns:", names(sp), "\n")
    cat("  Rows:", nrow(sp), "\n")
  },
  error = function(e) {
    cat("✗ Error:", e$message, "\n")
  }
)

cat("\nTesting pipe composition: fit_obj |> classify_stationarity_point()\n")
tryCatch(
  {
    csp <- classify_stationarity_point(fit_obj)
    cat("✓ classify_stationarity_point(fit_obj) works!\n")
    cat("  Unique classes:", unique(csp$classification), "\n")
  },
  error = function(e) {
    cat("✗ Error:", e$message, "\n")
  }
)

cat("\nTesting pipe composition: fit_obj |> posterior_ridge_analysis()\n")
tryCatch(
  {
    ridge <- posterior_ridge_analysis(fit_obj)
    cat("✓ posterior_ridge_analysis(fit_obj) works!\n")
    cat("  Dimensions:", dim(ridge), "\n")
  },
  error = function(e) {
    cat("✗ Error:", e$message, "\n")
  }
)

cat("\nTesting pipe composition: fit_obj |> credible_optimum_region()\n")
tryCatch(
  {
    cr <- credible_optimum_region(fit_obj)
    cat("✓ credible_optimum_region(fit_obj) works!\n")
    cat("  Columns:", names(cr), "\n")
  },
  error = function(e) {
    cat("✗ Error:", e$message, "\n")
  }
)

cat("\nTesting pipe composition: fit_obj |> steepest_ascent()\n")
# Note: steepest_ascent has default start = c(x1=0, x2=0) when not provided
tryCatch(
  {
    sa <- steepest_ascent(fit_obj)
    cat("✓ steepest_ascent(fit_obj) works!\n")
    cat("  Dimensions:", dim(sa), "\n")
  },
  error = function(e) {
    cat("✗ Error:", e$message, "\n")
  }
)

cat("\n=== Test 3: Backward compatibility - bare data frame + factor_names ===\n")
draws_df <- as_brsm_draws(draws_raw, factor_names = c("x1", "x2"))
cat("Testing backward compat: congruence_parameters(draws_df, factor_names = c('x1', 'x2'))\n")
tryCatch(
  {
    result_compat <- congruence_parameters(draws_df, factor_names = c("x1", "x2"))
    cat("✓ Backward compat works!\n")
    cat("  Columns:", names(result_compat), "\n")
  },
  error = function(e) {
    cat("✗ Error:", e$message, "\n")
  }
)

cat("\n✓ All pipe-friendly design tests completed!\n")
