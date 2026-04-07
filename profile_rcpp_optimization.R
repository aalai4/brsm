#!/usr/bin/env Rscript
#
# Profiling script to measure Rcpp optimization impact on predict_surface()
#

library(brsm)
library(brms)

set.seed(42)

# Generate realistic test data
n_obs <- 50
data <- data.frame(
  x1 = rnorm(n_obs, 0, 1),
  x2 = rnorm(n_obs, 0, 1),
  y = rnorm(n_obs, 0, 1)
)

# Fit a basic model (using weak priors for speed)
cat("Fitting test model...\n")
fit <- brm(
  y ~ x1 + x2 + I(x1^2) + I(x2^2) + x1:x2,
  data = data,
  chains = 2,
  iter = 500,
  warmup = 250,
  cores = 2,
  refresh = 0,
  backend = "rstan"
)

cat("Extracting draws...\n")
draws <- as_brsm_draws(fit)

cat(sprintf("Draws shape: %d × %d\n", nrow(draws), ncol(draws)))

# Create prediction grid
grid <- expand.grid(
  x1 = seq(-2, 2, length.out = 64),
  x2 = seq(-2, 2, length.out = 64)
)

cat(sprintf("Grid size: %d points\n", nrow(grid)))

# Benchmark: Full output with Rcpp vectorization
cat("\n=== BENCHMARK: Full output with Rcpp optimization ===\n")
time_rcpp <- system.time({
  for (rep in 1:3) {
    result <- predict_surface(draws, c("x1", "x2"), grid, summary = FALSE)
  }
})
cat(sprintf("Time for 3 reps: %.2f seconds\n", time_rcpp["elapsed"]))
cat(sprintf("Time per call: %.2f seconds\n", time_rcpp["elapsed"] / 3))
cat(sprintf("Output dimensions: %d rows × %d cols\n", nrow(result), ncol(result)))

# Benchmark: With max_draws (reduced cardinality)
cat("\n=== BENCHMARK: With draw_subset (first 500 draws) ===\n")
time_subset <- system.time({
  for (rep in 1:3) {
    result_subset <- predict_surface(
      draws, c("x1", "x2"), grid,
      summary = FALSE,
      draw_subset = 1:500
    )
  }
})
cat(sprintf("Time for 3 reps: %.2f seconds\n", time_subset["elapsed"]))
cat(sprintf("Time per call: %.2f seconds\n", time_subset["elapsed"] / 3))
cat(sprintf("Output dimensions: %d rows × %d cols\n", nrow(result_subset), ncol(result_subset)))

# Benchmark: With return_matrix (skip long-format assembly)
cat("\n=== BENCHMARK: With return_matrix=TRUE ===\n")
time_matrix <- system.time({
  for (rep in 1:3) {
    result_matrix <- predict_surface(
      draws, c("x1", "x2"), grid,
      summary = FALSE,
      return_matrix = TRUE
    )
  }
})
cat(sprintf("Time for 3 reps: %.2f seconds\n", time_matrix["elapsed"]))
cat(sprintf("Time per call: %.2f seconds\n", time_matrix["elapsed"] / 3))
cat(sprintf("Output dimensions: %d × %d (matrix)\n", nrow(result_matrix), ncol(result_matrix)))

# Summary statistics
cat("\n=== SUMMARY ===\n")
cat(sprintf("Full output: %.3f s per call\n", time_rcpp["elapsed"] / 3))
cat(sprintf(
  "With draw_subset: %.3f s per call (%.1f×)\n",
  time_subset["elapsed"] / 3,
  (time_rcpp["elapsed"] / 3) / (time_subset["elapsed"] / 3)
))
cat(sprintf(
  "With return_matrix: %.3f s per call (%.1f×)\n",
  time_matrix["elapsed"] / 3,
  (time_rcpp["elapsed"] / 3) / (time_matrix["elapsed"] / 3)
))
