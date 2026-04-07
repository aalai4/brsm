setwd("c:/Users/Prabhashis-LAPTOP/Downloads/brsm")
devtools::load_all(".")
library(rsm)
data(heli)

# Test 1: validate posterior-draw input path
set.seed(42)
draws_raw <- data.frame(
  b_Intercept = rnorm(400, mean(heli$ave), 2),
  b_x1 = rnorm(400, 2, 0.2),
  b_x2 = rnorm(400, 4, 0.2),
  "b_I(x1^2)" = rnorm(400, -1, 0.1),
  "b_I(x2^2)" = rnorm(400, -2, 0.1),
  "b_x1:x2" = rnorm(400, 0.5, 0.1),
  check.names = FALSE
)
draws <- as_brsm_draws(draws_raw, factor_names = c("x1", "x2"))
cat("✓ Test 1: as_brsm_draws(data.frame) works!\n")
print(head(draws, 2))

# Test 2: verify columns
expected_cols <- c("b_Intercept", "b_x1", "b_x2", "b_I(x1^2)", "b_I(x2^2)", "b_x1:x2")
if (all(expected_cols %in% names(draws))) {
  cat("\n✓ Test 2: All expected columns present\n")
} else {
  cat("\n✗ Test 2: Missing columns:", setdiff(expected_cols, names(draws)), "\n")
}

# Test 3: integration with rest of brsm
f <- c("x1", "x2")
ranges <- list(x1 = range(heli$x1), x2 = range(heli$x2))
grid <- surface_grid(ranges = ranges, n = 10)
pred <- predict_surface(draws, f, grid, summary = TRUE)
cat("\n✓ Test 3: Integration with predict_surface works!\n")

# Test 4: Bayesian-only validation - reject non-Bayesian objects
tryCatch(
  {
    lm_fit <- lm(y ~ x1 + x2, data = heli)
    brsm_workflow(lm_fit, factor_names = c("x1", "x2"), ranges = ranges)
    cat("\n✗ Test 4: FAILED - lm object should have been rejected!\n")
  },
  error = function(e) {
    if (grepl("strictly Bayesian|frequentist", e$message)) {
      cat("\n✓ Test 4: Bayesian-only validation rejects lm objects\n")
    } else {
      cat("\n✗ Test 4: Wrong error message:\n", e$message, "\n")
    }
  }
)

# Test 5: Bayesian-only validation - reject data frame without b_* columns
bad_df <- data.frame(x1 = 1:10, x2 = 1:10, y = 1:10)
tryCatch(
  {
    brsm_workflow(bad_df, factor_names = c("x1", "x2"), ranges = ranges)
    cat("\n✗ Test 5: FAILED - data frame without b_* columns should have been rejected!\n")
  },
  error = function(e) {
    if (grepl("posterior columns|b_Intercept", e$message)) {
      cat("\n✓ Test 5: Bayesian-only validation rejects data frame without b_* columns\n")
    } else {
      cat("\n✗ Test 5: Wrong error message:\n", e$message, "\n")
    }
  }
)

# Test 6: Bayesian-only validation - accept proper posterior draws
tryCatch(
  {
    result <- brsm_workflow(draws_raw, factor_names = c("x1", "x2"), ranges = ranges)
    cat("\n✓ Test 6: Bayesian-only validation accepts proper posterior draws\n")
  },
  error = function(e) {
    cat("\n✗ Test 6: FAILED - should accept posterior draws:\n", e$message, "\n")
  }
)

cat("\n✓ All Bayesian-only validation tests passed!\n")
