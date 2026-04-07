# Test setup and utilities

# Ensure testthat is loaded
library(testthat)

# Load the package
library(brsm)

# Helper function: check if object is a data frame
expect_data_frame <- function(object, nrows = NULL, ncols = NULL) {
  expect_s3_class(object, "data.frame")
  if (!is.null(nrows)) {
    expect_equal(nrow(object), nrows)
  }
  if (!is.null(ncols)) {
    expect_equal(ncol(object), ncols)
  }
}

# Helper function: suppress messages/warnings during tests
suppress_output <- function(expr) {
  suppressMessages(suppressWarnings(expr))
}

# Helper: Generate simple RSA simulation data
generate_simulation_data <- function(n = 50, seed = 123) {
  set.seed(seed)

  x1 <- runif(n, -2, 2)
  x2 <- runif(n, -2, 2)

  # True RSA model: y ~ 50 + 0.5*x1 + 0.5*x2 - 0.3*x1^2 - 0.25*x2^2 + noise
  y <- 50 + 0.5 * x1 + 0.5 * x2 - 0.3 * x1^2 - 0.25 * x2^2 + rnorm(n, 0, 5)

  data.frame(x1 = x1, x2 = x2, y = y)
}

# Helper: Generate mock brmsfit object (for testing without actually fitting)
# Note: This is a placeholder; real testing would use actual fitted models
create_mock_brmsfit <- function() {
  skip(
    "Mock brmsfit creation not fully implemented; ",
    "use actual models in integration tests"
  )
}

# Helper: gate expensive brms fitting tests behind dependency + env flag
brms_tests_enabled <- function() {
  flag <- tolower(Sys.getenv("BRSM_RUN_BRMS_TESTS", "false"))
  flag %in% c("1", "true", "yes", "y")
}

skip_if_no_brms_tests <- function() {
  skip_if_not_installed("brms")
  if (!brms_tests_enabled()) {
    skip("Set BRSM_RUN_BRMS_TESTS=true to run brms model-fitting tests.")
  }
}
