# Tests for gradient_quadratic (no brms needed - uses mock draws)

.make_grad_draws <- function(n = 10, seed = 123) {
  set.seed(seed)
  df <- data.frame(
    b_Intercept = rnorm(n, 0, 0.1),
    b_x1 = rnorm(n, 1, 0.1),
    b_x2 = rnorm(n, -1, 0.1),
    check.names = FALSE
  )
  df[["b_I(x1^2)"]] <- rnorm(n, -0.5, 0.05)
  df[["b_I(x2^2)"]] <- rnorm(n, -0.25, 0.05)
  df[["b_x1:x2"]] <- rnorm(n, 0.3, 0.05)
  df
}

test_that("gradient_quadratic errors when x is not numeric", {
  draws <- .make_grad_draws(5)
  expect_error(
    gradient_quadratic(draws, x = "bad", factor_names = c("x1", "x2")),
    "x must be numeric"
  )
})

test_that("gradient_quadratic errors when x contains NA", {
  draws <- .make_grad_draws(5)
  expect_error(
    gradient_quadratic(draws, x = c(NA_real_, 0), factor_names = c("x1", "x2")),
    "x cannot contain NA values"
  )
})

test_that("gradient_quadratic errors when linear columns are non-numeric", {
  draws <- .make_grad_draws(3)
  draws$b_x1 <- as.character(draws$b_x1)

  expect_error(
    gradient_quadratic(draws, x = c(0, 0), factor_names = c("x1", "x2")),
    "draws must contain numeric columns for linear and quadratic terms"
  )
})

test_that("gradient_quadratic errors on x dimension mismatch", {
  draws <- .make_grad_draws(5)
  expect_error(
    gradient_quadratic(draws, x = matrix(c(0, 1, 2), nrow = 1), factor_names = c("x1", "x2")),
    "number of columns of x must match number of factors"
  )
})

test_that("gradient_quadratic normalize handles zero-norm rows", {
  draws <- data.frame(
    b_Intercept = 0,
    b_x1 = 0,
    b_x2 = 0,
    check.names = FALSE
  )
  draws[["b_I(x1^2)"]] <- 0
  draws[["b_I(x2^2)"]] <- 0
  draws[["b_x1:x2"]] <- 0

  out <- gradient_quadratic(
    draws,
    x = c(0, 0),
    factor_names = c("x1", "x2"),
    normalize = TRUE
  )

  expect_equal(out$d_dx1, 0)
  expect_equal(out$d_dx2, 0)
  expect_false(any(is.nan(out$d_dx1)))
  expect_false(any(is.nan(out$d_dx2)))
})
