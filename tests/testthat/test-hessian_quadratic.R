# Tests for hessian_quadratic (no brms needed — uses mock draws)

.make_hess_draws <- function(n = 30, seed = 42) {
  set.seed(seed)
  df <- data.frame(
    b_Intercept = rnorm(n, 5, 0.2),
    b_x1        = rnorm(n, 2, 0.1),
    b_x2        = rnorm(n, -1, 0.1),
    check.names = FALSE
  )
  df[["b_I(x1^2)"]] <- rnorm(n, -1.5, 0.05)
  df[["b_I(x2^2)"]] <- rnorm(n, -2.0, 0.05)
  df[["b_x1:x2"]]   <- rnorm(n,  0.3, 0.05)
  df
}

test_that("hessian_quadratic returns correct column names", {
  draws  <- .make_hess_draws(n = 10)
  result <- hessian_quadratic(draws, factor_names = c("x1", "x2"))
  expect_s3_class(result, "data.frame")
  expect_named(result, c("draw", "row_factor", "col_factor", "value"))
})

test_that("hessian_quadratic returns p^2 rows per draw", {
  n      <- 20
  draws  <- .make_hess_draws(n = n)
  result <- hessian_quadratic(draws, factor_names = c("x1", "x2"))
  # 2 factors -> 2x2 Hessian -> 4 entries per draw
  expect_equal(nrow(result), n * 4L)
  expect_equal(sort(unique(result$draw)), seq_len(n))
})

test_that("hessian_quadratic diagonal = 2 * quadratic coefficients", {
  draws <- data.frame(b_Intercept = 5, b_x1 = 2, b_x2 = -1,
                      check.names = FALSE)
  draws[["b_I(x1^2)"]] <- -1.5
  draws[["b_I(x2^2)"]] <- -2.0
  draws[["b_x1:x2"]]   <-  0.4

  result <- hessian_quadratic(draws, factor_names = c("x1", "x2"))

  h11 <- result$value[result$row_factor == "x1" & result$col_factor == "x1"]
  h22 <- result$value[result$row_factor == "x2" & result$col_factor == "x2"]
  expect_equal(h11, 2 * (-1.5), tolerance = 1e-10)
  expect_equal(h22, 2 * (-2.0), tolerance = 1e-10)
})

test_that("hessian_quadratic off-diagonals are symmetric and equal b_x1:x2", {
  draws <- data.frame(b_Intercept = 5, b_x1 = 2, b_x2 = -1,
                      check.names = FALSE)
  draws[["b_I(x1^2)"]] <- -1.5
  draws[["b_I(x2^2)"]] <- -2.0
  draws[["b_x1:x2"]]   <-  0.4

  result <- hessian_quadratic(draws, factor_names = c("x1", "x2"))

  h12 <- result$value[result$row_factor == "x1" & result$col_factor == "x2"]
  h21 <- result$value[result$row_factor == "x2" & result$col_factor == "x1"]
  expect_equal(h12, h21, tolerance = 1e-10)
  expect_equal(h12, 0.4, tolerance = 1e-10)
})

test_that("hessian_quadratic works for a single factor", {
  draws <- data.frame(b_Intercept = 5, b_x1 = 2, check.names = FALSE)
  draws[["b_I(x1^2)"]] <- -1.5

  result <- hessian_quadratic(draws, factor_names = "x1")
  expect_equal(nrow(result), 1L)
  expect_equal(result$value, 2 * (-1.5), tolerance = 1e-10)
})

test_that("hessian_quadratic scales correctly across multiple draws", {
  n  <- 5
  b2 <- seq(-1, -2, length.out = n)   # known quadratic coefs for x1
  draws <- data.frame(b_Intercept = rep(5, n), b_x1 = rep(1, n),
                      check.names = FALSE)
  draws[["b_I(x1^2)"]] <- b2

  result <- hessian_quadratic(draws, factor_names = "x1")
  expect_equal(result$value, 2 * b2, tolerance = 1e-10)
})

test_that("hessian_quadratic errors on invalid factor_names", {
  draws <- .make_hess_draws(n = 5)
  expect_error(hessian_quadratic(draws, factor_names = character(0)))
  expect_error(hessian_quadratic(draws, factor_names = NULL))
})
