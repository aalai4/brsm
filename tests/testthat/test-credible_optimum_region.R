# Tests for credible_optimum_region (no brms needed — uses mock draws)

.make_cor_draws <- function(n = 50, seed = 77) {
  set.seed(seed)
  df <- data.frame(
    b_Intercept = rnorm(n, 5, 0.1),
    b_x1        = rnorm(n, 1, 0.05),
    b_x2        = rnorm(n, 0.5, 0.05),
    check.names = FALSE
  )
  # Strongly negative-definite Hessian -> guaranteed maximum -> stable stationary point
  df[["b_I(x1^2)"]] <- rnorm(n, -3, 0.1)
  df[["b_I(x2^2)"]] <- rnorm(n, -3, 0.1)
  df[["b_x1:x2"]]   <- rnorm(n,  0.2, 0.05)
  df
}

test_that("credible_optimum_region.default returns a data frame with one row per factor", {
  draws  <- .make_cor_draws()
  result <- suppressWarnings(
    credible_optimum_region(draws, factor_names = c("x1", "x2"))
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2L)
  expect_equal(rownames(result), c("x1", "x2"))
})

test_that("credible_optimum_region summary result has mean, sd, and quantile columns", {
  draws  <- .make_cor_draws()
  result <- suppressWarnings(
    credible_optimum_region(draws, factor_names = c("x1", "x2"),
                            probs = c(0.025, 0.975))
  )
  expect_true("mean"    %in% names(result))
  expect_true("sd"      %in% names(result))
  expect_true("q2.5"    %in% names(result))
  expect_true("q97.5"   %in% names(result))
})

test_that("credible_optimum_region summary=FALSE returns raw stationary points", {
  draws  <- .make_cor_draws()
  result <- suppressWarnings(
    credible_optimum_region(draws, factor_names = c("x1", "x2"),
                            summary = FALSE)
  )
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2L)   # one column per factor
})

test_that("credible_optimum_region errors when factor_names is NULL for data frame", {
  draws <- .make_cor_draws()
  expect_error(
    credible_optimum_region(draws),
    "factor_names"
  )
})

test_that("credible_optimum_region custom probs produce correct column names", {
  draws  <- .make_cor_draws()
  result <- suppressWarnings(
    credible_optimum_region(draws, factor_names = c("x1", "x2"),
                            probs = c(0.1, 0.9))
  )
  expect_true("q10.0" %in% names(result))
  expect_true("q90.0" %in% names(result))
})

test_that("credible_optimum_region 90% interval is narrower than 95%", {
  draws   <- .make_cor_draws()
  res_90  <- suppressWarnings(
    credible_optimum_region(draws, factor_names = c("x1", "x2"),
                            probs = c(0.05, 0.95))
  )
  res_95  <- suppressWarnings(
    credible_optimum_region(draws, factor_names = c("x1", "x2"),
                            probs = c(0.025, 0.975))
  )
  # Lower bound of 90% >= lower bound of 95%
  expect_true(all(res_90[["q5.0"]] >= res_95[["q2.5"]]))
  # Upper bound of 90% <= upper bound of 95%
  expect_true(all(res_90[["q95.0"]] <= res_95[["q97.5"]]))
})

test_that("credible_optimum_region errors on invalid probs", {
  draws <- .make_cor_draws()
  expect_error(
    credible_optimum_region(draws, factor_names = "x1", probs = c(-0.1, 0.9)),
    "probs"
  )
  expect_error(
    credible_optimum_region(draws, factor_names = "x1", probs = c(0.1, 1.1)),
    "probs"
  )
})
