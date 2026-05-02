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

test_that("credible_optimum_region dispatches through brsm_fit", {
  as.data.frame.fake_brmsfit_cor <<- function(x, ...) x$.draws_df
  registerS3method(
    "as.data.frame", "fake_brmsfit_cor",
    as.data.frame.fake_brmsfit_cor,
    envir = asNamespace("base")
  )

  draws <- .make_cor_draws(n = 25)
  fake_fit <- structure(
    list(
      fit = structure(list(.draws_df = draws), class = c("fake_brmsfit_cor", "brmsfit")),
      factor_names = c("x1", "x2"),
      model_terms = "second_order"
    ),
    class = "brsm_fit"
  )

  out <- suppressWarnings(credible_optimum_region(fake_fit))
  expect_s3_class(out, "data.frame")
  expect_equal(rownames(out), c("x1", "x2"))
})

test_that("credible_optimum_region handles stationary-point edge outputs", {
  draws <- .make_cor_draws(n = 12)
  ns <- asNamespace("brsm")

  .replace_stationary <- function(fun) {
    unlockBinding("stationary_point", ns)
    old <- get("stationary_point", envir = ns)
    assign("stationary_point", fun, envir = ns)
    lockBinding("stationary_point", ns)
    old
  }

  .restore_stationary <- function(old) {
    unlockBinding("stationary_point", ns)
    assign("stationary_point", old, envir = ns)
    lockBinding("stationary_point", ns)
  }

  old <- .replace_stationary(function(...) data.frame(x1 = numeric(0), x2 = numeric(0)))
  on.exit(.restore_stationary(old), add = TRUE)
  expect_error(
    credible_optimum_region(draws, factor_names = c("x1", "x2")),
    "zero rows"
  )

  .restore_stationary(old)
  old <- .replace_stationary(function(...) data.frame(a = 1:3, b = 2:4, c = 3:5))
  expect_error(
    credible_optimum_region(draws, factor_names = c("x1", "x2")),
    "unexpected number of columns"
  )

  .restore_stationary(old)
  old <- .replace_stationary(function(...) data.frame(x1 = rep(NA_real_, 4), x2 = rep(NA_real_, 4)))
  expect_warning(
    out_all_na <- credible_optimum_region(draws, factor_names = c("x1", "x2")),
    "All stationary points could not be computed"
  )
  expect_true(all(is.na(out_all_na$mean)))
  expect_true(all(is.na(out_all_na$sd)))

  .restore_stationary(old)
  old <- .replace_stationary(function(...) data.frame(
    x1 = c(NA, 0.1, 0.2, NA),
    x2 = c(0.0, 0.1, NA, 0.3)
  ))
  expect_warning(
    credible_optimum_region(draws, factor_names = c("x1", "x2")),
    "stationary points could not be computed"
  )
})
