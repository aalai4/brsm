context("Model term recommendation")

test_that("recommend_brsm_model_terms identifies linear data correctly", {
  # Pure linear relationship
  set.seed(123)
  dat_linear <- data.frame(
    x1 = runif(100, -1, 1),
    x2 = runif(100, -1, 1),
    y = 5 + 2*x1 + 3*x2 + rnorm(100, sd = 0.3)
  )

  result <- recommend_brsm_model_terms(
    data = dat_linear,
    response = "y",
    factor_names = c("x1", "x2"),
    criterion = "AIC",
    verbose = FALSE
  )

  expect_is(result, "brsm_term_recommendation")
  expect_equal(result$recommended_terms, "first_order")
  expect_true(result$fit_comparison$r_squared[1] > 0.95) # High R^2 for linear model
})

test_that("recommend_brsm_model_terms identifies quadratic data correctly", {
  # Quadratic relationship
  set.seed(456)
  dat_quad <- data.frame(
    x1 = runif(150, -1, 1),
    x2 = runif(150, -1, 1),
    y = 5 + 2*x1 - x1^2 + 0.5*x1*x2 + rnorm(150, sd = 0.5)
  )

  result <- recommend_brsm_model_terms(
    data = dat_quad,
    response = "y",
    factor_names = c("x1", "x2"),
    criterion = "AIC",
    verbose = FALSE
  )

  expect_is(result, "brsm_term_recommendation")
  expect_true(result$recommended_terms %in% c("second_order", "pure_quadratic"))
  expect_true(result$diagnostics$has_quadratic_value)
})

test_that("recommend_brsm_model_terms handles single factor", {
  set.seed(789)
  dat_single <- data.frame(
    x1 = seq(-1, 1, length.out = 50),
    y = 3 + 2*x1 + rnorm(50, sd = 0.2)
  )

  result <- recommend_brsm_model_terms(
    data = dat_single,
    response = "y",
    factor_names = "x1",
    criterion = "AIC"
  )

  expect_is(result, "brsm_term_recommendation")
  expect_equal(result$diagnostics$n_factors, 1L)
})

test_that("recommend_brsm_model_terms respects different criteria", {
  set.seed(111)
  dat <- data.frame(
    x1 = rnorm(80, 0, 1),
    x2 = rnorm(80, 0, 1),
    y = 5 + x1 + x2 + rnorm(80, sd = 0.5)
  )

  result_aic <- recommend_brsm_model_terms(
    data = dat, response = "y", factor_names = c("x1", "x2"),
    criterion = "AIC", verbose = FALSE
  )
  result_bic <- recommend_brsm_model_terms(
    data = dat, response = "y", factor_names = c("x1", "x2"),
    criterion = "BIC", verbose = FALSE
  )
  result_adjr2 <- recommend_brsm_model_terms(
    data = dat, response = "y", factor_names = c("x1", "x2"),
    criterion = "adjR2", verbose = FALSE
  )

  expect_is(result_aic, "brsm_term_recommendation")
  expect_is(result_bic, "brsm_term_recommendation")
  expect_is(result_adjr2, "brsm_term_recommendation")
})

test_that("recommend_brsm_model_terms validates input", {
  dat <- data.frame(x1 = rnorm(20), y = rnorm(20))

  expect_error(
    recommend_brsm_model_terms(data = "not_a_df", response = "y", factor_names = "x1"),
    "data must be a data.frame"
  )

  expect_error(
    recommend_brsm_model_terms(data = dat, response = c("y1", "y2"), factor_names = "x1"),
    "response must be a single character string"
  )

  expect_error(
    recommend_brsm_model_terms(data = dat, response = "missing", factor_names = "x1"),
    "response 'missing' not found"
  )

  expect_error(
    recommend_brsm_model_terms(data = dat, response = "y", factor_names = "x1", alpha = 2),
    "alpha must be a number between 0 and 1"
  )
})

test_that("recommend_brsm_model_terms fit_comparison has correct structure", {
  set.seed(222)
  dat <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    y = rnorm(50)
  )

  result <- recommend_brsm_model_terms(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    verbose = FALSE
  )

  expect_is(result$fit_comparison, "data.frame")
  expect_equal(nrow(result$fit_comparison), 4L) # 4 model types
  expected_cols <- c("model_terms", "n_params", "n_obs", "r_squared",
                     "adj_r_squared", "AIC", "BIC", "sigma")
  expect_true(all(expected_cols %in% names(result$fit_comparison)))
})

test_that("recommend_brsm_model_terms captures interaction effects", {
  # Strong interaction effect
  set.seed(333)
  dat_interact <- data.frame(
    x1 = runif(120, -1, 1),
    x2 = runif(120, -1, 1)
  )
  dat_interact$y <- 5 + x1 + x2 + 3*dat_interact$x1*dat_interact$x2 + rnorm(120, sd = 0.4)

  result <- recommend_brsm_model_terms(
    data = dat_interact,
    response = "y",
    factor_names = c("x1", "x2"),
    alpha = 0.05,
    verbose = FALSE
  )

  expect_true(result$diagnostics$has_interaction_value)
  expect_true(result$recommended_terms %in% c("first_order_twi", "second_order"))
})

test_that("recommend_brsm_model_terms print method works", {
  set.seed(444)
  dat <- data.frame(
    x1 = rnorm(40, 0, 1),
    y = 2*rnorm(40)
  )

  result <- recommend_brsm_model_terms(
    data = dat,
    response = "y",
    factor_names = "x1",
    verbose = FALSE
  )

  # Test that print method doesn't error
  expect_silent(print(result))
  expect_output(print(result), "BRSM Model Term Recommendation")
})

test_that("recommend_brsm_model_terms handles verbose output", {
  set.seed(555)
  dat <- data.frame(
    x1 = runif(60, -1, 1),
    x2 = runif(60, -1, 1),
    y = 3 + x1 + rnorm(60, sd = 0.5)
  )

  expect_output(
    recommend_brsm_model_terms(
      data = dat,
      response = "y",
      factor_names = c("x1", "x2"),
      verbose = TRUE
    ),
    "Model selection"
  )
})

test_that("recommend_brsm_model_terms diagnostics are reasonable", {
  set.seed(666)
  dat <- data.frame(
    x1 = rnorm(100, 0, 1),
    x2 = rnorm(100, 0, 1),
    y = 5 + 2*x1 + rnorm(100, sd = 0.5)
  )

  result <- recommend_brsm_model_terms(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    verbose = FALSE
  )

  diag <- result$diagnostics
  expect_is(diag$linear_rsq, "numeric")
  expect_true(diag$linear_rsq >= 0 && diag$linear_rsq <= 1)
  expect_is(diag$has_interaction_value, "logical")
  expect_is(diag$has_quadratic_value, "logical")
  expect_true(diag$mean_abs_residual >= 0)
})