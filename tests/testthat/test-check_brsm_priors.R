context("Prior predictive checks")

test_that("check_brsm_priors handles basic 1-factor case", {
  dat <- data.frame(
    x1 = seq(-1, 1, length.out = 30),
    y = rnorm(30, mean = 5, sd = 1)
  )

  result <- check_brsm_priors(
    data = dat,
    response = "y",
    factor_names = "x1",
    prior_profile = "regularized",
    model_terms = "first_order",
    n_prior_samples = 50L,
    n_grid_per_factor = 5L,
    plot = FALSE,
    seed = 123
  )

  expect_is(result, "list")
  expect_true("prior_samples" %in% names(result))
  expect_true("predictions" %in% names(result))
  expect_true("summary" %in% names(result))

  expect_is(result$prior_samples, "data.frame")
  expect_equal(nrow(result$prior_samples), 50L)

  expect_is(result$predictions, "data.frame")
  expect_true("pred" %in% names(result$predictions))

  expect_is(result$summary$stats, "list")
  expect_true(all(c("min", "max", "mean", "sd") %in% names(result$summary$stats)))
})

test_that("check_brsm_priors handles 2-factor case", {
  dat <- data.frame(
    x1 = runif(30, -1, 1),
    x2 = runif(30, -1, 1),
    y = rnorm(30, mean = 3, sd = 2)
  )

  result <- check_brsm_priors(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    prior_profile = "legacy",
    model_terms = "second_order",
    n_prior_samples = 50L,
    n_grid_per_factor = 4L,
    plot = FALSE,
    seed = 42
  )

  expect_is(result, "list")
  expect_equal(nrow(result$prior_samples), 50L)
  expect_true(nrow(result$predictions) > 0)
  expect_true(all(c("x1", "x2") %in% names(result$predictions)))
})

test_that("check_brsm_priors respects seed for reproducibility", {
  dat <- data.frame(
    x1 = rnorm(25, 0, 1),
    y = rnorm(25, 5, 1)
  )

  result1 <- check_brsm_priors(
    data = dat,
    response = "y",
    factor_names = "x1",
    model_terms = "first_order",
    n_prior_samples = 30L,
    n_grid_per_factor = 3L,
    plot = FALSE,
    seed = 999
  )

  result2 <- check_brsm_priors(
    data = dat,
    response = "y",
    factor_names = "x1",
    model_terms = "first_order",
    n_prior_samples = 30L,
    n_grid_per_factor = 3L,
    plot = FALSE,
    seed = 999
  )

  # With same seed, prior samples should be identical
  expect_equal(result1$prior_samples, result2$prior_samples)
})

test_that("check_brsm_priors validates input", {
  dat <- data.frame(x1 = rnorm(20), y = rnorm(20))

  expect_error(
    check_brsm_priors(data = "not_a_df", response = "y", factor_names = "x1"),
    "data must be a data.frame"
  )

  expect_error(
    check_brsm_priors(data = dat, response = c("y1", "y2"), factor_names = "x1"),
    "response must be a single character string"
  )

  expect_error(
    check_brsm_priors(data = dat, response = "missing", factor_names = "x1"),
    "response 'missing' not found"
  )

  expect_error(
    check_brsm_priors(data = dat, response = "y", factor_names = "x1",
                      n_prior_samples = 0L),
    "n_prior_samples must be >= 1"
  )

  expect_error(
    check_brsm_priors(data = dat, response = "y", factor_names = "x1",
                      n_grid_per_factor = 1L),
    "n_grid_per_factor must be >= 2"
  )
})

test_that("check_brsm_priors summary stats are computed correctly", {
  dat <- data.frame(
    x1 = seq(-1, 1, length.out = 20),
    y = seq(10, 20, length.out = 20)
  )

  result <- check_brsm_priors(
    data = dat,
    response = "y",
    factor_names = "x1",
    prior_profile = "regularized",
    model_terms = "first_order",
    n_prior_samples = 50L,
    n_grid_per_factor = 4L,
    plot = FALSE,
    seed = 123
  )

  stats <- result$summary$stats
  preds <- result$predictions$pred

  # Check that reported min/max match computed values
  expect_equal(stats$min, min(preds, na.rm = TRUE))
  expect_equal(stats$max, max(preds, na.rm = TRUE))
  expect_equal(stats$mean, mean(preds, na.rm = TRUE), tolerance = 0.001)

  # Check observed range is captured
  obs_range <- result$summary$obs_range
  expect_length(obs_range, 2)
  expect_true(obs_range[1] <= obs_range[2])
})

test_that("check_brsm_priors handles different prior profiles", {
  dat <- data.frame(
    x1 = rnorm(30, 0, 1),
    x2 = rnorm(30, 0, 1),
    y = rnorm(30, 5, 2)
  )

  for (profile in c("legacy", "regularized", "adaptive")) {
    result <- check_brsm_priors(
      data = dat,
      response = "y",
      factor_names = c("x1", "x2"),
      prior_profile = profile,
      model_terms = "second_order",
      n_prior_samples = 40L,
      n_grid_per_factor = 3L,
      plot = FALSE,
      seed = 111
    )

    expect_is(result, "list")
    expect_equal(nrow(result$prior_samples), 40L)
  }
})

test_that("check_brsm_priors handles custom prior", {
  if (!requireNamespace("brms", quietly = TRUE)) {
    skip("brms not installed")
  }

  dat <- data.frame(
    x1 = rnorm(20, 0, 1),
    y = rnorm(20, 5, 1)
  )

  custom_prior <- brms::set_prior("normal(0, 2)", class = "b") +
                  brms::set_prior("normal(5, 3)", class = "Intercept") +
                  brms::set_prior("exponential(1)", class = "sigma")

  result <- check_brsm_priors(
    data = dat,
    response = "y",
    factor_names = "x1",
    prior = custom_prior,
    model_terms = "first_order",
    n_prior_samples = 30L,
    n_grid_per_factor = 3L,
    plot = FALSE,
    seed = 456
  )

  expect_is(result, "list")
  expect_equal(nrow(result$prior_samples), 30L)
})