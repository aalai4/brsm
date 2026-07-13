# Tests for fit_brsm exact conjugate sampler backend

test_that("fit_brsm exact conjugate sampler fits and matches expected output", {
  dat <- generate_simulation_data(n = 50, seed = 123)

  # Fit using exact conjugate sampler
  fit <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    draws = 1000,
    seed = 123,
    coding_policy = "ignore"
  )

  expect_s3_class(fit, "brsm_fit")
  expect_s3_class(fit$fit, "brsm_conjugate_fit")

  # Draws conversion
  draws <- as_brsm_draws(fit)
  expect_true(is.data.frame(draws))
  expect_equal(nrow(draws), 1000)
  expect_true(all(c("b_Intercept", "b_x1", "b_x2", "b_x1:x2", "b_I(x1^2)", "b_I(x2^2)", "sigma") %in% names(draws)))

  # Downstream workflow checks
  # 1. Stationary Point
  sp <- stationary_point(fit)
  expect_true(is.matrix(sp) || is.data.frame(sp))
  expect_equal(ncol(sp), 2)
  expect_equal(colnames(sp), c("x1", "x2"))

  # 2. Curvature classification
  cls <- classify_stationary_point(fit)
  expect_s3_class(cls, "data.frame")
  expect_true("classification" %in% names(cls))

  # 3. Prediction surface
  grid <- surface_grid(fit$ranges, n = 5)
  pred <- predict_surface(fit, newdata = grid, summary = TRUE)
  expect_s3_class(pred, "data.frame")
  expect_true("mean" %in% names(pred))

  # 4. Posterior predictive check
  ppc <- check_brsm_ppc(fit, ndraws = 200)
  expect_s3_class(ppc$summary, "data.frame")
  expect_equal(ppc$summary$n_obs, 50)
  expect_equal(ppc$summary$ndraws, 200)

  # Compare posterior mean with classical OLS to ensure it aligns
  # (since priors are weakly informative, they should be very close)
  ols <- lm(y ~ x1 + x2 + x1:x2 + I(x1^2) + I(x2^2), data = dat)
  ols_coefs <- coef(ols)
  names(ols_coefs) <- c("b_Intercept", "b_x1", "b_x2", "b_x1:x2", "b_I(x1^2)", "b_I(x2^2)")
  
  post_means <- colMeans(draws[, names(ols_coefs)])
  
  for (name in names(ols_coefs)) {
    expect_equal(post_means[[name]], ols_coefs[[name]], tolerance = 0.5)
  }
})
