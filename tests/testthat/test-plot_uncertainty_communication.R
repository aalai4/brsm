test_that("plot_posterior_contours uncertainty type has explicit uncertainty legend", {
  skip_if_not_installed("ggplot2")

  set.seed(901)
  draws_raw <- data.frame(
    b_Intercept = rnorm(150, 5, 0.2),
    b_x1 = rnorm(150, 2, 0.15),
    b_x2 = rnorm(150, 3, 0.15),
    "b_I(x1^2)" = rnorm(150, -1.2, 0.08),
    "b_I(x2^2)" = rnorm(150, -1.5, 0.08),
    "b_x1:x2" = rnorm(150, 0.5, 0.08),
    check.names = FALSE
  )

  draws <- as_brsm_draws(draws_raw, factor_names = c("x1", "x2"))

  p <- plot_posterior_contours(
    draws = draws,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    type = "uncertainty",
    probs = c(0.1, 0.9)
  )

  expect_s3_class(p, "ggplot")
  expect_true(grepl("PI width", p$labels$fill, fixed = TRUE))
  expect_true(grepl("higher posterior predictive uncertainty", p$labels$subtitle, fixed = TRUE))
})

test_that("plot_posterior_contours can overlay uncertainty on mean surface", {
  skip_if_not_installed("ggplot2")

  set.seed(902)
  draws_raw <- data.frame(
    b_Intercept = rnorm(120, 4.8, 0.25),
    b_x1 = rnorm(120, 1.7, 0.12),
    b_x2 = rnorm(120, 2.4, 0.12),
    "b_I(x1^2)" = rnorm(120, -0.9, 0.07),
    "b_I(x2^2)" = rnorm(120, -1.1, 0.07),
    "b_x1:x2" = rnorm(120, 0.4, 0.06),
    check.names = FALSE
  )

  draws <- as_brsm_draws(draws_raw, factor_names = c("x1", "x2"))

  p <- plot_posterior_contours(
    draws = draws,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-1.5, 1.5), x2 = c(-1.5, 1.5)),
    type = "mean",
    probs = c(0.1, 0.9),
    uncertainty_overlay = TRUE
  )

  expect_s3_class(p, "ggplot")
  expect_true(grepl("Dashed contours show predictive interval width", p$labels$subtitle, fixed = TRUE))
  expect_gte(length(p$layers), 2)
})

test_that("plot_optimum_posterior reports uncertainty contour levels", {
  skip_if_not_installed("ggplot2")

  set.seed(903)
  draws_raw <- data.frame(
    b_Intercept = rnorm(160, 5, 0.2),
    b_x1 = rnorm(160, 2, 0.1),
    b_x2 = rnorm(160, 3, 0.1),
    "b_I(x1^2)" = rnorm(160, -1.1, 0.05),
    "b_I(x2^2)" = rnorm(160, -1.4, 0.05),
    "b_x1:x2" = rnorm(160, 0.45, 0.05),
    check.names = FALSE
  )

  draws <- as_brsm_draws(draws_raw, factor_names = c("x1", "x2"))

  p <- brsm:::plot_optimum_posterior(
    draws = draws,
    factor_names = c("x1", "x2"),
    show_uncertainty_contours = TRUE,
    uncertainty_levels = c(0.5, 0.9)
  )

  expect_s3_class(p, "ggplot")
  expect_true(grepl("Dashed ellipses", p$labels$subtitle, fixed = TRUE))
  expect_true(grepl("50%", p$labels$subtitle, fixed = TRUE))
  expect_true(grepl("90%", p$labels$subtitle, fixed = TRUE))
})

test_that("plot_ridge_path shows interval communication in subtitle", {
  skip_if_not_installed("ggplot2")

  set.seed(904)
  levels <- rep(seq(0, 1, length.out = 6), each = 40)
  ridge_draws <- data.frame(
    x1 = rnorm(length(levels), mean = rep(seq(-1, 1, length.out = 6), each = 40), sd = 0.2),
    x2 = rnorm(length(levels), mean = rep(seq(-0.8, 1.2, length.out = 6), each = 40), sd = 0.2),
    ridge_level = levels
  )

  p <- brsm:::plot_ridge_path(
    ridge_draws = ridge_draws,
    factor_names = c("x1", "x2"),
    response = "ridge_level",
    show_interval_crossbars = TRUE,
    interval_prob = 0.8
  )

  expect_s3_class(p, "ggplot")
  expect_true(grepl("crossbars show 80% pointwise posterior intervals", p$labels$subtitle, fixed = TRUE))
})