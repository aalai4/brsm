# Focused tests for remaining mid-range coverage gaps.

.make_minimal_brsm_fit <- function(model_terms = "second_order") {
  draws <- data.frame(
    b_Intercept = 1:3,
    b_x1 = 4:6,
    b_x2 = 7:9,
    "b_I(x1^2)" = 10:12,
    "b_I(x2^2)" = 13:15,
    "b_x1:x2" = 16:18,
    sigma = c(0.5, 0.6, 0.7),
    check.names = FALSE
  )

  fit_obj <- structure(
    list(
      data = draws,
      .draws_df = draws
    ),
    class = c("fake_brmsfit_pp", "brmsfit")
  )

  as.data.frame.fake_brmsfit_pp <<- function(x, ...) x$.draws_df
  registerS3method("as.data.frame", "fake_brmsfit_pp",
    as.data.frame.fake_brmsfit_pp, envir = asNamespace("base"))

  obj <- list(
    fit = fit_obj,
    factor_names = c("x1", "x2"),
    model_terms = model_terms,
    ranges = list(x1 = c(-1, 1), x2 = c(-1, 1))
  )
  class(obj) <- "brsm_fit"
  obj
}

.make_quad_draws <- function(n = 20, seed = 901) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 10, 1),
    b_x1 = rnorm(n, 0.8, 0.1),
    b_x2 = rnorm(n, -0.4, 0.1),
    "b_I(x1^2)" = rnorm(n, -0.3, 0.05),
    "b_I(x2^2)" = rnorm(n, -0.2, 0.05),
    "b_x1:x2" = rnorm(n, 0.03, 0.02),
    check.names = FALSE
  )
}

.make_linear_draws <- function(n = 15, seed = 902) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 5, 0.5),
    b_x1 = rnorm(n, 0.4, 0.1),
    b_x2 = rnorm(n, -0.3, 0.1),
    check.names = FALSE
  )
}

# ── classify_stationary_point ─────────────────────────────────────────────────

test_that("classify_stationary_point validates inputs and dispatches brsm_fit", {
  draws <- .make_quad_draws()
  result <- brsm::classify_stationary_point(draws, factor_names = c("x1", "x2"))
  expect_s3_class(result$classification, "factor")
  expect_true(all(levels(result$classification) %in%
    c("maximum", "minimum", "saddle", "indeterminate")))

  result_kappa <- brsm::classify_stationary_point(
    draws,
    factor_names = c("x1", "x2"),
    return_kappa = TRUE
  )
  expect_true("kappa" %in% names(result_kappa))

  expect_error(
    brsm::classify_stationary_point(draws, factor_names = c("x1", "x2"), tol = -1),
    "tol must be a non-negative"
  )
  expect_error(
    brsm::classify_stationary_point(draws, factor_names = c("x1", "x2"), return_kappa = NA),
    "return_kappa must be TRUE or FALSE"
  )

  fit_obj <- .make_minimal_brsm_fit()
  result_fit <- brsm::classify_stationary_point(fit_obj)
  expect_s3_class(result_fit, "data.frame")
  expect_true("classification" %in% names(result_fit))
})

test_that("classify_stationary_point covers indeterminate and warning branches", {
  expect_error(
    brsm::classify_stationary_point(.make_quad_draws(n = 3)),
    "factor_names must be supplied"
  )

  # Exactly zero curvature in one direction should produce indeterminate labels.
  draws_indeterminate <- data.frame(
    b_Intercept = 1:3,
    b_x1 = c(0.2, 0.3, 0.4),
    b_x2 = c(-0.1, -0.2, -0.3),
    `b_I(x1^2)` = c(0, 0, 0),
    `b_I(x2^2)` = c(-0.5, -0.5, -0.5),
    `b_x1:x2` = c(0, 0, 0),
    check.names = FALSE
  )
  out_ind <- brsm::classify_stationary_point(
    draws_indeterminate,
    factor_names = c("x1", "x2")
  )
  expect_true(any(as.character(out_ind$classification) == "indeterminate"))

  # Missing coefficients create NA Hessians, triggering the warning branch.
  draws_bad <- draws_indeterminate
  draws_bad$`b_I(x1^2)`[1] <- NA_real_
  expect_warning(
    brsm::classify_stationary_point(draws_bad, factor_names = c("x1", "x2")),
    "draws could not be classified"
  )
})

# ── plot_optimum_posterior ────────────────────────────────────────────────────

test_that("plot_optimum_posterior validates branches and renders plots", {
  skip_if_not_installed("ggplot2")

  draws <- .make_quad_draws(n = 20)
  stationary <- suppressWarnings(
    brsm::stationary_point(draws, factor_names = c("x1", "x2"))
  )
  stationary_df <- as.data.frame(stationary)
  names(stationary_df) <- c("x1", "x2")

  p <- brsm:::plot_optimum_posterior(
    draws = draws,
    factor_names = c("x1", "x2"),
    stationary_draws = stationary_df
  )
  expect_s3_class(p, "ggplot")

  p_points <- brsm:::plot_optimum_posterior(
    draws = draws,
    factor_names = c("x1", "x2"),
    stationary_draws = stationary_df,
    show_points = TRUE,
    point_alpha = 0.5
  )
  expect_equal(length(p_points$layers), 7)

  expect_error(
    brsm:::plot_optimum_posterior(draws, c("x1", "x2"), stationary_draws = stationary_df, bins = 0),
    "bins must be a positive integer"
  )
  expect_error(
    brsm:::plot_optimum_posterior(draws, c("x1", "x2"), stationary_draws = stationary_df, levels = 0),
    "levels must be a positive integer"
  )
  expect_error(
    brsm:::plot_optimum_posterior(draws, c("x1", "x2"), stationary_draws = stationary_df, alpha = 2),
    "alpha must be between 0 and 1"
  )
  expect_error(
    brsm:::plot_optimum_posterior(draws, c("x1", "x2"), stationary_draws = stationary_df, point_alpha = -1),
    "point_alpha must be between 0 and 1"
  )
  expect_error(
    brsm:::plot_optimum_posterior(draws, c("x1", "x2"), stationary_draws = stationary_df, seed = c(1, 2)),
    "seed must be a single numeric"
  )
  expect_error(
    suppressWarnings(brsm:::plot_optimum_posterior(draws, c("x1", "x2"),
      stationary_draws = data.frame(x1 = NA, x2 = NA)
    )),
    "no valid stationary points"
  )

  stationary_with_na <- rbind(stationary_df, c(NA, NA))
  expect_warning(
    brsm:::plot_optimum_posterior(
      draws, c("x1", "x2"), stationary_draws = stationary_with_na
    ),
    "stationary points removed due to missing values"
  )

  # Cover branch where stationary_draws are computed internally.
  p_auto <- suppressWarnings(brsm:::plot_optimum_posterior(
    draws = draws,
    factor_names = c("x1", "x2"),
    seed = 123
  ))
  expect_s3_class(p_auto, "ggplot")
})

# ── plot_ridge_path ───────────────────────────────────────────────────────────

test_that("plot_ridge_path validates and renders ridge path", {
  skip_if_not_installed("ggplot2")

  draws <- .make_quad_draws(n = 12)
  ridge_result <- suppressWarnings(
    brsm::posterior_ridge_analysis(
      draws,
      factor_names = c("x1", "x2"),
      radii = c(0, 0.5, 1),
      summary = FALSE
    )
  )

  expect_error(
    brsm:::plot_ridge_path(ridge_result, c("x1", "x2"), response = "missing_col"),
    "response must be supplied and correspond to a column"
  )

  non_numeric_response <- ridge_result
  non_numeric_response$radius <- as.character(non_numeric_response$radius)
  expect_error(
    brsm:::plot_ridge_path(non_numeric_response, c("x1", "x2"), response = "radius"),
    "response column must be numeric"
  )

  expect_error(
    brsm:::plot_ridge_path(ridge_result, c("x1", "x2"), response = "radius", alpha = 2),
    "alpha must be between 0 and 1"
  )

  p_path <- brsm:::plot_ridge_path(ridge_result, c("x1", "x2"), response = "radius")
  expect_s3_class(p_path, "ggplot")

  p_draws <- brsm:::plot_ridge_path(
    ridge_result, c("x1", "x2"),
    response = "radius", show_draws = TRUE
  )
  expect_true(length(p_draws$layers) > length(p_path$layers))

  na_ridge <- ridge_result
  na_ridge$x1[1] <- NA
  expect_warning(
    brsm:::plot_ridge_path(na_ridge, c("x1", "x2"), response = "radius"),
    "rows removed from ridge_draws due to missing values"
  )

  all_na <- ridge_result
  all_na$x1 <- NA
  all_na$x2 <- NA
  expect_error(
    suppressWarnings(brsm:::plot_ridge_path(all_na, c("x1", "x2"), response = "radius")),
    "no valid rows"
  )
})

# ── steepest_ascent ───────────────────────────────────────────────────────────

test_that("steepest_ascent validates all inputs", {
  draws <- .make_quad_draws(n = 8)
  fn <- c("x1", "x2")

  expect_error(
    brsm::steepest_ascent(draws),
    "factor_names must be supplied"
  )

  expect_error(
    brsm::steepest_ascent(draws, fn, start = c(x1 = 0, x2 = 0), step_size = -1),
    "step_size must be a positive"
  )
  expect_error(
    brsm::steepest_ascent(draws, fn, start = c(x1 = 0, x2 = 0), n_steps = 0),
    "n_steps must be a positive"
  )
  expect_error(
    brsm::steepest_ascent(draws, fn, start = c(x1 = 0, x2 = 0), tol = -1),
    "tol must be a non-negative"
  )
  expect_error(
    brsm::steepest_ascent(draws, fn, start = c(1, 2)),
    "start must be a named"
  )
  expect_error(
    brsm::steepest_ascent(draws, fn, start = c(x1 = 0)),
    "start must contain values for all factor_names"
  )
  expect_error(
    brsm::steepest_ascent(draws, fn, start = "bad"),
    "start must be a numeric"
  )

  result <- brsm::steepest_ascent(
    draws, fn,
    start = c(x1 = 0, x2 = 0),
    n_steps = 5, return_mean_path = TRUE
  )
  expect_type(result, "list")
  expect_true(all(c("directions", "paths", "mean_path") %in% names(result)))

  fit_obj <- .make_minimal_brsm_fit()
  result_fit <- brsm::steepest_ascent(fit_obj)
  expect_type(result_fit, "list")
  expect_true("directions" %in% names(result_fit))

  zero_grad_draws <- data.frame(
    b_Intercept = c(1, 1),
    b_x1 = c(0, 0),
    b_x2 = c(0, 0),
    `b_I(x1^2)` = c(0, 0),
    `b_I(x2^2)` = c(0, 0),
    `b_x1:x2` = c(0, 0),
    check.names = FALSE
  )
  expect_warning(
    out_zero <- brsm::steepest_ascent(
      zero_grad_draws,
      fn,
      start = c(x1 = 0, x2 = 0),
      drop_null_paths = FALSE
    ),
    "near-zero gradients"
  )
  expect_equal(nrow(out_zero$paths), 0)
  expect_equal(names(out_zero$paths), c("draw", "step", "x1", "x2"))

  one_valid_draw <- zero_grad_draws
  one_valid_draw$b_x1[1] <- 1
  expect_warning(
    out_one <- brsm::steepest_ascent(
      one_valid_draw,
      fn,
      start = c(x1 = 0, x2 = 0),
      n_steps = 3,
      return_mean_path = TRUE
    ),
    "near-zero gradients"
  )
  expect_false(is.null(out_one$mean_path))
  expect_equal(nrow(out_one$mean_path), 4)
})

# ── posterior_predict_brsm dispatch ──────────────────────────────────────────

test_that("posterior_predict_brsm dispatches through brsm_fit and brmsfit methods", {
  grid <- brsm::surface_grid(list(x1 = c(-1, 1), x2 = c(-1, 1)), n = 4)

  fit_obj <- .make_minimal_brsm_fit(model_terms = "second_order")

  out_fit <- brsm::posterior_predict_brsm(
    fit_obj,
    newdata = grid,
    summary = TRUE
  )
  expect_equal(nrow(out_fit), nrow(grid))
  expect_true("mean" %in% names(out_fit))

  linear_df <- data.frame(
    b_Intercept = 1:3,
    b_x1 = 4:6,
    b_x2 = 7:9,
    check.names = FALSE
  )
  fake_brmsfit <- structure(
    list(.draws_df = linear_df),
    class = c("fake_brmsfit_pp", "brmsfit")
  )

  out_brmsfit <- brsm::posterior_predict_brsm(
    fake_brmsfit,
    factor_names = c("x1", "x2"),
    newdata = grid,
    include_residual = FALSE,
    summary = TRUE
  )
  expect_equal(nrow(out_brmsfit), nrow(grid))
})

test_that("surface_grid covers remaining validation branches", {
  expect_error(
    brsm::surface_grid(ranges = list(x1 = 1, x2 = c(-1, 1)), n = 3),
    "range for factor 'x1' must be numeric of length 2"
  )

  expect_error(
    brsm::surface_grid(
      ranges = list(x1 = c(-1, 1), x2 = c(-1, 1)),
      n = 3,
      center = list(x1 = 0, x2 = 0)
    ),
    "center must be a numeric vector"
  )
})

test_that("stability summary fills missing status_count levels", {
  obj <- list(dummy = TRUE)
  attr(obj, "diagnostics") <- list(
    function_name = "custom_fn",
    status_code = c(0L, 0L, 1L),
    status_counts = c(ok = 2L)
  )

  out <- brsm::summarize_brsm_stability(obj)
  expect_equal(out$status_ok[[1]], 2)
  expect_equal(out$status_lapack_fail[[1]], 0)
  expect_equal(out$status_invalid_lu_diag[[1]], 0)
  expect_equal(out$status_kappa_exceeded[[1]], 0)
})
