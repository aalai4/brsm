# Focused coverage tests for helper-heavy branches.

.make_local_mock_brsm_fit <- function() {
  obj <- list(
    factor_names = c("x1", "x2"),
    fit = NULL
  )
  class(obj) <- "brsm_fit"
  obj
}

.create_fake_brmsfit <- function(fixed,
                                 nuts_params = NULL,
                                 metadata = NULL,
                                 draws = NULL) {
  fit_slot <- list()
  if (!is.null(metadata)) {
    fit_slot$metadata <- metadata
  }

  structure(
    list(
      .summary = list(fixed = fixed),
      .nuts_params = nuts_params,
      .draws = draws,
      fit = fit_slot
    ),
    class = c("fake_brmsfit", "brmsfit")
  )
}

summary.fake_brmsfit <- function(object, ...) {
  object$.summary
}

as.data.frame.fake_brmsfit <- function(x, ...) {
  x$.draws
}

registerS3method(
  "summary",
  "fake_brmsfit",
  summary.fake_brmsfit,
  envir = asNamespace("base")
)
registerS3method(
  "as.data.frame",
  "fake_brmsfit",
  as.data.frame.fake_brmsfit,
  envir = asNamespace("base")
)

.create_helper_draws_2d <- function(n = 30, seed = 802) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 10, 1),
    b_x1 = rnorm(n, 0.8, 0.1),
    b_x2 = rnorm(n, -0.4, 0.1),
    "b_I(x1^2)" = rnorm(n, -0.2, 0.05),
    "b_I(x2^2)" = rnorm(n, -0.1, 0.05),
    "b_x1:x2" = rnorm(n, 0.03, 0.02),
    b_extra = rnorm(n, 0, 0.01),
    check.names = FALSE
  )
}

.create_helper_draws_3d <- function(n = 20, seed = 803) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 5, 1),
    b_x1 = rnorm(n, 0.4, 0.1),
    b_x2 = rnorm(n, -0.3, 0.1),
    b_x3 = rnorm(n, 0.2, 0.1),
    check.names = FALSE
  )
}

test_that("surface_grid sorts ranges and handles center insertion", {
  grid <- brsm::surface_grid(
    ranges = list(x1 = c(1, -1), x2 = c(2, 0)),
    n = 3,
    center = c(x1 = 0.5, x2 = 1)
  )

  expect_equal(nrow(grid), 10)
  expect_equal(unname(as.numeric(grid[1, c("x1", "x2")])), c(0.5, 1))
  expect_equal(range(grid$x1), c(-1, 1))
  expect_equal(range(grid$x2), c(0, 2))

  deduped <- brsm::surface_grid(
    ranges = list(x1 = c(-1, 1), x2 = c(0, 2)),
    n = 3,
    center = c(x1 = 0, x2 = 1)
  )
  expect_equal(nrow(deduped), 9)
})

test_that("surface_grid validates ranges, n, and center inputs", {
  expect_error(brsm::surface_grid(NULL), "ranges must be supplied")
  expect_error(brsm::surface_grid(c(-1, 1)), "ranges must be a named list")

  unnamed_ranges <- list(c(-1, 1))
  names(unnamed_ranges) <- NULL
  expect_error(brsm::surface_grid(unnamed_ranges), "factor names")

  expect_error(
    brsm::surface_grid(list(x1 = c(-1, 1)), n = 1),
    "n must be a finite integer >= 2"
  )
  expect_error(
    brsm::surface_grid(list(x1 = c(-1, Inf))),
    "must contain finite values"
  )
  expect_error(
    brsm::surface_grid(list(x1 = c(-1, 1)), center = 1),
    "center must be a named numeric vector"
  )
  expect_error(
    brsm::surface_grid(
      list(x1 = c(-1, 1), x2 = c(-1, 1)),
      center = c(x1 = 0)
    ),
    "center must contain values for all factors"
  )
  expect_error(
    brsm::surface_grid(
      list(x1 = c(-1, 1), x2 = c(-1, 1)),
      center = c(x1 = 0, x2 = 0, x3 = 0)
    ),
    "same length"
  )
})

test_that("as_brsm_draws normalizes data frame coefficient names", {
  raw <- .create_helper_draws_2d(n = 8)
  names(raw)[names(raw) == "b_I(x1^2)"] <- "b_Ix1E2"
  names(raw)[names(raw) == "b_I(x2^2)"] <- "b_Ix2E2"
  names(raw)[names(raw) == "b_x1:x2"] <- "b_x1.x2"

  draws <- brsm::as_brsm_draws(raw, factor_names = c("x1", "x2"))

  expect_true(all(c(
    "b_Intercept",
    "b_x1",
    "b_x2",
    "b_I(x1^2)",
    "b_I(x2^2)",
    "b_x1:x2",
    "b_extra"
  ) %in% names(draws)))
})

test_that("as_brsm_draws validates fake brmsfit and dispatches by model terms", {
  fake_linear <- .create_fake_brmsfit(
    fixed = data.frame(Estimate = 1, row.names = "b_x1"),
    draws = data.frame(
      b_Intercept = 1:3,
      b_x1 = 4:6,
      b_x2 = 7:9,
      check.names = FALSE
    )
  )

  expect_error(as_brsm_draws(fake_linear), "factor_names must be supplied")
  expect_error(
    as_brsm_draws(fake_linear, factor_names = c("x1", "x2")),
    "missing required posterior coefficient columns"
  )

  first_order_fit <- .make_local_mock_brsm_fit()
  first_order_fit$fit <- fake_linear
  first_order_fit$model_terms <- "first_order"
  out_first <- brsm::as_brsm_draws(first_order_fit)
  expect_true(all(c("b_Intercept", "b_x1", "b_x2") %in% names(out_first)))

  twi_fit <- .make_local_mock_brsm_fit()
  twi_fit$fit <- .create_fake_brmsfit(
    fixed = data.frame(Estimate = 1, row.names = "b_x1"),
    draws = data.frame(
      b_Intercept = 1:3,
      b_x1 = 4:6,
      b_x2 = 7:9,
      "b_x1:x2" = 10:12,
      check.names = FALSE
    )
  )
  twi_fit$model_terms <- "first_order_twi"
  out_twi <- brsm::as_brsm_draws(twi_fit)
  expect_true("b_x1:x2" %in% names(out_twi))

  quad_fit <- .make_local_mock_brsm_fit()
  quad_fit$fit <- .create_fake_brmsfit(
    fixed = data.frame(Estimate = 1, row.names = "b_x1"),
    draws = data.frame(
      b_Intercept = 1:3,
      b_x1 = 4:6,
      b_x2 = 7:9,
      "b_I(x1^2)" = 10:12,
      "b_I(x2^2)" = 13:15,
      check.names = FALSE
    )
  )
  quad_fit$model_terms <- "pure_quadratic"
  out_quad <- brsm::as_brsm_draws(quad_fit)
  expect_true(all(c("b_I(x1^2)", "b_I(x2^2)") %in% names(out_quad)))

  expect_error(
    brsm::as_brsm_draws(structure(list(), class = "not_supported"), c("x1", "x2")),
    "does not support objects of class"
  )
})

test_that("as_brsm_draws covers remaining edge branches", {
  # Exercise generic UseMethod line via direct generic call.
  one_factor <- data.frame(
    b_Intercept = 1:3,
    b_x1 = 4:6,
    `b_I(x1^2)` = 7:9,
    check.names = FALSE
  )
  out_one <- brsm::as_brsm_draws(one_factor, factor_names = "x1")
  expect_true(all(c("b_Intercept", "b_x1", "b_I(x1^2)") %in% names(out_one)))

  # Invalid brsm_fit path.
  bad_fit <- structure(list(fit = NULL, factor_names = c("x1", "x2")), class = "brsm_fit")
  expect_error(
    brsm::as_brsm_draws(bad_fit),
    "must contain a valid brmsfit model"
  )

  # Missing interaction when require_interactions = TRUE.
  missing_int <- data.frame(
    b_Intercept = 1:3,
    b_x1 = 4:6,
    b_x2 = 7:9,
    `b_I(x1^2)` = 10:12,
    `b_I(x2^2)` = 13:15,
    check.names = FALSE
  )
  expect_error(
    brsm:::as_brsm_draws.data.frame(
      missing_int,
      factor_names = c("x1", "x2"),
      require_quadratic = TRUE,
      require_interactions = TRUE
    ),
    "b_x1:x2"
  )

  # model_terms NULL fallback in brsm_fit method.
  fallback_fit <- .make_local_mock_brsm_fit()
  fallback_fit$fit <- .create_fake_brmsfit(
    fixed = data.frame(Estimate = 1, row.names = "b_x1"),
    draws = data.frame(
      b_Intercept = 1:3,
      b_x1 = 4:6,
      b_x2 = 7:9,
      `b_I(x1^2)` = 10:12,
      `b_I(x2^2)` = 13:15,
      `b_x1:x2` = 16:18,
      check.names = FALSE
    )
  )
  fallback_fit$model_terms <- NULL
  out_fallback <- brsm::as_brsm_draws(fallback_fit)
  expect_true(all(c("b_I(x1^2)", "b_I(x2^2)", "b_x1:x2") %in% names(out_fallback)))
})


test_that("plot_posterior_contours covers uncertainty and quantile surfaces", {
  skip_if_not_installed("ggplot2")

  draws <- brsm::as_brsm_draws(
    .create_helper_draws_2d(),
    factor_names = c("x1", "x2")
  )
  ranges <- list(x1 = c(-2, 2), x2 = c(-2, 2))

  p_uncertainty <- brsm::plot_posterior_contours(
    draws = draws,
    factor_names = c("x1", "x2"),
    ranges = ranges,
    type = "uncertainty",
    probs = c(0.1, 0.9)
  )
  p_quantile <- brsm::plot_posterior_contours(
    draws = draws,
    factor_names = c("x1", "x2"),
    ranges = ranges,
    type = "quantile",
    quantile = 0.9
  )

  expect_s3_class(p_uncertainty, "ggplot")
  expect_identical(p_uncertainty$labels$title, "Posterior Uncertainty Surface")
  expect_identical(p_quantile$labels$title, "Posterior Quantile Surface (0.9)")
})

test_that("plot_posterior_contours handles conditioning and overlays", {
  skip_if_not_installed("ggplot2")

  draws_3d <- brsm:::as_brsm_draws.data.frame(
    .create_helper_draws_3d(),
    factor_names = c("x1", "x2", "x3"),
    require_quadratic = FALSE,
    require_interactions = FALSE
  )
  ranges_3d <- list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2))

  p_user <- brsm::plot_posterior_contours(
    draws = draws_3d,
    factor_names = c("x1", "x2", "x3"),
    ranges = ranges_3d,
    vary_factors = c("x1", "x2"),
    conditioning = "user",
    fixed = c(x3 = 0.25)
  )
  expect_true(all(abs(p_user$data$x3 - 0.25) < 1e-8))

  expect_warning(
    p_fallback <- brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      vary_factors = c("x1", "x2"),
      conditioning = "optimum_mean"
    ),
    "falling back to center conditioning"
  )
  expect_true(all(abs(p_fallback$data$x3) < 1e-8))

  draws_2d <- brsm::as_brsm_draws(
    .create_helper_draws_2d(),
    factor_names = c("x1", "x2")
  )
  stationary_mean <- data.frame(x1 = c(-0.5, 0.5, NA), x2 = c(0.2, -0.2, NA))
  p_mean <- brsm::plot_posterior_contours(
    draws = draws_2d,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    overlay_stationary = TRUE,
    stationary_draws = stationary_mean,
    overlay_type = "mean"
  )
  expect_equal(length(p_mean$layers), 2)

  stationary_posterior <- data.frame(
    x1 = seq(-0.5, 0.5, length.out = 6),
    x2 = seq(0.5, -0.5, length.out = 6)
  )
  p_posterior <- brsm::plot_posterior_contours(
    draws = draws_2d,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    overlay_stationary = TRUE,
    stationary_draws = stationary_posterior,
    overlay_type = "posterior",
    overlay_max_draws = 3,
    seed = 11
  )
  built <- ggplot2::ggplot_build(p_posterior)
  expect_equal(nrow(built$data[[2]]), 3)

  expect_warning(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1", "x2"),
      ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
      pairwise = TRUE,
      overlay_stationary = TRUE
    ),
    "Stationary overlays are only available"
  )
})

test_that("plot_posterior_contours validates branch-specific inputs", {
  skip_if_not_installed("ggplot2")

  draws_2d <- brsm::as_brsm_draws(
    .create_helper_draws_2d(),
    factor_names = c("x1", "x2")
  )
  ranges_2d <- list(x1 = c(-2, 2), x2 = c(-2, 2))
  draws_3d <- brsm:::as_brsm_draws.data.frame(
    .create_helper_draws_3d(),
    factor_names = c("x1", "x2", "x3"),
    require_quadratic = FALSE,
    require_interactions = FALSE
  )
  ranges_3d <- list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2))

  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1", "x2"),
      ranges = ranges_2d,
      quantile = 1.5
    ),
    "quantile must be between 0 and 1"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1", "x2"),
      ranges = ranges_2d,
      bins = 0
    ),
    "bins must be a positive integer"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      vary_factors = c("x1", "x2"),
      slice = list(x1 = 0)
    ),
    "slice factor must differ from vary_factors"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1", "x2"),
      ranges = ranges_2d,
      overlay_stationary = TRUE,
      overlay_alpha = 2
    ),
    "overlay_alpha must be between 0 and 1"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1", "x2"),
      ranges = ranges_2d,
      overlay_stationary = TRUE,
      overlay_max_draws = 0
    ),
    "overlay_max_draws must be a positive integer"
  )
})

test_that("plot_posterior_contours covers remaining validation branches", {
  skip_if_not_installed("ggplot2")

  draws_2d <- brsm::as_brsm_draws(
    .create_helper_draws_2d(),
    factor_names = c("x1", "x2")
  )
  ranges_2d <- list(x1 = c(-2, 2), x2 = c(-2, 2))

  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1"),
      ranges = list(x1 = c(-1, 1))
    ),
    "at least two factors"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1", "x2"),
      ranges = NULL
    ),
    "ranges must be a named list"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1", "x2"),
      ranges = list(x1 = c(-1, 1), x2 = 0)
    ),
    "numeric vector of length 2"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1", "x2"),
      ranges = ranges_2d,
      n = 0
    ),
    "n must be a positive integer"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1", "x2"),
      ranges = ranges_2d,
      pairwise = c(TRUE, FALSE)
    ),
    "pairwise must be a single logical"
  )

  draws_3d <- brsm:::as_brsm_draws.data.frame(
    .create_helper_draws_3d(),
    factor_names = c("x1", "x2", "x3"),
    require_quadratic = FALSE,
    require_interactions = FALSE
  )
  ranges_3d <- list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2))

  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      vary_factors = c("x1", "x9")
    ),
    "subset of factor_names"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      slice = list(c(-1, 1))
    ),
    "named list"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      slice = list(z = c(-1, 1))
    ),
    "slice factor must be present"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      slice = list(x3 = c(NA_real_))
    ),
    "non-empty finite numeric vector"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      vary_factors = c("x1", "x2"),
      conditioning = "user",
      fixed = c(0.1)
    ),
    "fixed must be a named numeric"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      vary_factors = c("x1", "x2"),
      conditioning = "user",
      fixed = c(x9 = 0.1)
    ),
    "fixed is missing values"
  )
})

test_that("plot_posterior_contours covers mocked prediction and optimum branches", {
  skip_if_not_installed("ggplot2")

  draws_3d <- brsm:::as_brsm_draws.data.frame(
    .create_helper_draws_3d(),
    factor_names = c("x1", "x2", "x3"),
    require_quadratic = FALSE,
    require_interactions = FALSE
  )
  ranges_3d <- list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2))

  ns <- asNamespace("brsm")

  unlockBinding("stationary_point", ns)
  old_sp <- get("stationary_point", envir = ns)
  assign("stationary_point", function(...) {
    data.frame(x1 = c(0.1, 0.2), x2 = c(-0.1, 0.2), x3 = c(NA_real_, NA_real_))
  }, envir = ns)
  lockBinding("stationary_point", ns)

  on.exit({
    unlockBinding("stationary_point", ns)
    assign("stationary_point", old_sp, envir = ns)
    lockBinding("stationary_point", ns)
  }, add = TRUE)

  expect_warning(
    p_opt <- brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      vary_factors = c("x1", "x2"),
      conditioning = "optimum_mean"
    ),
    "Optimum mean unavailable"
  )
  expect_s3_class(p_opt, "ggplot")

  unlockBinding("predict_surface", ns)
  old_ps <- get("predict_surface", envir = ns)
  assign("predict_surface", function(draws, factor_names, newdata, summary, probs) {
    data.frame(newdata, mean = 0)
  }, envir = ns)
  lockBinding("predict_surface", ns)
  on.exit({
    unlockBinding("predict_surface", ns)
    assign("predict_surface", old_ps, envir = ns)
    lockBinding("predict_surface", ns)
  }, add = TRUE)

  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      type = "uncertainty"
    ),
    "expected quantile columns"
  )
  expect_error(
    brsm::plot_posterior_contours(
      draws = draws_3d,
      factor_names = c("x1", "x2", "x3"),
      ranges = ranges_3d,
      type = "quantile",
      quantile = 0.9
    ),
    "requested quantile not found"
  )
})

test_that("plot_posterior_contours uses optimum_mean conditioning when available", {
  skip_if_not_installed("ggplot2")

  draws_3d <- brsm:::as_brsm_draws.data.frame(
    .create_helper_draws_3d(),
    factor_names = c("x1", "x2", "x3"),
    require_quadratic = FALSE,
    require_interactions = FALSE
  )
  ranges_3d <- list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2))

  ns <- asNamespace("brsm")
  unlockBinding("stationary_point", ns)
  old_sp <- get("stationary_point", envir = ns)
  assign("stationary_point", function(...) {
    data.frame(
      x1 = c(0.1, 0.2),
      x2 = c(-0.1, 0.2),
      x3 = c(0.3, 0.5)
    )
  }, envir = ns)
  lockBinding("stationary_point", ns)
  on.exit({
    unlockBinding("stationary_point", ns)
    assign("stationary_point", old_sp, envir = ns)
    lockBinding("stationary_point", ns)
  }, add = TRUE)

  p <- brsm::plot_posterior_contours(
    draws = draws_3d,
    factor_names = c("x1", "x2", "x3"),
    ranges = ranges_3d,
    vary_factors = c("x1", "x2"),
    conditioning = "optimum_mean"
  )
  expect_s3_class(p, "ggplot")
  expect_true(all(abs(p$data$x3 - mean(c(0.3, 0.5))) < 1e-8))
})

test_that("plot_posterior_contours computes/validates stationary overlays internally", {
  skip_if_not_installed("ggplot2")

  draws_2d <- brsm::as_brsm_draws(
    .create_helper_draws_2d(),
    factor_names = c("x1", "x2")
  )

  p_overlay <- suppressWarnings(brsm::plot_posterior_contours(
    draws = draws_2d,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    overlay_stationary = TRUE,
    overlay_type = "mean"
  ))
  expect_s3_class(p_overlay, "ggplot")

  expect_warning(
    brsm::plot_posterior_contours(
      draws = draws_2d,
      factor_names = c("x1", "x2"),
      ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
      overlay_stationary = TRUE,
      stationary_draws = data.frame(x1 = NA_real_, x2 = NA_real_),
      overlay_type = "posterior"
    ),
    "No valid stationary points"
  )
})