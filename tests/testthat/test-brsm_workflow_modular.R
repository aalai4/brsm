# Tests for modular brsm_workflow behavior

.create_workflow_draws <- function(n = 80, seed = 777) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 10, 1),
    b_x1 = rnorm(n, 1, 0.2),
    b_x2 = rnorm(n, -0.5, 0.2),
    "b_I(x1^2)" = rnorm(n, -0.2, 0.05),
    "b_I(x2^2)" = rnorm(n, -0.1, 0.05),
    "b_x1:x2" = rnorm(n, 0.05, 0.02),
    check.names = FALSE
  )
}

.create_workflow_draws_3d <- function(n = 80, seed = 784) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 10, 1),
    b_x1 = rnorm(n, 0.8, 0.2),
    b_x2 = rnorm(n, -0.5, 0.2),
    b_x3 = rnorm(n, 0.3, 0.2),
    check.names = FALSE
  )
}

.create_workflow_draws_3d_quadratic <- function(n = 80, seed = 788) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 10, 1),
    b_x1 = rnorm(n, 0.8, 0.2),
    b_x2 = rnorm(n, -0.5, 0.2),
    b_x3 = rnorm(n, 0.3, 0.2),
    "b_I(x1^2)" = rnorm(n, -0.15, 0.05),
    "b_I(x2^2)" = rnorm(n, -0.10, 0.05),
    "b_I(x3^2)" = rnorm(n, -0.08, 0.05),
    "b_x1:x2" = rnorm(n, 0.04, 0.02),
    "b_x1:x3" = rnorm(n, -0.03, 0.02),
    "b_x2:x3" = rnorm(n, 0.02, 0.02),
    check.names = FALSE
  )
}

test_that("brsm_workflow steps can skip ridge and predictions", {
  draws <- .create_workflow_draws()

  out <- brsm_workflow(
    object = draws,
    factor_names = c("x1", "x2"),
    steps = c("stationary", "classification", "steepest_ascent"),
    include_plots = FALSE
  )

  expect_true("draws" %in% names(out))
  expect_true("stationary_points" %in% names(out))
  expect_true("classification" %in% names(out))
  expect_true("steepest_ascent" %in% names(out))
  expect_false("ridge_analysis" %in% names(out))
  expect_false("predictions" %in% names(out))
  expect_false("grid" %in% names(out))
})

test_that("brsm_workflow predictions step requires ranges", {
  draws <- .create_workflow_draws()

  expect_error(
    brsm_workflow(
      object = draws,
      factor_names = c("x1", "x2"),
      steps = c("predictions"),
      include_plots = FALSE
    ),
    "ranges must be supplied"
  )

  out <- brsm_workflow(
    object = draws,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    steps = c("predictions"),
    include_plots = FALSE
  )

  expect_true("grid" %in% names(out))
  expect_true("predictions" %in% names(out))
  expect_true("surface" %in% names(out))
})

test_that("brsm_workflow fit_mode forwards fit args", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 778)

  out <- brsm_workflow(
    object = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    fit_mode = "fit_brsm",
    fit_args = list(
      chains = 1,
      iter = 250,
      warmup = 125,
      seed = 778,
      sampling_preset = "fast",
      refresh = 0,
      silent = 2,
      prior = brms::prior(normal(0, 1), class = "b")
    ),
    steps = c("stationary"),
    include_plots = FALSE
  )

  expect_true("fit" %in% names(out))
  expect_s3_class(out$fit, "brsm_fit")
  expect_true("stationary_points" %in% names(out))
})

test_that("brsm_workflow supports first-order draws for prediction", {
  set.seed(779)
  draws <- data.frame(
    b_Intercept = rnorm(40, 10, 1),
    b_x1 = rnorm(40, 1, 0.2),
    b_x2 = rnorm(40, -0.5, 0.2),
    check.names = FALSE
  )

  out <- brsm_workflow(
    object = draws,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    steps = c("predictions", "steepest_ascent"),
    include_plots = FALSE
  )

  expect_true("predictions" %in% names(out))
  expect_true("steepest_ascent" %in% names(out))
  expect_false("stationary_points" %in% names(out))
})

test_that("brsm_workflow skips geometry steps without quadratics", {
  set.seed(780)
  draws <- data.frame(
    b_Intercept = rnorm(30, 10, 1),
    b_x1 = rnorm(30, 1, 0.2),
    b_x2 = rnorm(30, -0.5, 0.2),
    check.names = FALSE
  )

  expect_warning(
    out <- brsm_workflow(
      object = draws,
      factor_names = c("x1", "x2"),
      ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
      steps = c("stationary", "classification", "ridge", "predictions"),
      include_plots = FALSE
    ),
    "Skipping geometry-dependent step"
  )

  expect_true("predictions" %in% names(out))
  expect_false("stationary_points" %in% names(out))
  expect_false("classification" %in% names(out))
  expect_false("ridge_analysis" %in% names(out))
})

test_that("brsm_workflow supports fit_mode first_order models", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 781)

  out <- brsm_workflow(
    object = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    fit_mode = "fit_brsm",
    fit_args = list(
      model_terms = "first_order",
      chains = 1,
      iter = 250,
      warmup = 125,
      seed = 781,
      sampling_preset = "fast",
      refresh = 0,
      silent = 2
    ),
    steps = c("predictions", "steepest_ascent", "stationary"),
    include_plots = FALSE
  )

  expect_s3_class(out$fit, "brsm_fit")
  expect_equal(out$fit$model_terms, "first_order")
  expect_true("predictions" %in% names(out))
  expect_true("steepest_ascent" %in% names(out))
  expect_false("stationary_points" %in% names(out))
})

test_that("brsm_workflow uses direction-only plot mode for first-order draws", {
  set.seed(782)
  draws <- data.frame(
    b_Intercept = rnorm(40, 10, 1),
    b_x1 = rnorm(40, 1, 0.2),
    b_x2 = rnorm(40, -0.5, 0.2),
    check.names = FALSE
  )

  out <- suppressWarnings(brsm_workflow(
    object = draws,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    steps = c("plots"),
    include_plots = TRUE
  ))

  expect_true("plots" %in% names(out))
  expect_equal(out$plots$plot_mode, "direction_only")
  expect_false(out$plots$geometry_available)
  expect_true("contours" %in% names(out$plots))
  expect_true("direction" %in% names(out$plots))
  expect_false("optimum" %in% names(out$plots))
})

test_that("brsm_workflow marks optimization plot mode for quadratic draws", {
  draws <- .create_workflow_draws(n = 40, seed = 783)

  out <- brsm_workflow(
    object = draws,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    steps = c("plots"),
    include_plots = TRUE
  )

  expect_true("plots" %in% names(out))
  expect_equal(out$plots$plot_mode, "optimization_geometry")
  expect_true(out$plots$geometry_available)
  expect_true("optimum" %in% names(out$plots))
})

test_that("plot_posterior_contours supports pairwise mode for 3 factors", {
  skip_if_not_installed("ggplot2")

  draws <- .create_workflow_draws_3d(n = 30, seed = 785)

  p <- plot_posterior_contours(
    draws = as_brsm_draws(draws, c("x1", "x2", "x3"),
      require_quadratic = FALSE,
      require_interactions = FALSE
    ),
    factor_names = c("x1", "x2", "x3"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2)),
    pairwise = TRUE
  )

  expect_s3_class(p, "ggplot")
  expect_true(".panel" %in% names(p$data))
  expect_gt(length(unique(p$data$.panel)), 1)
})

test_that("plot_posterior_contours supports slice panels for 3 factors", {
  skip_if_not_installed("ggplot2")

  draws <- .create_workflow_draws_3d(n = 30, seed = 786)

  p <- plot_posterior_contours(
    draws = as_brsm_draws(draws, c("x1", "x2", "x3"),
      require_quadratic = FALSE,
      require_interactions = FALSE
    ),
    factor_names = c("x1", "x2", "x3"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2)),
    vary_factors = c("x1", "x2"),
    slice = list(x3 = c(-1, 1))
  )

  expect_s3_class(p, "ggplot")
  expect_equal(length(unique(p$data$.panel)), 2)
})

test_that("brsm_workflow forwards pairwise contour controls", {
  skip_if_not_installed("ggplot2")

  draws <- .create_workflow_draws_3d(n = 30, seed = 787)

  out <- suppressWarnings(brsm_workflow(
    object = draws,
    factor_names = c("x1", "x2", "x3"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2)),
    steps = c("plots"),
    include_plots = TRUE,
    plot_pairwise = TRUE
  ))

  expect_true("plots" %in% names(out))
  expect_s3_class(out$plots$contours, "ggplot")
  expect_true(".panel" %in% names(out$plots$contours$data))
  expect_gt(length(unique(out$plots$contours$data$.panel)), 1)
})

test_that("brsm_workflow projects optimum and ridge plots for selected pair", {
  skip_if_not_installed("ggplot2")

  draws <- .create_workflow_draws_3d_quadratic(n = 40, seed = 789)

  out <- suppressWarnings(brsm_workflow(
    object = draws,
    factor_names = c("x1", "x2", "x3"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2)),
    steps = c("plots"),
    include_plots = TRUE,
    plot_pairwise = FALSE,
    plot_vary_factors = c("x2", "x3")
  ))

  expect_true("plots" %in% names(out))
  expect_equal(out$plots$plot_mode, "optimization_geometry")
  expect_true(out$plots$geometry_available)
  expect_s3_class(out$plots$contours, "ggplot")
  expect_s3_class(out$plots$optimum, "ggplot")
  expect_s3_class(out$plots$ridge_path, "ggplot")
  expect_equal(out$plots$optimum$labels$x, "x2")
  expect_equal(out$plots$optimum$labels$y, "x3")
})

test_that("brsm_workflow creates pairwise projected optimum and ridge plots", {
  skip_if_not_installed("ggplot2")

  draws <- .create_workflow_draws_3d_quadratic(n = 40, seed = 790)

  out <- suppressWarnings(brsm_workflow(
    object = draws,
    factor_names = c("x1", "x2", "x3"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2)),
    steps = c("plots"),
    include_plots = TRUE,
    plot_pairwise = TRUE
  ))

  expect_true("plots" %in% names(out))
  expect_true("optimum_pairwise" %in% names(out$plots))
  expect_true("ridge_path_pairwise" %in% names(out$plots))
  expect_equal(length(out$plots$optimum_pairwise), 3)
  expect_equal(length(out$plots$ridge_path_pairwise), 3)

  expect_true(all(vapply(out$plots$optimum_pairwise, inherits, logical(1), "ggplot")))
  expect_true(all(vapply(out$plots$ridge_path_pairwise, inherits, logical(1), "ggplot")))

  expected_pairs <- c("x1_x2", "x1_x3", "x2_x3")
  expect_true(all(expected_pairs %in% names(out$plots$optimum_pairwise)))
  expect_true(all(expected_pairs %in% names(out$plots$ridge_path_pairwise)))
})

test_that("brsm_workflow validates plot_vary_factors", {
  draws <- .create_workflow_draws_3d_quadratic(n = 30, seed = 791)

  expect_error(
    brsm_workflow(
      object = draws,
      factor_names = c("x1", "x2", "x3"),
      ranges = list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2)),
      steps = c("plots"),
      include_plots = TRUE,
      plot_vary_factors = c("x1", "x4")
    ),
    "plot_vary_factors must be a subset of factor_names"
  )
})

test_that("brsm_workflow validates fit args, steps, and fit-mode inputs", {
  draws <- .create_workflow_draws(n = 20, seed = 792)

  expect_error(
    brsm_workflow(
      object = draws,
      factor_names = c("x1", "x2"),
      fit_args = 1,
      include_plots = FALSE
    ),
    "fit_args must be a named list"
  )

  expect_error(
    brsm_workflow(
      object = draws,
      factor_names = c("x1", "x2"),
      steps = character(0),
      include_plots = FALSE
    ),
    "steps must be a non-empty"
  )

  expect_error(
    brsm_workflow(
      object = draws,
      factor_names = c("x1", "x2"),
      steps = c("bad_step"),
      include_plots = FALSE
    ),
    "Unknown step"
  )

  expect_error(
    brsm_workflow(
      object = structure(list(a = 1), class = "not_df"),
      factor_names = c("x1", "x2"),
      fit_mode = "fit_brsm",
      response = "y",
      include_plots = FALSE
    ),
    "object must be a raw data.frame"
  )

  expect_error(
    brsm_workflow(
      object = generate_simulation_data(n = 8, seed = 793),
      factor_names = c("x1", "x2"),
      fit_mode = "fit_brsm",
      response = NULL,
      include_plots = FALSE
    ),
    "response must be supplied"
  )

  expect_error(
    brsm_workflow(
      object = draws,
      factor_names = c("x1", "x2"),
      grid_n = 1,
      include_plots = FALSE
    ),
    "grid_n must be a finite integer >= 2"
  )
})

test_that("brsm_workflow supports mocked fit_mode output and ranges fallback", {
  dat <- generate_simulation_data(n = 10, seed = 794)
  fit_draws <- .create_workflow_draws(n = 25, seed = 795)

  as.data.frame.fake_brmsfit_wf <<- function(x, ...) x$.draws_df
  registerS3method(
    "as.data.frame", "fake_brmsfit_wf",
    as.data.frame.fake_brmsfit_wf,
    envir = asNamespace("base")
  )

  fake_fit <- structure(
    list(
      fit = structure(list(.draws_df = fit_draws), class = c("fake_brmsfit_wf", "brmsfit")),
      factor_names = c("x1", "x2"),
      ranges = list(x1 = c(-1, 1), x2 = c(-1, 1)),
      model_terms = "second_order"
    ),
    class = "brsm_fit"
  )

  ns <- asNamespace("brsm")
  unlockBinding("fit_brsm", ns)
  old_fit_brsm <- get("fit_brsm", envir = ns)
  assign("fit_brsm", function(...) fake_fit, envir = ns)
  lockBinding("fit_brsm", ns)
  on.exit({
    unlockBinding("fit_brsm", ns)
    assign("fit_brsm", old_fit_brsm, envir = ns)
    lockBinding("fit_brsm", ns)
  }, add = TRUE)

  out <- brsm_workflow(
    object = dat,
    factor_names = c("x1", "x2"),
    fit_mode = "fit_brsm",
    response = "y",
    ranges = NULL,
    steps = c("predictions"),
    include_plots = FALSE,
    fit_args = list(chains = 1)
  )

  expect_true("fit" %in% names(out))
  expect_true("ranges" %in% names(out))
  expect_true("predictions" %in% names(out))
  expect_true("grid" %in% names(out))
})

test_that("brsm_workflow handles draw filtering, geometry outputs, and plot guards", {
  draws_some_na <- .create_workflow_draws(n = 20, seed = 796)
  draws_some_na$b_x1[1] <- NA_real_
  expect_warning(
    out_filtered <- brsm_workflow(
      object = draws_some_na,
      factor_names = c("x1", "x2"),
      ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
      steps = c("credible_region", "ridge"),
      include_plots = FALSE
    ),
    "Dropped"
  )
  expect_true("credible_region" %in% names(out_filtered))
  expect_true("ridge_analysis" %in% names(out_filtered))

  draws_all_na <- .create_workflow_draws(n = 8, seed = 797)
  draws_all_na$b_x1 <- NA_real_
  expect_error(
    suppressWarnings(brsm_workflow(
      object = draws_all_na,
      factor_names = c("x1", "x2"),
      ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
      steps = c("predictions"),
      include_plots = FALSE
    )),
    "No valid draw rows remain"
  )

  # Non-quadratic + precomputed steepest ascent triggers the out$steepest_ascent branch.
  set.seed(798)
  linear_draws <- data.frame(
    b_Intercept = rnorm(30, 10, 1),
    b_x1 = rnorm(30, 1, 0.2),
    b_x2 = rnorm(30, -0.5, 0.2),
    check.names = FALSE
  )
  out_dir <- suppressWarnings(brsm_workflow(
    object = linear_draws,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    steps = c("plots", "steepest_ascent"),
    include_plots = TRUE
  ))
  expect_true("direction" %in% names(out_dir$plots))

  # plots step ignored branch for single-factor workflows.
  one_factor <- data.frame(
    b_Intercept = 1:6,
    b_x1 = 2:7,
    check.names = FALSE
  )
  expect_warning(
    brsm_workflow(
      object = one_factor,
      factor_names = c("x1"),
      ranges = list(x1 = c(-1, 1)),
      steps = c("plots"),
      include_plots = TRUE
    ),
    "plots step ignored"
  )
})

test_that("brsm_workflow covers remaining internal branch edges", {
  skip_if_not_installed("ggplot2")

  draws2 <- .create_workflow_draws(n = 24, seed = 799)
  ranges2 <- list(x1 = c(-2, 2), x2 = c(-2, 2))

  # Line 146: include_plots=TRUE auto-adds plots even when omitted from steps.
  out_auto <- suppressWarnings(brsm_workflow(
    object = draws2,
    factor_names = c("x1", "x2"),
    ranges = ranges2,
    steps = c("stationary"),
    include_plots = TRUE
  ))
  expect_true("plots" %in% names(out_auto))

  # Line 507: non-quadratic branch with unavailable mean path.
  ns <- asNamespace("brsm")
  unlockBinding("steepest_ascent", ns)
  old_sa <- get("steepest_ascent", envir = ns)
  assign("steepest_ascent", function(...) list(mean_path = NULL), envir = ns)
  lockBinding("steepest_ascent", ns)

  on.exit({
    unlockBinding("steepest_ascent", ns)
    assign("steepest_ascent", old_sa, envir = ns)
    lockBinding("steepest_ascent", ns)
  }, add = TRUE)

  set.seed(800)
  linear_draws <- data.frame(
    b_Intercept = rnorm(20, 10, 1),
    b_x1 = rnorm(20, 1, 0.2),
    b_x2 = rnorm(20, -0.5, 0.2),
    check.names = FALSE
  )
  out_linear <- suppressWarnings(brsm_workflow(
    object = linear_draws,
    factor_names = c("x1", "x2"),
    ranges = ranges2,
    steps = c("plots"),
    include_plots = TRUE
  ))
  expect_true("plots" %in% names(out_linear))

  # Restore steepest_ascent before geometry-path mocks.
  unlockBinding("steepest_ascent", ns)
  assign("steepest_ascent", old_sa, envir = ns)
  lockBinding("steepest_ascent", ns)

  # Lines 385 and 531 via controlled mocks in quadratic 3-factor plotting path.
  draws3 <- .create_workflow_draws_3d_quadratic(n = 24, seed = 801)
  ranges3 <- list(x1 = c(-2, 2), x2 = c(-2, 2), x3 = c(-2, 2))

  unlockBinding("stationary_point", ns)
  unlockBinding("plot_posterior_contours", ns)
  unlockBinding("plot_optimum_posterior", ns)
  unlockBinding("posterior_ridge_analysis", ns)
  unlockBinding("plot_ridge_path", ns)

  old_sp <- get("stationary_point", envir = ns)
  old_contours <- get("plot_posterior_contours", envir = ns)
  old_opt <- get("plot_optimum_posterior", envir = ns)
  old_ridge <- get("posterior_ridge_analysis", envir = ns)
  old_ridge_plot <- get("plot_ridge_path", envir = ns)

  assign("stationary_point", function(...) NULL, envir = ns)
  assign("plot_posterior_contours", function(...) ggplot2::ggplot(), envir = ns)
  assign("plot_optimum_posterior", function(...) ggplot2::ggplot(), envir = ns)
  assign("posterior_ridge_analysis", function(object, factor_names, radii, summary = FALSE, ...) {
    data.frame(
      x1 = NA_real_,
      x2 = NA_real_,
      x3 = NA_real_,
      radius = NA_real_
    )
  }, envir = ns)
  assign("plot_ridge_path", function(...) ggplot2::ggplot(), envir = ns)

  lockBinding("stationary_point", ns)
  lockBinding("plot_posterior_contours", ns)
  lockBinding("plot_optimum_posterior", ns)
  lockBinding("posterior_ridge_analysis", ns)
  lockBinding("plot_ridge_path", ns)

  on.exit({
    unlockBinding("stationary_point", ns)
    unlockBinding("plot_posterior_contours", ns)
    unlockBinding("plot_optimum_posterior", ns)
    unlockBinding("posterior_ridge_analysis", ns)
    unlockBinding("plot_ridge_path", ns)
    assign("stationary_point", old_sp, envir = ns)
    assign("plot_posterior_contours", old_contours, envir = ns)
    assign("plot_optimum_posterior", old_opt, envir = ns)
    assign("posterior_ridge_analysis", old_ridge, envir = ns)
    assign("plot_ridge_path", old_ridge_plot, envir = ns)
    lockBinding("stationary_point", ns)
    lockBinding("plot_posterior_contours", ns)
    lockBinding("plot_optimum_posterior", ns)
    lockBinding("posterior_ridge_analysis", ns)
    lockBinding("plot_ridge_path", ns)
  }, add = TRUE)

  out_mocked <- suppressWarnings(brsm_workflow(
    object = draws3,
    factor_names = c("x1", "x2", "x3"),
    ranges = ranges3,
    steps = c("plots"),
    include_plots = TRUE,
    plot_pairwise = TRUE
  ))
  expect_true("plots" %in% names(out_mocked))
})