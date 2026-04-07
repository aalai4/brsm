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
