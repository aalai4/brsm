# Tests for brsm_fit S3 dispatch in downstream analysis

test_that("predict_surface dispatches on brsm_fit", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 780)

  fit_obj <- suppressWarnings(fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 780,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2,
    coding_policy = "ignore"
  ))

  grid <- surface_grid(ranges = fit_obj$ranges, n = 5)

  out <- predict_surface(
    draws = fit_obj,
    newdata = grid,
    summary = TRUE
  )

  expect_true(is.data.frame(out))
  expect_true(all(c("x1", "x2", "mean") %in% names(out)))
})

test_that("brsm_workflow uses dispatch object cleanly for brsm_fit", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 781)

  fit_obj <- suppressWarnings(fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 781,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2,
    coding_policy = "ignore"
  ))

  out <- brsm_workflow(
    object = fit_obj,
    factor_names = c("x1", "x2"),
    ranges = fit_obj$ranges,
    include_plots = FALSE,
    steps = c(
      "predictions", "stationary", "classification",
      "credible_region", "steepest_ascent"
    )
  )

  expect_true("predictions" %in% names(out))
  expect_true("stationary_points" %in% names(out))
  expect_true("classification" %in% names(out))
  expect_true("credible_region" %in% names(out))
  expect_true("steepest_ascent" %in% names(out))
})
