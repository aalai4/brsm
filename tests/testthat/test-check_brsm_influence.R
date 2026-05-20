# Tests for check_brsm_influence internal helpers and public function

test_that(".brsm_residual_scale returns sd for regular residuals", {
  x <- c(-2, -1, 0, 1, 2)
  s <- brsm:::.brsm_residual_scale(x)
  expect_true(is.numeric(s))
  expect_true(is.finite(s))
  expect_true(s > 0)
  expect_equal(s, stats::sd(x), tolerance = 1e-10)
})

test_that(".brsm_residual_scale falls back when sd is zero", {
  x <- rep(5, 8)
  s <- brsm:::.brsm_residual_scale(x)
  expect_true(is.numeric(s))
  expect_true(is.finite(s))
  expect_true(s > 0)
})

test_that(".brsm_extract_pareto_k returns empty vector for invalid n_obs", {
  out <- brsm:::.brsm_extract_pareto_k(fit = NULL, n_obs = 0)
  expect_equal(length(out), 0L)
})

test_that(".brsm_extract_pareto_k returns NA vector when fit cannot be processed", {
  skip_if_not_installed("loo")
  k <- brsm:::.brsm_extract_pareto_k(fit = list(), n_obs = 5)
  expect_equal(length(k), 5L)
  expect_true(all(is.na(k)))
})

test_that("check_brsm_influence validates scalar arguments", {
  if (!requireNamespace("brms", quietly = TRUE)) {
    skip("brms not installed")
  }

  fake <- structure(list(), class = "brmsfit")

  expect_error(
    check_brsm_influence(fake, ndraws = 0),
    "ndraws must be a positive finite integer"
  )
  expect_error(
    check_brsm_influence(fake, outlier_sd_threshold = 0),
    "outlier_sd_threshold must be a positive finite numeric scalar"
  )
  expect_error(
    check_brsm_influence(fake, pareto_k_threshold = 0),
    "pareto_k_threshold must be a positive finite numeric scalar"
  )
})

test_that("check_brsm_influence returns expected structure on real fit", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 30, seed = 777)
  fit <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 777,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  out <- suppressMessages(check_brsm_influence(
    fit,
    ndraws = 100,
    verbose = TRUE,
    include_plot = FALSE,
    seed = 777
  ))

  expect_s3_class(out$overview, "data.frame")
  expect_s3_class(out$observations, "data.frame")
  expect_equal(nrow(out$overview), 1L)
  expect_equal(nrow(out$observations), nrow(dat))

  expect_true(all(c(
    "obs_id", "observed", "predicted_mean", "residual", "std_residual",
    "outlier", "pareto_k", "influential"
  ) %in% names(out$observations)))

  expect_true(is.logical(out$passed))
})

test_that("check_brsm_influence can return plot", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 20, seed = 778)
  fit <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 778,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  out <- suppressMessages(check_brsm_influence(
    fit,
    ndraws = 80,
    include_plot = TRUE,
    verbose = FALSE,
    seed = 778
  ))

  expect_true("plot" %in% names(out))
  expect_s3_class(out$plot, "ggplot")
})