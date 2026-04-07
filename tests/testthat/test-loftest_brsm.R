# Tests for loftest_brsm

.loftest_cache <- new.env(parent = emptyenv())

.get_loftest_baseline <- function(seed = 301) {
  key <- paste0("baseline_", seed)
  if (!exists(key, envir = .loftest_cache, inherits = FALSE)) {
    dat <- prepare_brsm_data(
      generate_simulation_data(n = 25, seed = seed),
      factor_names = c("x1", "x2"),
      method = "zscore"
    )
    fit <- fit_brsm(
      data = dat,
      response = "y",
      factor_names = c("x1", "x2"),
      coding_policy = "ignore",
      chains = 1,
      iter = 250,
      warmup = 125,
      seed = seed,
      sampling_preset = "fast",
      refresh = 0,
      silent = 2
    )
    assign(key, fit, envir = .loftest_cache)
  }
  get(key, envir = .loftest_cache, inherits = FALSE)
}

test_that("loftest_brsm with provided reference model", {
  skip_if_no_brms_tests()

  baseline <- .get_loftest_baseline(seed = 301)
  reference <- fit_brsm(
    data = baseline$fit$data,
    response = "y",
    factor_names = c("x1", "x2"),
    coding_policy = "ignore",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 302,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result <- loftest_brsm(
    object = baseline,
    reference_model = reference,
    include_ppc = FALSE,
    criterion = "loo"
  )

  expect_type(result, "list")
  expect_true("comparison" %in% names(result))
  expect_false(isTRUE(result$reference_fitted))
})

test_that("loftest_brsm with auto-fit cubic reference model", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 303)
  baseline <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    coding_policy = "ignore",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 303,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result <- loftest_brsm(
    object = baseline,
    reference_type = "cubic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 303,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_type(result, "list")
  expect_true(isTRUE(result$reference_fitted))
  expect_true("comparison" %in% names(result))
  expect_true("comparison" %in% names(result$comparison))
})

test_that("loftest_brsm with auto-fit extended reference model", {
  skip_if_no_brms_tests()

  baseline <- .get_loftest_baseline(seed = 304)
  result <- loftest_brsm(
    object = baseline,
    reference_type = "extended",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 304,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  ref_formula <- paste(deparse(stats::formula(result$reference_model)), collapse = " ")
  expect_type(result, "list")
  expect_true(isTRUE(result$reference_fitted))
  expect_true(grepl("I\\(x1\\^3\\)", ref_formula))
  expect_true(grepl("I\\(x2\\^3\\)", ref_formula))
})

test_that("loftest_brsm reference_type='cubic' adds cubic terms", {
  f <- brsm:::.brsm_build_reference_formula(
    response = "y",
    factor_names = c("x1", "x2"),
    reference_type = "cubic"
  )
  txt <- paste(deparse(f), collapse = " ")

  expect_true(grepl("I\\(x1\\^3\\)", txt))
  expect_true(grepl("I\\(x2\\^3\\)", txt))
  expect_true(grepl("x1:x2", txt))
})

test_that("loftest_brsm reference_type='extended' adds cubic and interaction terms", {
  f <- brsm:::.brsm_build_reference_formula(
    response = "y",
    factor_names = c("x1", "x2"),
    reference_type = "extended"
  )
  txt <- paste(deparse(f), collapse = " ")

  expect_true(grepl("I\\(x1\\^3\\)", txt))
  expect_true(grepl("I\\(x2\\^3\\)", txt))
  expect_true(grepl("I\\(x1\\^2\\):x2", txt))
  expect_true(grepl("x1:I\\(x2\\^2\\)", txt))
})

test_that("loftest_brsm without include_ppc", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 305),
    reference_type = "cubic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 305,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_false("ppc" %in% names(result))
})

test_that("loftest_brsm with include_ppc=TRUE", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 306),
    reference_type = "cubic",
    include_ppc = TRUE,
    ppc_ndraws = 50,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 306,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_true("ppc" %in% names(result))
  expect_true("summaries" %in% names(result$ppc))
})

test_that("loftest_brsm returns comparison with loo_diff", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 307),
    reference_type = "cubic",
    include_ppc = FALSE,
    criterion = "loo",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 307,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_true("comparison" %in% names(result$comparison))
  expect_true("elpd_diff" %in% names(result$comparison$comparison))
})

test_that("loftest_brsm richer model has lower LOO (better fit)", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 308),
    reference_type = "extended",
    include_ppc = FALSE,
    criterion = "loo",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 308,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_true(is.data.frame(result$comparison$comparison))
  expect_true(nrow(result$comparison$comparison) >= 2)
})

test_that("loftest_brsm errors when reference model is simpler than baseline", {
  skip_if_no_brms_tests()

  baseline <- .get_loftest_baseline(seed = 309)
  simpler_reference <- brms::brm(
    formula = y ~ x1 + x2,
    data = generate_simulation_data(n = 25, seed = 310),
    family = stats::gaussian(),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 310,
    control = list(adapt_delta = 0.8, max_treedepth = 10),
    refresh = 0,
    silent = 2
  )

  expect_no_error(
    loftest_brsm(object = baseline, reference_model = simpler_reference, include_ppc = FALSE)
  )
})

test_that("loftest_brsm handles single factor case", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 311)
  dat$x2 <- NULL
  baseline <- brms::brm(
    formula = y ~ x1 + I(x1^2),
    data = dat,
    family = stats::gaussian(),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 311,
    control = list(adapt_delta = 0.8, max_treedepth = 10),
    refresh = 0,
    silent = 2
  )

  result <- loftest_brsm(
    object = baseline,
    data = dat,
    response = "y",
    factor_names = "x1",
    reference_type = "cubic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 311,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_type(result, "list")
  expect_true(grepl("I\\(x1\\^3\\)", paste(deparse(stats::formula(result$reference_model)), collapse = " ")))
})

test_that("loftest_brsm handles two factor case", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 312),
    reference_type = "cubic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 312,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  ref_formula <- paste(deparse(stats::formula(result$reference_model)), collapse = " ")
  expect_true(grepl("I\\(x1\\^3\\)", ref_formula))
  expect_true(grepl("I\\(x2\\^3\\)", ref_formula))
})

test_that("loftest_brsm reference model inherits priors from baseline", {
  skip_if_no_brms_tests()

  pr <- brms::prior(normal(0, 1), class = "b")
  dat <- generate_simulation_data(n = 25, seed = 313)
  baseline <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    prior = pr,
    coding_policy = "ignore",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 313,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result <- loftest_brsm(
    object = baseline,
    reference_type = "cubic",
    prior = pr,
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 313,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_type(result, "list")
  expect_true(isTRUE(result$reference_fitted))
})

test_that("loftest_brsm comparison returns LOO and/or WAIC", {
  skip_if_no_brms_tests()

  result_loo <- loftest_brsm(
    object = .get_loftest_baseline(seed = 314),
    reference_type = "cubic",
    criterion = "loo",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 314,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result_waic <- loftest_brsm(
    object = .get_loftest_baseline(seed = 315),
    reference_type = "cubic",
    criterion = "waic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 315,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_equal(result_loo$criterion, "loo")
  expect_equal(result_waic$criterion, "waic")
  expect_true(is.data.frame(result_loo$comparison$comparison))
  expect_true(is.data.frame(result_waic$comparison$comparison))
})
