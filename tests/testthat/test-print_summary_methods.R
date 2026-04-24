# Tests for print.brsm_fit, summary.brsm_fit, and print.summary.brsm_fit

.psm_cache <- new.env(parent = emptyenv())

.psm_fit <- function(key = "plain", seed = 401) {
  cache_key <- paste0(key, "_", seed)
  if (!exists(cache_key, envir = .psm_cache, inherits = FALSE)) {
    dat <- generate_simulation_data(n = 25, seed = seed)
    fit <- switch(key,
      plain = fit_brsm(
        data = dat,
        response = "y",
        factor_names = c("x1", "x2"),
        chains = 1,
        iter = 250,
        warmup = 125,
        seed = seed,
        sampling_preset = "fast",
        refresh = 0,
        silent = 2
      ),
      coded_z = {
        prepared <- prepare_brsm_data(
          dat,
          factor_names = c("x1", "x2"),
          method = "zscore"
        )
        fit_brsm(
          data = prepared,
          response = "y",
          factor_names = c("x1", "x2"),
          chains = 1,
          iter = 250,
          warmup = 125,
          seed = seed,
          sampling_preset = "fast",
          refresh = 0,
          silent = 2
        )
      },
      coded_range = {
        prepared <- prepare_brsm_data(
          dat,
          factor_names = c("x1", "x2"),
          method = "range"
        )
        fit_brsm(
          data = prepared,
          response = "y",
          factor_names = c("x1", "x2"),
          chains = 1,
          iter = 250,
          warmup = 125,
          seed = seed,
          sampling_preset = "fast",
          refresh = 0,
          silent = 2
        )
      },
      stop("Unknown key: ", key)
    )
    assign(cache_key, fit, envir = .psm_cache)
  }
  get(cache_key, envir = .psm_cache, inherits = FALSE)
}

test_that("print.brsm_fit displays formula", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 401)
  fit <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 401,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result <- capture.output(print(fit))
  expect_true(any(grepl("Model Formula", result, fixed = TRUE)))
  expect_true(any(grepl("x1", result, fixed = TRUE)))
  expect_true(any(grepl("x2", result, fixed = TRUE)))
})

test_that("print.brsm_fit displays response variable name", {
  skip_if_no_brms_tests()

  result <- capture.output(print(.psm_fit("plain", seed = 406)))
  expect_true(any(grepl("Response Variable", result, fixed = TRUE)))
  expect_true(any(grepl("y", result, fixed = TRUE)))
})

test_that("print.brsm_fit displays factor names", {
  skip_if_no_brms_tests()

  result <- capture.output(print(.psm_fit("plain", seed = 407)))
  expect_true(any(grepl("Factor Variables", result, fixed = TRUE)))
  expect_true(any(grepl("x1", result, fixed = TRUE)))
  expect_true(any(grepl("x2", result, fixed = TRUE)))
})

test_that("print.brsm_fit displays factor ranges", {
  skip_if_no_brms_tests()

  result <- capture.output(print(.psm_fit("plain", seed = 408)))
  expect_true(any(grepl("Factor Ranges", result, fixed = TRUE)))
  expect_true(any(grepl("x1", result, fixed = TRUE)))
})

test_that("print.brsm_fit displays coding method if applicable", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 403)
  prepared <- prepare_brsm_data(
    dat,
    factor_names = c("x1", "x2"),
    method = "zscore"
  )
  fit <- fit_brsm(
    data = prepared,
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 403,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result <- capture.output(print(fit))
  expect_true(any(grepl("Coding Method", result, fixed = TRUE)))
  expect_true(any(grepl("zscore", result, fixed = TRUE)))
})

test_that("print.brsm_fit displays MCMC configuration", {
  skip_if_no_brms_tests()

  result <- capture.output(print(.psm_fit("plain", seed = 409)))
  expect_true(any(grepl("MCMC Sampling", result, fixed = TRUE)))
  expect_true(any(grepl("Chains", result, fixed = TRUE)))
})

test_that("print.brsm_fit displays convergence diagnostics (Rhat)", {
  skip_if_no_brms_tests()

  result <- capture.output(print(.psm_fit("plain", seed = 410)))
  expect_true(any(grepl("Rhat", result, fixed = TRUE)))
})

test_that("print.brsm_fit produces readable output without errors", {
  skip_if_no_brms_tests()

  fit <- .psm_fit("plain", seed = 411)
  output <- capture.output(expect_no_error(print(fit)))
  expect_true(length(output) > 0)
})

test_that("summary.brsm_fit delegates to brms::summary", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 404)
  fit <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 404,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result <- summary(fit)
  expect_s3_class(result, "summary.brsm_fit")
  expect_true("brmsfit_summary" %in% names(result))
  expect_equal(result$response, "y")
})

test_that("summary.brsm_fit includes brsm metadata in output", {
  skip_if_no_brms_tests()

  fit <- .psm_fit("plain", seed = 412)
  result <- summary(fit)

  expect_identical(result$response, fit$response)
  expect_identical(result$factor_names, fit$factor_names)
  expect_true(!is.null(result$formula))
})

test_that("print.summary.brsm_fit produces readable output", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 405)
  fit <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 405,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )
  summary_obj <- summary(fit)

  result <- capture.output(print(summary_obj))
  expect_true(any(
    grepl("Bayesian Response Surface Model Summary", result, fixed = TRUE)
  ))
  expect_true(any(grepl("Coefficient Summary", result, fixed = TRUE)))
})

test_that("print.summary.brsm_fit includes coefficient table", {
  skip_if_no_brms_tests()

  summary_obj <- summary(.psm_fit("plain", seed = 414))
  result <- capture.output(print(summary_obj))
  expect_true(any(grepl("Coefficient Summary", result, fixed = TRUE)))
})

test_that("S3 method dispatch works correctly", {
  skip_if_no_brms_tests()

  fit <- .psm_fit("plain", seed = 415)
  expect_s3_class(fit, "brsm_fit")
  expect_s3_class(summary(fit), "summary.brsm_fit")
  expect_no_error(print(fit))
})

test_that("print output with different coding methods differs appropriately", {
  skip_if_no_brms_tests()

  output_z <- capture.output(print(.psm_fit("coded_z", seed = 419)))
  output_r <- capture.output(print(.psm_fit("coded_range", seed = 420)))

  expect_true(any(grepl("zscore", output_z, fixed = TRUE)))
  expect_true(any(grepl("range", output_r, fixed = TRUE)))
  expect_false(identical(output_z, output_r))
})

test_that("print handles missing metadata gracefully", {
  skip_if_no_brms_tests()

  fit <- .psm_fit("plain", seed = 421)
  fit_missing <- fit
  fit_missing$ranges <- NULL
  fit_missing$coding <- NULL

  expect_no_error(print(fit_missing))
})

test_that("summary.brsm_fit with different chains/iterations", {
  skip_if_no_brms_tests()

  fit_small <- fit_brsm(
    data = generate_simulation_data(n = 25, seed = 422),
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 200,
    warmup = 100,
    seed = 422,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  fit_large <- fit_brsm(
    data = generate_simulation_data(n = 25, seed = 423),
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 300,
    warmup = 150,
    seed = 423,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result_small <- summary(fit_small)
  result_large <- summary(fit_large)

  expect_s3_class(result_small, "summary.brsm_fit")
  expect_s3_class(result_large, "summary.brsm_fit")
  expect_equal(fit_small$sampling$iter, 200)
  expect_equal(fit_large$sampling$iter, 300)
})
