# Integration tests: full workflow end-to-end

.integration_cache <- new.env(parent = emptyenv())

.get_integration_data <- function() {
  if (!exists("data", envir = .integration_cache, inherits = FALSE)) {
    .integration_cache$data <- generate_simulation_data(n = 30, seed = 601)
  }
  .integration_cache$data
}

.fit_integration_model <- function(
  coded = FALSE,
  seed = 601) {
  dat <- .get_integration_data()
  input_data <- if (coded) {
    prepare_brsm_data(dat, factor_names = c("x1", "x2"), method = "zscore")
  } else {
    dat
  }

  fit_brsm(
    data = input_data,
    response = "y",
    factor_names = c("x1", "x2"),
    draws = 250,
    seed = seed,
    coding_policy = "ignore"
  )
}

test_that("Full workflow: prepare -> fit -> draw conversion", {
  prepared <- prepare_brsm_data(
    .get_integration_data(),
    factor_names = c("x1", "x2"),
    method = "zscore"
  )
  fit <- .fit_integration_model(coded = TRUE, seed = 602)
  posterior <- as_brsm_draws(fit)

  expect_s3_class(prepared, "brsm_coded_data")
  expect_s3_class(fit, "brsm_fit")
  expect_true(is.data.frame(posterior))
  expect_gt(nrow(posterior), 0)
})

test_that("Workflow: decode predictions to original scale", {
  dat <- .get_integration_data()
  prepared <- prepare_brsm_data(
    dat,
    factor_names = c("x1", "x2"),
    method = "zscore"
  )
  fit <- .fit_integration_model(coded = TRUE, seed = 607)
  decoded <- decode_brsm_data(
    prepared[, c("x1", "x2")],
    coding = get_brsm_coding(fit)
  )

  expect_true(all(abs(decoded$x1 - dat$x1) < 1e-8))
  expect_true(all(abs(decoded$x2 - dat$x2) < 1e-8))
})

test_that("Print and summary workflow", {
  fit <- .fit_integration_model(seed = 608)

  expect_no_error(print(fit))
  expect_no_error(summary(fit))
  print_output <- capture.output(print(fit))
  summary_output <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Model Formula", print_output, fixed = TRUE)))
  expect_true(any(grepl("Coefficient Summary", summary_output, fixed = TRUE)))
})

test_that("Coding metadata persists through full workflow", {
  prepared <- prepare_brsm_data(
    .get_integration_data(),
    c("x1", "x2"),
    method = "zscore"
  )
  coding_before <- get_brsm_coding(prepared)
  fit <- .fit_integration_model(coded = TRUE, seed = 613)
  coding_after <- get_brsm_coding(fit)

  expect_equal(coding_before$method, coding_after$method)
  expect_equal(names(coding_before$factors), names(coding_after$factors))
  expect_equal(coding_before$factors$x1$center, coding_after$factors$x1$center)
  expect_equal(coding_before$factors$x2$scale, coding_after$factors$x2$scale)
})
