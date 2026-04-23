# Integration tests: full workflow end-to-end

.integration_cache <- new.env(parent = emptyenv())

.get_integration_data <- function() {
  if (!exists("data", envir = .integration_cache, inherits = FALSE)) {
    .integration_cache$data <- generate_simulation_data(n = 30, seed = 601)
  }
  .integration_cache$data
}

.fit_integration_model <- function(
  congruence_type = NULL,
  coded = FALSE,
  seed = 601) {
  dat <- .get_integration_data()
  input_data <- if (coded) {
    prepare_brsm_data(dat, factor_names = c("x1", "x2"), method = "zscore")
  } else {
    dat
  }

  if (is.null(congruence_type)) {
    fit_brsm(
      data = input_data,
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
  } else {
    fit_brsm_congruence(
      data = input_data,
      response = "y",
      factor_names = c("x1", "x2"),
      congruence_type = congruence_type,
      chains = 1,
      iter = 250,
      warmup = 125,
      seed = seed,
      sampling_preset = "fast",
      refresh = 0,
      silent = 2
    )
  }
}

test_that("Full workflow: prepare -> fit -> congruence parameters -> rope", {
  skip_if_no_brms_tests()

  prepared <- prepare_brsm_data(
    .get_integration_data(),
    factor_names = c("x1", "x2"),
    method = "zscore"
  )
  fit <- .fit_integration_model(coded = TRUE, seed = 602)
  posterior <- as_brsm_draws(fit)
  params <- congruence_parameters(posterior, factor_names = c("x1", "x2"))
  summary_df <- summarize_congruence(posterior, factor_names = c("x1", "x2"))
  rope <- rope_congruence(
    posterior,
    factor_names = c("x1", "x2"),
    rope = c(-0.1, 0.1)
  )

  expect_s3_class(prepared, "brsm_coded_data")
  expect_s3_class(fit, "brsm_fit")
  expect_data_frame(params, ncols = 5)
  expect_data_frame(summary_df, nrows = 5)
  expect_data_frame(rope, nrows = 5)
})

test_that("Workflow with constrained congruence model", {
  skip_if_no_brms_tests()

  fit_broad <- .fit_integration_model(congruence_type = "broad", seed = 603)
  fit_strict <- .fit_integration_model(congruence_type = "strict", seed = 604)

  expect_s3_class(fit_broad, "brsm_fit")
  expect_s3_class(fit_strict, "brsm_fit")
  expect_equal(fit_broad$congruence_type, "broad")
  expect_equal(fit_strict$congruence_type, "strict")
})

test_that("Workflow: fit -> LOF test -> comparison", {
  skip_if_no_brms_tests()

  baseline <- .fit_integration_model(seed = 605)
  result <- loftest_brsm(
    object = baseline,
    reference_type = "cubic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 606,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_type(result, "list")
  expect_true("comparison" %in% names(result))
  expect_true(is.data.frame(result$comparison$comparison))
})

test_that("Workflow: decode predictions to original scale", {
  skip_if_no_brms_tests()

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
  skip_if_no_brms_tests()

  fit <- .fit_integration_model(seed = 608)

  expect_no_error(print(fit))
  expect_no_error(summary(fit))
  print_output <- capture.output(print(fit))
  summary_output <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Model Formula", print_output, fixed = TRUE)))
  expect_true(any(grepl("Coefficient Summary", summary_output, fixed = TRUE)))
})

test_that(
  "Backward compatibility: existing fit_brsm models work with new functions",
  {
  skip_if_no_brms_tests()

  fit <- .fit_integration_model(seed = 614)
  old_style_fit <- fit
  old_style_fit$coding <- NULL
  old_style_fit$congruence_type <- NULL

  draws <- as_brsm_draws(old_style_fit)
  params <- congruence_parameters(draws, factor_names = c("x1", "x2"))
  summary_df <- summarize_congruence(draws, factor_names = c("x1", "x2"))

  expect_data_frame(params, ncols = 5)
  expect_data_frame(summary_df, nrows = 5)
  expect_no_error(print(old_style_fit))
}
)

test_that("Workflow: multiple congruence types compared via LOF", {
  skip_if_no_brms_tests()

  fit_unconstrained <- .fit_integration_model(
    congruence_type = "unconstrained",
    seed = 609
  )
  fit_broad <- .fit_integration_model(congruence_type = "broad", seed = 610)
  fit_strict <- .fit_integration_model(congruence_type = "strict", seed = 611)

  comp <- compare_brsm_models(
    models = list(
      unconstrained = fit_unconstrained,
      broad = fit_broad,
      strict = fit_strict
    ),
    criterion = "loo"
  )

  params_uncon <- congruence_parameters(
    as_brsm_draws(fit_unconstrained),
    factor_names = c("x1", "x2")
  )
  rope_uncon <- rope_congruence(
    as_brsm_draws(fit_unconstrained),
    factor_names = c("x1", "x2"),
    rope = c(-0.1, 0.1)
  )

  expect_true(is.data.frame(comp$comparison))
  expect_data_frame(params_uncon, ncols = 5)
  expect_data_frame(rope_uncon, nrows = 5)
  expect_equal(
    sort(unique(c(fit_broad$congruence_type, fit_strict$congruence_type))),
    c("broad", "strict")
  )
})

test_that(
  "Integration: prepare -> fit_congruence -> summarize",
  {
  skip_if_no_brms_tests()

  prepared <- prepare_brsm_data(
    .get_integration_data(),
    c("x1", "x2"),
    method = "zscore"
  )
  fit <- .fit_integration_model(
    congruence_type = "unconstrained",
    coded = TRUE,
    seed = 612
  )
  posterior <- as_brsm_draws(fit)
  params <- congruence_parameters(posterior, factor_names = c("x1", "x2"))
  summary_df <- summarize_congruence(posterior, factor_names = c("x1", "x2"))

  expect_s3_class(prepared, "brsm_coded_data")
  expect_s3_class(fit, "brsm_fit")
  expect_data_frame(params, ncols = 5)
  expect_data_frame(summary_df, nrows = 5)
}
)

test_that("Coding metadata persists through full workflow", {
  skip_if_no_brms_tests()

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
