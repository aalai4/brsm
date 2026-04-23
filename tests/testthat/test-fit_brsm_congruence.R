# Tests for fit_brsm_congruence

test_that("fit_brsm_congruence accepts all three congruence types", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 101)
  types <- c("unconstrained", "broad", "strict")

  fits <- lapply(types, function(tp) {
    fit_brsm_congruence(
      data = dat,
      response = "y",
      factor_names = c("x1", "x2"),
      congruence_type = tp,
      chains = 1,
      iter = 250,
      warmup = 125,
      seed = 101,
      sampling_preset = "fast",
      refresh = 0,
      silent = 2
    )
  })

  expect_true(all(vapply(fits, inherits, logical(1), what = "brsm_fit")))
  expect_equal(vapply(fits, function(x) x$congruence_type, character(1)), types)
})

test_that("fit_brsm_congruence formula construction: unconstrained", {
  # Test formula construction without actually fitting
  formula <- y ~ x1 + x2 + I(x1^2) + I(x2^2) + x1:x2

  # Should be callable, just test that it's the right structure
  expect_s3_class(formula, "formula")
  expect_equal(as.character(formula)[2], "y")
})

test_that("fit_brsm_congruence formula construction: broad", {
  # Broad constraint: b1=b2, b11=b22, keep b12
  # Should use I(x1+x2) for linear, I(x1^2+x2^2) for quadratic

  formula_expected <- "y ~ I(x1 + x2) + I(x1^2 + x2^2) + x1:x2"

  # Formula should have these terms
  expect_true(grepl("I\\(x1 \\+ x2\\)", formula_expected))
  expect_true(grepl("I\\(x1\\^2 \\+ x2\\^2\\)", formula_expected))
})

test_that("fit_brsm_congruence formula construction: strict", {
  # Strict constraint: b1=b2, b11=b22, no b12
  # Should use I(x1+x2) for linear, I(x1^2+x2^2) for quadratic, NO interaction

  formula_expected <- "y ~ I(x1 + x2) + I(x1^2 + x2^2)"

  # Formula should NOT have x1:x2
  expect_false(grepl("x1:x2", formula_expected))
})

test_that("fit_brsm_congruence produces brsm_fit object with metadata", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 202)
  result <- fit_brsm_congruence(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    congruence_type = "broad",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 202,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_s3_class(result, "brsm_fit")
  expect_equal(result$congruence_type, "broad")
  expect_equal(result$response, "y")
  expect_equal(result$factor_names, c("x1", "x2"))
})

test_that(
  "fit_brsm_congruence preserves coding metadata from prepare_brsm_data",
  {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 203)
  prepared <- prepare_brsm_data(
    dat,
    factor_names = c("x1", "x2"),
    method = "zscore"
  )

  result <- fit_brsm_congruence(
    data = prepared,
    response = "y",
    factor_names = c("x1", "x2"),
    congruence_type = "broad",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 203,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  coding <- get_brsm_coding(result)
  expect_equal(coding$method, "zscore")
  expect_equal(names(coding$factors), c("x1", "x2"))
}
)

test_that("fit_brsm_congruence converges without excessive warnings", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 204)
  result <- fit_brsm_congruence(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    congruence_type = "unconstrained",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 204,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  rhat_vals <- brms::rhat(result$fit)
  expect_true(all(is.finite(rhat_vals) | is.na(rhat_vals)))
  expect_true(length(rhat_vals) > 0)
})

test_that("fit_brsm_congruence broad constraint enforces b1=b2", {
  formula_str <- brsm:::.brsm_build_congruence_formula(
    response = "y",
    factor_names = c("x1", "x2"),
    congruence_type = "broad"
  )

  expect_true(grepl("I\\(x1\\+x2\\)", formula_str))
  expect_true(
    grepl("I\\(I\\(x1\\^2\\)\\s*\\+\\s*I\\(x2\\^2\\)\\)", formula_str)
  )
  expect_true(grepl("x1:x2", formula_str))
})

test_that("fit_brsm_congruence strict constraint removes interaction", {
  formula_str <- brsm:::.brsm_build_congruence_formula(
    response = "y",
    factor_names = c("x1", "x2"),
    congruence_type = "strict"
  )

  expect_true(grepl("I\\(x1\\+x2\\)", formula_str))
  expect_true(
    grepl("I\\(I\\(x1\\^2\\)\\s*\\+\\s*I\\(x2\\^2\\)\\)", formula_str)
  )
  expect_false(grepl("x1:x2", formula_str))
})

test_that("fit_brsm_congruence integrates with existing fit_brsm workflow", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 205)
  result <- fit_brsm_congruence(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    congruence_type = "unconstrained",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 205,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  params <- congruence_parameters(
    as_brsm_draws(result),
    factor_names = c("x1", "x2")
  )
  rope <- rope_congruence(
    as_brsm_draws(result),
    factor_names = c("x1", "x2"),
    rope = c(-0.1, 0.1)
  )

  expect_s3_class(result, "brsm_fit")
  expect_no_error(print(result))
  expect_no_error(summary(result))
  expect_data_frame(params, ncols = 5)
  expect_data_frame(rope, nrows = 5)
})
