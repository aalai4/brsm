# Tests for print.brsm_fit using a mock brsm_fit object (no brms needed)

.make_mock_brsm_fit <- function(include_fit = FALSE,
                                 include_ranges = TRUE,
                                 include_sampling = TRUE,
                                 include_coding = TRUE) {
  obj <- list(
    formula      = y ~ x1 + x2 + I(x1^2) + I(x2^2) + x1:x2,
    response     = "y",
    factor_names = c("x1", "x2"),
    coding       = if (include_coding) list(method = "coded") else NULL,
    ranges       = if (include_ranges)
                     list(x1 = c(-1, 1), x2 = c(-2, 2)) else NULL,
    sampling     = if (include_sampling)
                     list(chains = 4, iter = 2000, warmup = 1000,
                          sampling_preset = "standard") else NULL,
    fit          = NULL,   # NULL avoids brms dependency
    call         = quote(fit_brsm(dat, response = "y",
                                  factor_names = c("x1", "x2")))
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

nuts_params.fake_brmsfit <- function(object, ...) {
  object$.nuts_params
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

if (requireNamespace("brms", quietly = TRUE)) {
  registerS3method(
    "nuts_params",
    "fake_brmsfit",
    nuts_params.fake_brmsfit,
    envir = asNamespace("brms")
  )
}

# ── print.brsm_fit ────────────────────────────────────────────────────────────

test_that("print.brsm_fit outputs header and formula", {
  mock <- .make_mock_brsm_fit()
  out  <- capture.output(print(mock))
  expect_true(any(grepl("Bayesian Response Surface Model", out)))
  expect_true(any(grepl("y ~", out)))
})

test_that("print.brsm_fit shows factor names", {
  mock <- .make_mock_brsm_fit()
  out  <- capture.output(print(mock))
  expect_true(any(grepl("x1", out)))
  expect_true(any(grepl("x2", out)))
})

test_that("print.brsm_fit shows coding method when present", {
  mock <- .make_mock_brsm_fit(include_coding = TRUE)
  out  <- capture.output(print(mock))
  expect_true(any(grepl("coded", out)))
})

test_that("print.brsm_fit skips coding section when NULL", {
  mock <- .make_mock_brsm_fit(include_coding = FALSE)
  out  <- capture.output(print(mock))
  # Should not crash — just omit the coding block
  expect_true(any(grepl("Bayesian Response Surface Model", out)))
})

test_that("print.brsm_fit shows factor ranges when present", {
  mock <- .make_mock_brsm_fit(include_ranges = TRUE)
  out  <- capture.output(print(mock))
  expect_true(any(grepl("Factor Ranges", out)))
})

test_that("print.brsm_fit skips ranges section when NULL", {
  mock <- .make_mock_brsm_fit(include_ranges = FALSE)
  out  <- capture.output(print(mock))
  expect_false(any(grepl("Factor Ranges", out)))
})

test_that("print.brsm_fit shows MCMC sampling info when present", {
  mock <- .make_mock_brsm_fit(include_sampling = TRUE)
  out  <- capture.output(print(mock))
  expect_true(any(grepl("MCMC Sampling", out)))
  expect_true(any(grepl("4", out)))   # chains
})

test_that("print.brsm_fit skips sampling section when NULL", {
  mock <- .make_mock_brsm_fit(include_sampling = FALSE)
  out  <- capture.output(print(mock))
  expect_false(any(grepl("MCMC Sampling", out)))
})

test_that("print.brsm_fit invisibly returns the object", {
  mock <- .make_mock_brsm_fit()
  ret  <- suppressMessages(capture.output(res <- print(mock)))
  expect_identical(res, mock)
})

# ── summary.brsm_fit (error path when brms absent) ───────────────────────────

test_that("summary.brsm_fit errors when brms is not installed", {
  skip_if(requireNamespace("brms", quietly = TRUE),
          "brms is installed; skipping absence test")
  mock      <- .make_mock_brsm_fit()
  mock$fit  <- structure(list(), class = "brmsfit")
  expect_error(summary(mock), "brms")
})

# ── check_brsm_fit input validation (no real fit needed) ─────────────────────

test_that("check_brsm_fit errors on non-model objects or missing brms", {
  # check_brsm_fit() checks requireNamespace("brms") before object type,
  # so without brms the error is about the package, not the object class.
  expect_error(brsm::check_brsm_fit(list(a = 1)))
  expect_error(brsm::check_brsm_fit(42))
  expect_error(brsm::check_brsm_fit("not_a_model"))
  expect_error(brsm::check_brsm_fit(data.frame(x = 1)))
})

test_that("low-information helpers compute complexity and uncertainty summaries", {
  expect_equal(
    brsm:::.brsm_model_fixed_effect_count(c("x1", "x2"), "first_order"),
    3L
  )
  expect_equal(
    brsm:::.brsm_model_fixed_effect_count(c("x1", "x2"), "second_order"),
    6L
  )

  expect_null(brsm:::.brsm_low_information_message(60, 6, threshold = 10))
  expect_match(
    brsm:::.brsm_low_information_message(18, 6, threshold = 10),
    "regularized"
  )

  fake <- .create_fake_brmsfit(
    fixed = data.frame(),
    draws = data.frame(
      b_x1 = c(1, 2, 3, 4),
      b_x2 = c(-1, -2, -3, -4),
      check.names = FALSE
    )
  )
  uncertainty <- brsm:::.brsm_fixed_effect_uncertainty(fake)
  expect_equal(uncertainty$term, c("x1", "x2"))
  expect_equal(uncertainty$pd, c(1, 1))
  expect_false(any(uncertainty$overlap_zero))
})

test_that("internal check_brsm_fit helpers handle fake brmsfit objects", {
  fixed <- data.frame(
    Estimate = 1,
    Rhat = 1,
    Bulk_ESS = 500,
    Tail_ESS = 500,
    row.names = "b_x1"
  )
  fake <- .create_fake_brmsfit(
    fixed = fixed,
    metadata = function() list(max_treedepth = 12)
  )
  wrapped <- .make_mock_brsm_fit()
  wrapped$fit <- fake

  expect_identical(brsm:::.brsm_extract_fit(fake, caller = "test"), fake)
  expect_identical(brsm:::.brsm_extract_fit(wrapped, caller = "test"), fake)
  expect_error(
    brsm:::.brsm_extract_fit(.make_mock_brsm_fit(), caller = "test"),
    "object must be a brsm_fit or brmsfit"
  )

  expect_equal(brsm:::.brsm_get_max_treedepth(fake), 12L)
  expect_true(is.na(brsm:::.brsm_get_max_treedepth(.create_fake_brmsfit(fixed))))

  expect_equal(nrow(brsm:::.brsm_compute_bfmi(NULL)), 0)
  bfmi <- brsm:::.brsm_compute_bfmi(data.frame(
    Chain = c(1, 1, 1, 2, 2, 2),
    Parameter = rep("energy__", 6),
    Value = c(1, 2, 4, 5, 5, 5)
  ))
  expect_equal(bfmi$chain, c(1L, 2L))
  expect_true(is.finite(bfmi$bfmi[1]))
  expect_true(is.na(bfmi$bfmi[2]))
})

test_that("check_brsm_fit summarizes fake brmsfit diagnostics", {
  skip_if_not_installed("brms")

  fixed <- data.frame(
    Estimate = c(1, 2),
    Rhat = c(1.00, 1.03),
    Bulk_ESS = c(600, 300),
    Tail_ESS = c(700, 200),
    row.names = c("b_x1", "b_x2")
  )
  np <- data.frame(
    Chain = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2),
    Parameter = c(
      "divergent__", "divergent__", "treedepth__", "treedepth__",
      "energy__", "energy__", "energy__",
      "divergent__", "treedepth__", "treedepth__",
      "energy__", "energy__", "energy__", "energy__"
    ),
    Value = c(0, 1, 8, 9, 1, 2, 4, 0, 7, 8, 2, 4, 7, 11)
  )
  fake <- .create_fake_brmsfit(
    fixed = fixed,
    nuts_params = np,
    metadata = function() list(max_treedepth = 8)
  )

  expect_message(
    result <- brsm::check_brsm_fit(
      fake,
      bfmi_threshold = 1.5,
      verbose = TRUE
    ),
    "passed=FALSE"
  )

  expect_false(result$passed)
  expect_equal(result$overview$n_parameters, 2)
  expect_equal(result$overview$n_rhat_over_threshold, 1)
  expect_equal(result$overview$n_ess_bulk_below_min, 1)
  expect_equal(result$overview$n_ess_tail_below_min, 1)
  expect_equal(result$overview$divergences, 1)
  expect_equal(result$overview$treedepth_limit, 8)
  expect_equal(result$overview$n_max_treedepth_hits, 3)
  expect_equal(result$overview$n_bfmi_below_threshold, 2)
  expect_equal(result$parameters$parameter, c("b_x1", "b_x2"))
})

test_that("check_brsm_fit supports Eff.Sample fallback and validates inputs", {
  skip_if_not_installed("brms")

  fixed <- data.frame(
    Estimate = 1,
    Rhat = 1,
    Eff.Sample = 500,
    row.names = "b_x1"
  )
  fake <- .create_fake_brmsfit(fixed = fixed)

  result <- brsm::check_brsm_fit(fake, verbose = FALSE)
  expect_equal(result$overview$ess_bulk_min_observed, 500)
  expect_true(is.na(result$overview$ess_tail_min_observed))
  expect_true(is.na(result$overview$n_ess_tail_below_min))

  expect_error(
    brsm::check_brsm_fit(fake, treedepth_limit = "bad", verbose = FALSE),
    "treedepth_limit"
  )
  expect_error(
    brsm::check_brsm_fit(fake, bfmi_threshold = c(0.3, 0.4), verbose = FALSE),
    "bfmi_threshold"
  )
  expect_error(
    brsm::check_brsm_fit(.create_fake_brmsfit(fixed = NULL), verbose = FALSE),
    "No fixed-effect diagnostics"
  )
})
