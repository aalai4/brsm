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
  expect_error(check_brsm_fit(list(a = 1)))
  expect_error(check_brsm_fit(42))
  expect_error(check_brsm_fit("not_a_model"))
  expect_error(check_brsm_fit(data.frame(x = 1)))
})
