# Tests for print.brsm_fit using a mock brsm_fit object (no brms needed)

.make_mock_brsm_fit <- function(include_fit = FALSE,
                                 include_ranges = TRUE,
                                 include_sampling = TRUE,
                                 include_coding = TRUE) {
  # Create a mock draws data frame
  draws_df <- data.frame(
    b_Intercept = rnorm(10),
    b_x1 = rnorm(10),
    b_x2 = rnorm(10),
    b_x1_x2 = rnorm(10),
    b_I_x1_2 = rnorm(10),
    b_I_x2_2 = rnorm(10),
    sigma = runif(10)
  )

  fit_obj <- list(
    draws = draws_df,
    data = data.frame(x1 = 1:10, x2 = 1:10, y = 1:10)
  )
  class(fit_obj) <- "brsm_conjugate_fit"

  obj <- list(
    formula      = y ~ x1 + x2 + I(x1^2) + I(x2^2) + x1:x2,
    response     = "y",
    factor_names = c("x1", "x2"),
    coding       = if (include_coding) list(method = "coded") else NULL,
    ranges       = if (include_ranges)
                     list(x1 = c(-1, 1), x2 = c(-2, 2)) else NULL,
    sampling     = if (include_sampling)
                     list(draws = 2000) else NULL,
    fit          = if (include_fit) fit_obj else NULL,
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

test_that("print.brsm_fit shows Conjugate Draws info when present", {
  mock <- .make_mock_brsm_fit(include_sampling = TRUE)
  out  <- capture.output(print(mock))
  expect_true(any(grepl("Conjugate Draws", out)))
  expect_true(any(grepl("2000", out)))
})

test_that("print.brsm_fit skips sampling section when NULL", {
  mock <- .make_mock_brsm_fit(include_sampling = FALSE)
  out  <- capture.output(print(mock))
  expect_false(any(grepl("Conjugate Draws", out)))
})

test_that("print.brsm_fit invisibly returns the object", {
  mock <- .make_mock_brsm_fit()
  ret  <- suppressMessages(capture.output(res <- print(mock)))
  expect_identical(res, mock)
})
