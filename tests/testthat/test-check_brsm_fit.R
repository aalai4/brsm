# Tests for check_brsm_fit internal helpers and public function branches

# ---------------------------------------------------------------------------
# Fake class + summary method registered in .GlobalEnv so S3 dispatch finds it
# ---------------------------------------------------------------------------

# A minimal brmsfit subclass whose summary() is controlled via an attribute.
# nuts_params.brmsfit from brms will error on this fake → np = NULL.
.make_cbf_fake <- function(fixed_df) {
  obj <- structure(
    list(),
    class = c("cbf_fake_brmsfit", "brmsfit")
  )
  attr(obj, ".__fixed") <- fixed_df
  obj
}

# Register the summary method in the global environment so R's S3 dispatch
# finds it ahead of brms:::summary.brmsfit.
assign(
  "summary.cbf_fake_brmsfit",
  function(object, ...) list(fixed = attr(object, ".__fixed")),
  envir = .GlobalEnv
)

# ---------------------------------------------------------------------------
# .brsm_compute_bfmi
# ---------------------------------------------------------------------------

test_that(".brsm_compute_bfmi returns empty data frame when no energy__ rows", {
  np_no_energy <- data.frame(
    Chain     = c(1L, 1L),
    Parameter = c("divergent__", "treedepth__"),
    Value     = c(0.0, 5.0),
    stringsAsFactors = FALSE
  )
  result <- brsm:::.brsm_compute_bfmi(np_no_energy)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
  expect_true("chain" %in% names(result))
  expect_true("bfmi" %in% names(result))
})

test_that(".brsm_compute_bfmi computes BFMI for chains with energy__ rows", {
  np_with_energy <- data.frame(
    Chain     = c(1L, 1L, 1L, 2L, 2L, 2L),
    Parameter = rep("energy__", 6L),
    Value     = c(100.0, 102.0, 104.0, 200.0, 203.0, 207.0),
    stringsAsFactors = FALSE
  )
  result <- brsm:::.brsm_compute_bfmi(np_with_energy)
  expect_equal(nrow(result), 2L)
  expect_true(all(is.finite(result$bfmi)))
})

test_that(".brsm_compute_bfmi returns empty data frame for NULL or bad np", {
  expect_equal(nrow(brsm:::.brsm_compute_bfmi(NULL)), 0L)
  expect_equal(nrow(brsm:::.brsm_compute_bfmi(data.frame())), 0L)
})

# ---------------------------------------------------------------------------
# .brsm_get_max_treedepth
# ---------------------------------------------------------------------------

test_that(".brsm_get_max_treedepth returns NA when fit has no usable info", {
  # fit$fit is a plain list — @sim errors and no metadata function present
  fake_fit <- list(fit = list())
  td <- brsm:::.brsm_get_max_treedepth(fake_fit)
  expect_true(is.na(td))
})

test_that(".brsm_get_max_treedepth uses metadata() function when available", {
  fake_fit <- list(
    fit = list(
      metadata = function() list(max_treedepth = 12L)
    )
  )
  td <- brsm:::.brsm_get_max_treedepth(fake_fit)
  expect_equal(td, 12L)
})

test_that(".brsm_get_max_treedepth returns NA when metadata() lacks max_treedepth", {
  fake_fit <- list(
    fit = list(
      metadata = function() list(other_param = 1L)
    )
  )
  td <- brsm:::.brsm_get_max_treedepth(fake_fit)
  expect_true(is.na(td))
})

test_that(".brsm_get_max_treedepth returns NA when metadata() returns NA for the key", {
  fake_fit <- list(
    fit = list(
      metadata = function() list(max_treedepth = NA_integer_)
    )
  )
  td <- brsm:::.brsm_get_max_treedepth(fake_fit)
  expect_true(is.na(td))
})

# ---------------------------------------------------------------------------
# check_brsm_fit — public function branches
# ---------------------------------------------------------------------------

test_that("check_brsm_fit errors when fixed summary is NULL", {
  skip_if_not_installed("brms")
  fake <- .make_cbf_fake(NULL)
  expect_error(
    check_brsm_fit(fake),
    "No fixed-effect diagnostics found"
  )
})

test_that("check_brsm_fit errors when fixed summary has zero rows", {
  skip_if_not_installed("brms")
  fake <- .make_cbf_fake(data.frame())
  expect_error(
    check_brsm_fit(fake),
    "No fixed-effect diagnostics found"
  )
})

test_that("check_brsm_fit sets n_ess_tail_below_min to NA when Tail_ESS column absent", {
  skip_if_not_installed("brms")
  # Fixed df has Rhat and Bulk_ESS but NO Tail_ESS column
  fixed_df <- data.frame(
    Rhat     = c(1.001, 1.002),
    Bulk_ESS = c(500L, 600L),
    row.names = c("b_Intercept", "b_x1")
  )
  fake <- .make_cbf_fake(fixed_df)
  out <- suppressMessages(check_brsm_fit(fake))
  expect_true(is.na(out$overview$n_ess_tail_below_min))
})

test_that("check_brsm_fit uses Eff.Sample as fallback for bulk ESS column", {
  skip_if_not_installed("brms")
  fixed_df <- data.frame(
    Rhat       = c(1.001),
    Eff.Sample = c(550L),
    row.names  = c("b_Intercept")
  )
  fake <- .make_cbf_fake(fixed_df)
  out <- suppressMessages(check_brsm_fit(fake))
  # Eff.Sample column should be detected and used; 550 > default 400 → 0 bad
  expect_false(is.na(out$overview$n_ess_bulk_below_min))
  expect_equal(out$overview$n_ess_bulk_below_min, 0L)
})

test_that("check_brsm_fit sets treedepth_limit to NA when nuts_params returns NULL", {
  skip_if_not_installed("brms")
  fixed_df <- data.frame(
    Rhat     = c(1.001),
    Bulk_ESS = c(500L),
    Tail_ESS = c(500L),
    row.names = c("b_Intercept")
  )
  fake <- .make_cbf_fake(fixed_df)
  # nuts_params will error on cbf_fake_brmsfit → np = NULL
  # → treedepth_limit remains NULL → set to NA_real_ at line 212
  out <- suppressMessages(check_brsm_fit(fake))
  expect_true(is.na(out$overview$treedepth_limit))
})

test_that("check_brsm_fit emits verbose message when verbose = TRUE", {
  skip_if_not_installed("brms")
  fixed_df <- data.frame(
    Rhat     = c(1.001),
    Bulk_ESS = c(500L),
    Tail_ESS = c(500L),
    row.names = c("b_Intercept")
  )
  fake <- .make_cbf_fake(fixed_df)
  expect_message(
    check_brsm_fit(fake, verbose = TRUE),
    "check_brsm_fit:"
  )
})

test_that("check_brsm_fit returns passed=TRUE for clean diagnostics", {
  skip_if_not_installed("brms")
  fixed_df <- data.frame(
    Rhat     = c(1.001, 1.001),
    Bulk_ESS = c(800L, 900L),
    Tail_ESS = c(700L, 750L),
    row.names = c("b_Intercept", "b_x1")
  )
  fake <- .make_cbf_fake(fixed_df)
  out <- suppressMessages(check_brsm_fit(fake))
  expect_true(out$passed)
  expect_s3_class(out$overview, "data.frame")
  expect_equal(nrow(out$overview), 1L)
})

test_that("check_brsm_fit returns passed=FALSE when Rhat is too high", {
  skip_if_not_installed("brms")
  fixed_df <- data.frame(
    Rhat     = c(1.02),
    Bulk_ESS = c(800L),
    Tail_ESS = c(700L),
    row.names = c("b_Intercept")
  )
  fake <- .make_cbf_fake(fixed_df)
  out <- suppressMessages(check_brsm_fit(fake))
  expect_false(out$passed)
  expect_equal(out$overview$n_rhat_over_threshold, 1L)
})

test_that("check_brsm_fit validates treedepth_limit argument", {
  skip_if_not_installed("brms")
  fixed_df <- data.frame(
    Rhat = c(1.001),
    row.names = c("b_Intercept")
  )
  fake <- .make_cbf_fake(fixed_df)
  expect_error(
    check_brsm_fit(fake, treedepth_limit = "bad"),
    "treedepth_limit must be NULL or a finite numeric scalar"
  )
  expect_error(
    check_brsm_fit(fake, treedepth_limit = Inf),
    "treedepth_limit must be NULL or a finite numeric scalar"
  )
})
