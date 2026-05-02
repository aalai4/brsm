test_that("posterior_ridge_analysis attaches solver diagnostics", {
  draws <- data.frame(
    b_x1 = c(1, 1, 1),
    b_x2 = c(1, 1, 1),
    `b_I(x1^2)` = c(-1, -1, 0),
    `b_I(x2^2)` = c(-1, -1e-12, -1),
    `b_x1:x2` = c(0, 0, 0)
  )

  result <- suppressWarnings(
    posterior_ridge_analysis(
      draws,
      factor_names = c("x1", "x2"),
      radii = c(0, 1),
      summary = TRUE
    )
  )

  diag_info <- attr(result, "diagnostics", exact = TRUE)

  expect_false(is.null(diag_info))
  expect_true(all(c(
    "status_code",
    "status_label",
    "kappa_proxy",
    "status_counts",
    "function_name"
  ) %in% names(diag_info)))
  expect_identical(diag_info$function_name, "posterior_ridge_analysis")
  expect_equal(diag_info$n_draws, 3)
})

test_that("canonical_analysis attaches solver diagnostics when computing scores", {
  draws <- data.frame(
    b_x1 = c(1, 1, 1),
    b_x2 = c(1, 1, 1),
    `b_I(x1^2)` = c(-1, -1, 0),
    `b_I(x2^2)` = c(-1, -1e-12, -1),
    `b_x1:x2` = c(0, 0, 0)
  )

  result <- suppressWarnings(
    canonical_analysis(
      draws,
      factor_names = c("x1", "x2"),
      include_scores = TRUE,
      kappa_thresh = 1e10,
      summary = TRUE
    )
  )

  diag_info <- attr(result, "diagnostics", exact = TRUE)

  expect_false(is.null(diag_info))
  expect_true(all(c(
    "status_code",
    "status_label",
    "kappa_proxy",
    "status_counts",
    "function_name"
  ) %in% names(diag_info)))
  expect_identical(diag_info$function_name, "canonical_analysis")
  expect_equal(diag_info$n_draws, 3)
})

test_that("summarize_brsm_stability extracts and aggregates diagnostics", {
  draws <- data.frame(
    b_x1 = c(1, 1, 1),
    b_x2 = c(1, 1, 1),
    `b_I(x1^2)` = c(-1, -1, 0),
    `b_I(x2^2)` = c(-1, -1e-12, -1),
    `b_x1:x2` = c(0, 0, 0)
  )

  sp_result <- suppressWarnings(
    stationary_point(
      draws,
      factor_names = c("x1", "x2"),
      kappa_thresh = 1e10,
      diagnostics = "basic",
      auto_guidance = FALSE
    )
  )

  ridge_result <- suppressWarnings(
    posterior_ridge_analysis(
      draws,
      factor_names = c("x1", "x2"),
      radii = c(0, 1),
      summary = TRUE
    )
  )

  summary_df <- summarize_brsm_stability(sp_result, ridge_result)

  expect_s3_class(summary_df, "data.frame")
  expect_equal(nrow(summary_df), 2)
  expect_true(all(c(
    "function_name",
    "n_draws",
    "n_excluded",
    "pct_excluded",
    "status_ok",
    "status_lapack_fail",
    "status_invalid_lu_diag",
    "status_kappa_exceeded"
  ) %in% names(summary_df)))
  expect_equal(summary_df$function_name, c("stationary_point", "posterior_ridge_analysis"))
  expect_true(all(summary_df$n_draws == 3))
})

test_that("summarize_brsm_stability.default errors when no diagnostics found", {
  expect_error(
    summarize_brsm_stability(data.frame(x = 1)),
    "No diagnostics found in any supplied objects"
  )
  expect_error(
    summarize_brsm_stability(list(a = 1), list(b = 2)),
    "No diagnostics found"
  )
})

test_that("summarize_brsm_stability.brsm_fit handles no diagnostics and with extras", {
  fake_fit <- structure(list(factor_names = "x1"), class = "brsm_fit")
  expect_message(
    result <- brsm::summarize_brsm_stability(fake_fit),
    "No diagnostics found"
  )
  expect_null(result)

  draws <- data.frame(
    b_x1 = c(1, 1, 1),
    b_x2 = c(1, 1, 1),
    `b_I(x1^2)` = c(-1, -1, 0),
    `b_I(x2^2)` = c(-1, -1e-12, -1),
    `b_x1:x2` = c(0, 0, 0)
  )
  ridge_result <- suppressWarnings(
    posterior_ridge_analysis(draws, factor_names = c("x1", "x2"), radii = c(0, 1), summary = TRUE)
  )
  summary_df <- brsm::summarize_brsm_stability(fake_fit, ridge_result)
  expect_s3_class(summary_df, "data.frame")
  expect_equal(nrow(summary_df), 1)

  attr(fake_fit, "diagnostics") <- attr(ridge_result, "diagnostics", exact = TRUE)
  result_with_diag <- brsm::summarize_brsm_stability(fake_fit)
  expect_s3_class(result_with_diag, "data.frame")
  expect_equal(nrow(result_with_diag), 1)
})

test_that("posterior_ridge_analysis dispatch and validation branches are covered", {
  # brsm_fit dispatch branch (factor_names pulled from object).
  draws <- data.frame(
    b_Intercept = c(1, 1, 1),
    b_x1 = c(0.5, 0.6, 0.7),
    b_x2 = c(-0.2, -0.3, -0.1),
    `b_I(x1^2)` = c(-1, -0.9, -1.1),
    `b_I(x2^2)` = c(-1, -1.1, -0.8),
    `b_x1:x2` = c(0.1, 0.1, 0.05),
    check.names = FALSE
  )

  fit_obj <- structure(
    list(.draws_df = draws),
    class = c("fake_brmsfit_ridge", "brmsfit")
  )
  as.data.frame.fake_brmsfit_ridge <- function(x, ...) x$.draws_df
  registerS3method(
    "as.data.frame",
    "fake_brmsfit_ridge",
    as.data.frame.fake_brmsfit_ridge,
    envir = asNamespace("base")
  )

  mock <- structure(
    list(fit = fit_obj, factor_names = c("x1", "x2"), model_terms = "second_order"),
    class = "brsm_fit"
  )
  out_fit <- suppressWarnings(posterior_ridge_analysis(mock, radii = c(0, 0.5), summary = FALSE))
  expect_s3_class(out_fit, "data.frame")
  expect_true(all(c("draw", "radius", "x1", "x2") %in% names(out_fit)))

  # Input validation branches.
  expect_error(
    posterior_ridge_analysis(data.frame(b_x1 = 1), radii = c(0, 1)),
    "factor_names must be supplied"
  )
  expect_error(
    posterior_ridge_analysis(draws, factor_names = c("x1", "x2"), radii = c(-1, 0)),
    "radii must be a non-empty"
  )
  expect_error(
    posterior_ridge_analysis(draws, factor_names = c("x1", "x2"), tol = 0),
    "tol must be a positive"
  )
})

test_that("posterior_ridge_analysis missing-coef and NA-solver branches warn as expected", {
  # Missing quadratic and interaction columns triggers missing coefficient warning.
  sparse_draws <- data.frame(
    b_x1 = c(0.2, 0.3),
    b_x2 = c(-0.1, -0.2),
    check.names = FALSE
  )
  expect_warning(
    expect_error(
      posterior_ridge_analysis(
        sparse_draws,
        factor_names = c("x1", "x2"),
        radii = c(0, 1),
        summary = FALSE
      ),
      "missing quadratic columns"
    ),
    "Missing coefficients"
  )

  # All-NA draws exercise NA root-solving path and all-NA summary branch.
  na_draws <- data.frame(
    b_x1 = c(NA_real_, NA_real_),
    b_x2 = c(NA_real_, NA_real_),
    `b_I(x1^2)` = c(NA_real_, NA_real_),
    `b_I(x2^2)` = c(NA_real_, NA_real_),
    `b_x1:x2` = c(NA_real_, NA_real_),
    check.names = FALSE
  )

  out_na <- suppressWarnings(
    posterior_ridge_analysis(
      na_draws,
      factor_names = c("x1", "x2"),
      radii = c(0, 1),
      summary = TRUE
    )
  )
  expect_true(all(is.na(out_na$mean)))
  expect_true(all(is.na(out_na$q2.5)))
})