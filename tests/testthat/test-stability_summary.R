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