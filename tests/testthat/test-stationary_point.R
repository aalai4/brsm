test_that("stationary_point diagnostics include per-draw solver details", {
  draws <- data.frame(
    b_x1 = c(1, 1, 1),
    b_x2 = c(1, 1, 1),
    `b_I(x1^2)` = c(-1, -1, 0),
    `b_I(x2^2)` = c(-1, -1e-12, -1),
    `b_x1:x2` = c(0, 0, 0)
  )

  expect_warning(
    stationary_point(
      draws,
      factor_names = c("x1", "x2"),
      kappa_thresh = 1e10,
      diagnostics = "basic",
      auto_guidance = FALSE
    ),
    "near-singular Hessians"
  )

  result <- suppressWarnings(
    stationary_point(
      draws,
      factor_names = c("x1", "x2"),
      kappa_thresh = 1e10,
      diagnostics = "basic",
      auto_guidance = FALSE
    )
  )

  diag_info <- attr(result, "diagnostics", exact = TRUE)

  expect_false(is.null(diag_info))
  expect_true(all(c(
    "status_code",
    "status_label",
    "kappa_proxy",
    "status_counts"
  ) %in% names(diag_info)))

  expect_length(diag_info$status_code, 3)
  expect_length(diag_info$status_label, 3)
  expect_length(diag_info$kappa_proxy, 3)
  expect_named(
    diag_info$status_counts,
    c("ok", "lapack_fail", "invalid_lu_diag", "kappa_exceeded")
  )

  expect_identical(diag_info$status_label[[1]], "ok")
  expect_identical(diag_info$status_label[[2]], "kappa_exceeded")
  expect_identical(diag_info$status_label[[3]], "lapack_fail")

  expect_true(is.finite(diag_info$kappa_proxy[[1]]))
  expect_true(is.finite(diag_info$kappa_proxy[[2]]))
  expect_true(is.na(diag_info$kappa_proxy[[3]]))

  expect_equal(unname(diag_info$status_counts), c(1L, 1L, 0L, 1L))
  expect_equal(diag_info$n_excluded, 2)
})

test_that("decode_stationary_status returns ordered labels", {
  decoded <- decode_stationary_status(c(0L, 1L, 2L, 3L, NA_integer_, 99L))

  expect_s3_class(decoded, "factor")
  expect_true(is.ordered(decoded))
  expect_identical(
    levels(decoded),
    c("ok", "lapack_fail", "invalid_lu_diag", "kappa_exceeded")
  )
  expect_identical(as.character(decoded[[1]]), "ok")
  expect_identical(as.character(decoded[[4]]), "kappa_exceeded")
  expect_true(is.na(decoded[[5]]))
  expect_true(is.na(decoded[[6]]))
})

test_that("stationary_point full diagnostics include status breakdown by threshold", {
  draws <- data.frame(
    b_x1 = c(1, 1, 1),
    b_x2 = c(1, 1, 1),
    `b_I(x1^2)` = c(-1, -1, 0),
    `b_I(x2^2)` = c(-1, -1e-12, -1),
    `b_x1:x2` = c(0, 0, 0)
  )

  result <- suppressWarnings(
    stationary_point(
      draws,
      factor_names = c("x1", "x2"),
      kappa_thresh = 1e10,
      diagnostics = "full",
      auto_guidance = FALSE,
      sensitivity_thresholds = c(1e8, 1e12)
    )
  )

  diag_info <- attr(result, "diagnostics", exact = TRUE)
  status_df <- diag_info$status_by_threshold

  expect_false(is.null(status_df))
  expect_identical(
    names(status_df),
    c(
      "threshold",
      "ok",
      "lapack_fail",
      "invalid_lu_diag",
      "kappa_exceeded",
      "n_excluded",
      "pct_excluded"
    )
  )
  expect_equal(status_df$threshold, c(1e8, 1e10, 1e12))
  expect_equal(status_df$ok, c(1, 1, 2))
  expect_equal(status_df$lapack_fail, c(1, 1, 1))
  expect_equal(status_df$invalid_lu_diag, c(0, 0, 0))
  expect_equal(status_df$kappa_exceeded, c(1, 1, 0))
  expect_equal(status_df$n_excluded, c(2, 2, 1))
  expect_equal(status_df$pct_excluded, c(2 / 3, 2 / 3, 1 / 3))
  expect_equal(diag_info$stability_flag, "unstable")
})