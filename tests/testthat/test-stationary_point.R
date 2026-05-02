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

test_that("stationary_point validates inputs and brsm_fit dispatch", {
  draws <- data.frame(
    b_x1 = c(1, 1, 1),
    b_x2 = c(1, 1, 1),
    `b_I(x1^2)` = c(-1, -1, -1),
    `b_I(x2^2)` = c(-1, -1, -1),
    `b_x1:x2` = c(0, 0, 0)
  )

  expect_error(
    stationary_point(draws),
    "factor_names must be supplied"
  )
  expect_error(
    stationary_point(draws, factor_names = c("x1", "x2"), kappa_thresh = 0),
    "kappa_thresh must be a finite positive"
  )
  expect_error(
    stationary_point(draws, factor_names = c("x1", "x2"), auto_guidance = NA),
    "auto_guidance must be TRUE or FALSE"
  )
  expect_error(
    stationary_point(
      draws,
      factor_names = c("x1", "x2"),
      sensitivity_thresholds = numeric(0)
    ),
    "sensitivity_thresholds must be a non-empty"
  )

  missing_linear <- draws
  missing_linear$b_x2 <- NULL
  expect_error(
    stationary_point(missing_linear, factor_names = c("x1", "x2")),
    "missing linear columns"
  )

  non_numeric_linear <- draws
  non_numeric_linear$b_x1 <- as.character(non_numeric_linear$b_x1)
  expect_error(
    stationary_point(non_numeric_linear, factor_names = c("x1", "x2")),
    "numeric columns for linear terms"
  )

  as.data.frame.fake_brmsfit_sp <<- function(x, ...) x$.draws_df
  registerS3method(
    "as.data.frame", "fake_brmsfit_sp",
    as.data.frame.fake_brmsfit_sp,
    envir = asNamespace("base")
  )

  fake_fit <- structure(
    list(.draws_df = data.frame(
      b_Intercept = 1:3,
      b_x1 = c(1, 1, 1),
      b_x2 = c(1, 1, 1),
      `b_I(x1^2)` = c(-1, -1, -1),
      `b_I(x2^2)` = c(-1, -1, -1),
      `b_x1:x2` = c(0, 0, 0),
      check.names = FALSE
    )),
    class = c("fake_brmsfit_sp", "brmsfit")
  )
  fit_obj <- structure(
    list(fit = fake_fit, factor_names = c("x1", "x2"), model_terms = "second_order"),
    class = "brsm_fit"
  )

  out <- stationary_point(fit_obj)
  expect_s3_class(out, "data.frame")
  expect_equal(names(out), c("x1", "x2"))
})

test_that("stationary_point guidance warnings cover moderate/high and stability flags", {
  # Moderate exclusion path: one excluded out of six (~16.7%).
  draws_mod <- data.frame(
    b_x1 = rep(1, 6),
    b_x2 = rep(1, 6),
    `b_I(x1^2)` = c(-1, -1, -1, -1, -1, 0),
    `b_I(x2^2)` = c(-1, -1, -1, -1, -1, -1),
    `b_x1:x2` = rep(0, 6)
  )
  mod_out <- suppressWarnings(stationary_point(
    draws_mod,
    factor_names = c("x1", "x2"),
    diagnostics = "basic",
    auto_guidance = TRUE
  ))
  expect_s3_class(mod_out, "data.frame")

  # Full diagnostics with one near-threshold draw produces moderate sensitivity.
  draws_stable <- data.frame(
    b_x1 = rep(1, 10),
    b_x2 = rep(1, 10),
    `b_I(x1^2)` = c(rep(-1, 9), -1),
    `b_I(x2^2)` = c(rep(-1, 9), -1e-11),
    `b_x1:x2` = rep(0, 10)
  )
  diag_mod <- suppressWarnings(stationary_point(
    draws_stable,
    factor_names = c("x1", "x2"),
    diagnostics = "full",
    auto_guidance = FALSE,
    sensitivity_thresholds = c(1e8, 1e12)
  ))
  expect_equal(attr(diag_mod, "diagnostics")$stability_flag, "moderately_sensitive")

  # Existing unstable geometry should trigger low-information warning in auto guidance.
  draws_unstable <- data.frame(
    b_x1 = c(1, 1, 1),
    b_x2 = c(1, 1, 1),
    `b_I(x1^2)` = c(-1, -1, 0),
    `b_I(x2^2)` = c(-1, -1e-12, -1),
    `b_x1:x2` = c(0, 0, 0)
  )
  unstable_out <- suppressWarnings(stationary_point(
    draws_unstable,
    factor_names = c("x1", "x2"),
    diagnostics = "full",
    auto_guidance = TRUE,
    sensitivity_thresholds = c(1e8, 1e12)
  ))
  expect_identical(attr(unstable_out, "diagnostics")$stability_flag, "unstable")

  # Stable sensitivity branch: no exclusions across thresholds.
  draws_stable_flag <- data.frame(
    b_x1 = rep(1, 8),
    b_x2 = rep(1, 8),
    `b_I(x1^2)` = rep(-1, 8),
    `b_I(x2^2)` = rep(-1, 8),
    `b_x1:x2` = rep(0, 8)
  )
  stable_out <- suppressWarnings(stationary_point(
    draws_stable_flag,
    factor_names = c("x1", "x2"),
    diagnostics = "full",
    auto_guidance = FALSE,
    sensitivity_thresholds = c(1e8, 1e12)
  ))
  expect_identical(attr(stable_out, "diagnostics")$stability_flag, "stable")
})