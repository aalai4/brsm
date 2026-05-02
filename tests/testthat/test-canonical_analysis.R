# Tests for canonical_analysis

.create_ca_draws <- function(n = 100, seed = 401) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 10, 1),
    b_x1 = rnorm(n, 0.5, 0.1),
    b_x2 = rnorm(n, -0.3, 0.1),
    "b_I(x1^2)" = rnorm(n, -0.8, 0.05),
    "b_I(x2^2)" = rnorm(n, -0.5, 0.05),
    "b_x1:x2" = rnorm(n, 0.2, 0.03),
    check.names = FALSE
  )
}

.create_ca_draws_3d <- function(n = 100, seed = 402) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 10, 1),
    b_x1 = rnorm(n, 0.5, 0.1),
    b_x2 = rnorm(n, -0.3, 0.1),
    b_x3 = rnorm(n, 0.1, 0.1),
    "b_I(x1^2)" = rnorm(n, -0.8, 0.05),
    "b_I(x2^2)" = rnorm(n, -0.5, 0.05),
    "b_I(x3^2)" = rnorm(n, -0.4, 0.05),
    "b_x1:x2" = rnorm(n, 0.2, 0.03),
    "b_x1:x3" = rnorm(n, -0.1, 0.03),
    "b_x2:x3" = rnorm(n, 0.05, 0.03),
    check.names = FALSE
  )
}

test_that("canonical_analysis returns expected structure for 2 factors", {
  draws <- .create_ca_draws()
  ca <- canonical_analysis(draws, factor_names = c("x1", "x2"))

  expect_type(ca, "list")
  expect_true("eigenvalues" %in% names(ca))
  expect_true("eigenvectors" %in% names(ca))
  expect_true("scores" %in% names(ca))
})

test_that("canonical_analysis eigenvalue summary has correct dimensions", {
  draws <- .create_ca_draws()
  ca <- canonical_analysis(draws, factor_names = c("x1", "x2"))

  # 2 factors -> 2 canonical axes
  expect_equal(nrow(ca$eigenvalues), 2)
  expect_true("axis" %in% names(ca$eigenvalues))
  expect_true("label" %in% names(ca$eigenvalues))
  expect_true("mean" %in% names(ca$eigenvalues))
  expect_true("sd" %in% names(ca$eigenvalues))
})

test_that("canonical_analysis eigenvector summary has correct dimensions", {
  draws <- .create_ca_draws()
  ca <- canonical_analysis(draws, factor_names = c("x1", "x2"))

  # 2 factors × 2 axes = 4 rows
  expect_equal(nrow(ca$eigenvectors), 4)
  expect_true("axis" %in% names(ca$eigenvectors))
  expect_true("factor" %in% names(ca$eigenvectors))
  expect_true("mean" %in% names(ca$eigenvectors))
  expect_equal(sort(unique(ca$eigenvectors$factor)), c("x1", "x2"))
})

test_that("canonical_analysis scores summary has correct dimensions", {
  draws <- .create_ca_draws()
  ca <- canonical_analysis(draws, factor_names = c("x1", "x2"))

  # 2 canonical scores
  expect_equal(nrow(ca$scores), 2)
  expect_true("axis" %in% names(ca$scores))
  expect_true("label" %in% names(ca$scores))
  expect_equal(ca$scores$label, c("z_1", "z_2"))
})

test_that("canonical_analysis works for 3 factors", {
  draws <- .create_ca_draws_3d()
  ca <- canonical_analysis(draws, factor_names = c("x1", "x2", "x3"))

  expect_equal(nrow(ca$eigenvalues), 3)
  expect_equal(nrow(ca$eigenvectors), 9)
  expect_equal(nrow(ca$scores), 3)
})

test_that("canonical_analysis summary=FALSE returns raw matrices", {
  draws <- .create_ca_draws(n = 50)
  ca <- canonical_analysis(draws, factor_names = c("x1", "x2"),
                           summary = FALSE)

  expect_type(ca, "list")
  expect_true(is.matrix(ca$eigenvalues))
  expect_equal(dim(ca$eigenvalues), c(50, 2))
  expect_true(is.array(ca$eigenvectors))
  expect_equal(dim(ca$eigenvectors), c(50, 2, 2))
  expect_true(is.matrix(ca$scores))
  expect_equal(dim(ca$scores), c(50, 2))
})

test_that("canonical_analysis include_scores=FALSE omits scores", {
  draws <- .create_ca_draws()
  ca <- canonical_analysis(draws, factor_names = c("x1", "x2"),
                           include_scores = FALSE)

  expect_false("scores" %in% names(ca))
})

test_that("canonical_analysis eigenvalues are ordered descending on average", {
  # For a pure maximum (all negative curvatures), lambda_1 should be the
  # least negative (largest). eigen() orders by decreasing value.
  draws <- .create_ca_draws()
  ca <- canonical_analysis(draws, factor_names = c("x1", "x2"),
                           summary = FALSE)

  mean_ev <- colMeans(ca$eigenvalues, na.rm = TRUE)
  expect_true(mean_ev[1] >= mean_ev[2])
})

test_that("canonical_analysis eigenvalues are negative for a clear maximum", {
  # The mock draws have strongly negative diagonal Hessian so all
  # eigenvalues should be negative on average
  draws <- .create_ca_draws()
  ca <- canonical_analysis(draws, factor_names = c("x1", "x2"))

  expect_true(all(ca$eigenvalues$mean < 0))
})

test_that("canonical_analysis validates include_scores argument", {
  draws <- .create_ca_draws(n = 10)

  expect_error(
    canonical_analysis(draws, factor_names = c("x1", "x2"),
                       include_scores = NA),
    "include_scores must be TRUE or FALSE"
  )
})

test_that("canonical_analysis validates kappa_thresh argument", {
  draws <- .create_ca_draws(n = 10)

  expect_error(
    canonical_analysis(draws, factor_names = c("x1", "x2"),
                       kappa_thresh = -1),
    "kappa_thresh must be a finite positive numeric scalar"
  )
})

test_that("canonical_analysis validates summary argument", {
  draws <- .create_ca_draws(n = 10)

  expect_error(
    canonical_analysis(draws, factor_names = c("x1", "x2"),
                       summary = NA),
    "summary must be TRUE or FALSE"
  )
})

test_that("canonical_analysis validates probs argument", {
  draws <- .create_ca_draws(n = 10)

  expect_error(
    canonical_analysis(draws, factor_names = c("x1", "x2"),
                       probs = c(-0.1, 0.5)),
    "probs must be finite values between 0 and 1"
  )
})

test_that("canonical_analysis respects custom probs", {
  draws <- .create_ca_draws()
  ca <- canonical_analysis(draws, factor_names = c("x1", "x2"),
                           probs = c(0.1, 0.9))

  expect_true("q10.0" %in% names(ca$eigenvalues))
  expect_true("q90.0" %in% names(ca$eigenvalues))
})

test_that("canonical_analysis errors on missing factor_names for data frame", {
  draws <- .create_ca_draws(n = 10)

  expect_error(
    canonical_analysis(draws),
    "factor_names must be supplied"
  )
})

test_that("canonical_analysis works via brsm_fit dispatch", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 403)
  fit <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    coding_policy = "ignore",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 403,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  ca <- canonical_analysis(fit)

  expect_type(ca, "list")
  expect_equal(nrow(ca$eigenvalues), 2)
  expect_equal(nrow(ca$eigenvectors), 4)
  expect_equal(nrow(ca$scores), 2)
})

.replace_brsm_ns_binding <- function(name, new_fn) {
  ns <- asNamespace("brsm")
  old <- get(name, envir = ns)
  unlockBinding(name, ns)
  assign(name, new_fn, envir = ns)
  lockBinding(name, ns)
  old
}

.restore_brsm_ns_binding <- function(name, old_fn) {
  ns <- asNamespace("brsm")
  unlockBinding(name, ns)
  assign(name, old_fn, envir = ns)
  lockBinding(name, ns)
}

test_that("canonical_analysis.brsm_fit dispatches with object factor_names", {
  fake_draws <- .create_ca_draws(n = 4)
  fake_fit <- structure(list(factor_names = c("x1", "x2")), class = "brsm_fit")

  old_as <- .replace_brsm_ns_binding("as_brsm_draws", function(object, ...) fake_draws)
  old_default <- .replace_brsm_ns_binding("canonical_analysis.default", function(object, factor_names = NULL, include_scores = TRUE, kappa_thresh = 1e10, probs = c(0.025, 0.5, 0.975), summary = TRUE) {
    list(
      factor_names = factor_names,
      include_scores = include_scores,
      kappa_thresh = kappa_thresh,
      probs = probs,
      summary = summary,
      n_draws = nrow(object)
    )
  })
  on.exit({
    .restore_brsm_ns_binding("canonical_analysis.default", old_default)
    .restore_brsm_ns_binding("as_brsm_draws", old_as)
  }, add = TRUE)

  out <- brsm:::canonical_analysis.brsm_fit(
    fake_fit,
    include_scores = FALSE,
    kappa_thresh = 123,
    probs = c(0.1, 0.9),
    summary = FALSE
  )

  expect_identical(out$factor_names, c("x1", "x2"))
  expect_false(out$include_scores)
  expect_identical(out$kappa_thresh, 123)
  expect_identical(out$probs, c(0.1, 0.9))
  expect_false(out$summary)
  expect_identical(out$n_draws, 4L)
})

test_that("canonical_analysis.brmsfit validates and dispatches with supplied factors", {
  fake_fit <- structure(list(), class = "brmsfit")
  expect_error(
    brsm:::canonical_analysis.brmsfit(fake_fit),
    "factor_names must be supplied"
  )

  fake_draws <- .create_ca_draws(n = 3)
  old_as <- .replace_brsm_ns_binding("as_brsm_draws", function(object, factor_names = NULL, ...) {
    attr(fake_draws, "seen_factor_names") <- factor_names
    fake_draws
  })
  old_default <- .replace_brsm_ns_binding("canonical_analysis.default", function(object, factor_names = NULL, ...) {
    list(n_draws = nrow(object), factor_names = factor_names)
  })
  on.exit({
    .restore_brsm_ns_binding("canonical_analysis.default", old_default)
    .restore_brsm_ns_binding("as_brsm_draws", old_as)
  }, add = TRUE)

  out <- brsm:::canonical_analysis.brmsfit(fake_fit, factor_names = c("x1", "x2"))
  expect_identical(out$n_draws, 3L)
  expect_identical(out$factor_names, c("x1", "x2"))
})

test_that("canonical_analysis handles decomposition failures with NA summaries", {
  draws <- .create_ca_draws(n = 5)

  old_h <- .replace_brsm_ns_binding(".brsm_hessian_array", function(draws, factor_names) {
    array(NA_real_, dim = c(nrow(draws), length(factor_names), length(factor_names)))
  })
  old_sp <- .replace_brsm_ns_binding("stationary_points_batch_details", function(h_array, b_matrix, kappa_thresh) {
    n <- dim(h_array)[1]
    p <- dim(h_array)[2]
    list(
      x_star = matrix(0, nrow = n, ncol = p),
      status_code = rep(1L, n),
      kappa_proxy = rep(Inf, n)
    )
  })
  on.exit({
    .restore_brsm_ns_binding("stationary_points_batch_details", old_sp)
    .restore_brsm_ns_binding(".brsm_hessian_array", old_h)
  }, add = TRUE)

  out <- expect_warning(
    expect_warning(
      brsm:::canonical_analysis.default(
        draws,
        factor_names = c("x1", "x2"),
        include_scores = TRUE,
        summary = TRUE
      ),
      "draw\\(s\\) could not be decomposed"
    ),
    "draw\\(s\\) produced NA canonical scores"
  )

  expect_true(all(is.na(out$eigenvalues$mean)))
  expect_true(all(is.na(out$eigenvectors$mean)))
  expect_true(all(is.na(out$scores$mean)))
})