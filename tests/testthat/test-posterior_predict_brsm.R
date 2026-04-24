# Tests for posterior_predict_brsm

.create_postpred_draws <- function(n = 60, seed = 921) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 10, 1),
    b_x1 = rnorm(n, 1, 0.2),
    b_x2 = rnorm(n, -0.5, 0.2),
    "b_I(x1^2)" = rnorm(n, -0.2, 0.05),
    "b_I(x2^2)" = rnorm(n, -0.1, 0.05),
    "b_x1:x2" = rnorm(n, 0.05, 0.02),
    sigma = abs(rnorm(n, 1.2, 0.1)),
    check.names = FALSE
  )
}

.create_postpred_grid <- function() {
  surface_grid(
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    n = 6
  )
}

test_that("posterior_predict_brsm returns matrix when requested", {
  draws <- .create_postpred_draws(n = 25)
  grid <- .create_postpred_grid()

  mat <- posterior_predict_brsm(
    object = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    return_matrix = TRUE
  )

  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(25, nrow(grid)))
})

test_that("posterior_predict_brsm include_residual=FALSE matches mean surface", {
  draws <- .create_postpred_draws(n = 30)
  grid <- .create_postpred_grid()

  pp <- posterior_predict_brsm(
    object = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    include_residual = FALSE,
    return_matrix = TRUE
  )

  mu <- predict_surface(
    draws = as_brsm_draws.data.frame(
      draws,
      factor_names = c("x1", "x2"),
      require_quadratic = FALSE,
      require_interactions = FALSE
    ),
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = FALSE,
    return_matrix = TRUE
  )

  expect_equal(pp, mu)
})

test_that("posterior_predict_brsm summary returns expected columns", {
  draws <- .create_postpred_draws(n = 35)
  grid <- .create_postpred_grid()

  out <- posterior_predict_brsm(
    object = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = TRUE,
    probs = c(0.1, 0.9)
  )

  expect_equal(nrow(out), nrow(grid))
  expect_true("mean" %in% names(out))
  expect_true("q10.0" %in% names(out))
  expect_true("q90.0" %in% names(out))
})

test_that("posterior_predict_brsm supports sigma override scalar", {
  draws <- .create_postpred_draws(n = 20)
  grid <- .create_postpred_grid()

  set.seed(1)
  out1 <- posterior_predict_brsm(
    object = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    sigma = 0.5,
    return_matrix = TRUE,
    seed = 101
  )
  out2 <- posterior_predict_brsm(
    object = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    sigma = 0.5,
    return_matrix = TRUE,
    seed = 101
  )

  expect_equal(out1, out2)
})

test_that("posterior_predict_brsm supports sigma override vector", {
  draws <- .create_postpred_draws(n = 22)
  grid <- .create_postpred_grid()
  sigma_vec <- rep(0.25, 22)

  out <- posterior_predict_brsm(
    object = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    sigma = sigma_vec,
    return_matrix = TRUE,
    seed = 5
  )

  expect_true(is.matrix(out))
  expect_equal(nrow(out), 22)
})

test_that("posterior_predict_brsm errors when sigma missing with residual", {
  draws <- .create_postpred_draws(n = 15)
  draws$sigma <- NULL
  grid <- .create_postpred_grid()

  expect_error(
    posterior_predict_brsm(
      object = draws,
      factor_names = c("x1", "x2"),
      newdata = grid,
      include_residual = TRUE
    ),
    "No posterior sigma column found"
  )
})

test_that("posterior_predict_brsm validates sigma argument", {
  draws <- .create_postpred_draws(n = 15)
  grid <- .create_postpred_grid()

  expect_error(
    posterior_predict_brsm(
      object = draws,
      factor_names = c("x1", "x2"),
      newdata = grid,
      sigma = -1
    ),
    "sigma must be a finite non-negative"
  )

  expect_error(
    posterior_predict_brsm(
      object = draws,
      factor_names = c("x1", "x2"),
      newdata = grid,
      sigma = rep(0.2, 3)
    ),
    "sigma vector length must be 1 or match"
  )
})

test_that("posterior_predict_brsm respects draw controls", {
  draws <- .create_postpred_draws(n = 40)
  grid <- .create_postpred_grid()

  out <- posterior_predict_brsm(
    object = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    draw_subset = c(2, 3, 9, 11),
    max_draws = 3,
    summary = FALSE
  )

  expect_equal(sort(unique(out$draw)), c(2, 3, 9))
  expect_equal(nrow(out), 3 * nrow(grid))
})

test_that("posterior_predict_brsm validates include_residual", {
  draws <- .create_postpred_draws(n = 10)
  grid <- .create_postpred_grid()

  expect_error(
    posterior_predict_brsm(
      object = draws,
      factor_names = c("x1", "x2"),
      newdata = grid,
      include_residual = NA
    ),
    "include_residual must be TRUE or FALSE"
  )
})

test_that("posterior_predict_brsm works with non-quadratic draws", {
  draws <- .create_postpred_draws(n = 20)
  draws[["b_I(x1^2)"]] <- NULL
  draws[["b_I(x2^2)"]] <- NULL
  draws[["b_x1:x2"]] <- NULL
  grid <- .create_postpred_grid()

  out <- posterior_predict_brsm(
    object = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = TRUE
  )

  expect_equal(nrow(out), nrow(grid))
  expect_true("mean" %in% names(out))
})