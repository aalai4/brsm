# Tests for predict_surface

.create_predict_draws <- function(n = 50, seed = 901) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 10, 1),
    b_x1 = rnorm(n, 1, 0.2),
    b_x2 = rnorm(n, -0.5, 0.2),
    "b_I(x1^2)" = rnorm(n, -0.2, 0.05),
    "b_I(x2^2)" = rnorm(n, -0.1, 0.05),
    "b_x1:x2" = rnorm(n, 0.05, 0.02),
    check.names = FALSE
  )
}

.create_predict_grid <- function() {
  surface_grid(
    ranges = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    n = 8
  )
}

test_that("predict_surface max_draws limits output rows", {
  draws <- .create_predict_draws(n = 40)
  grid <- .create_predict_grid()

  result <- predict_surface(
    draws = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = FALSE,
    max_draws = 10
  )

  expect_equal(length(unique(result$draw)), 10)
  expect_equal(nrow(result), 10 * nrow(grid))
})

test_that("predict_surface draw_subset supports numeric and logical", {
  draws <- .create_predict_draws(n = 20)
  grid <- .create_predict_grid()

  idx <- c(2, 5, 9)
  res_idx <- predict_surface(
    draws = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = FALSE,
    draw_subset = idx
  )

  mask <- rep(FALSE, nrow(draws))
  mask[idx] <- TRUE
  res_mask <- predict_surface(
    draws = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = FALSE,
    draw_subset = mask
  )

  expect_equal(sort(unique(res_idx$draw)), idx)
  expect_equal(nrow(res_idx), length(idx) * nrow(grid))
  expect_equal(res_idx$estimate, res_mask$estimate)
})

test_that("predict_surface return_matrix returns draw-by-point matrix", {
  draws <- .create_predict_draws(n = 15)
  grid <- .create_predict_grid()

  mat <- predict_surface(
    draws = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = FALSE,
    return_matrix = TRUE
  )

  df <- predict_surface(
    draws = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = FALSE
  )

  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(nrow(draws), nrow(grid)))
  expect_equal(as.vector(t(mat)), df$estimate)
})

test_that("predict_surface chunked output matches unchunked output", {
  draws <- .create_predict_draws(n = 23)
  grid <- .create_predict_grid()

  chunked <- predict_surface(
    draws = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = FALSE,
    output_chunk_size = 4
  )

  unchunked <- predict_surface(
    draws = draws,
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = FALSE,
    output_chunk_size = 1000
  )

  expect_equal(chunked$draw, unchunked$draw)
  expect_equal(chunked$point_id, unchunked$point_id)
  expect_equal(chunked$estimate, unchunked$estimate)
})

test_that("predict_surface validates new draw controls", {
  draws <- .create_predict_draws(n = 10)
  grid <- .create_predict_grid()

  expect_error(
    predict_surface(draws, c("x1", "x2"), grid, draw_subset = c(TRUE, FALSE)),
    "Logical draw_subset must have length"
  )

  expect_error(
    predict_surface(draws, c("x1", "x2"), grid, draw_subset = c(0, 2)),
    "out of bounds"
  )

  expect_error(
    predict_surface(draws, c("x1", "x2"), grid, max_draws = 0),
    "max_draws"
  )

  expect_error(
    predict_surface(draws, c("x1", "x2"), grid, output_chunk_size = 0),
    "output_chunk_size"
  )
})
