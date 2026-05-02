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
    predict_surface(draws, c("x1", "x2"), grid, draw_subset = numeric(0)),
    "one or more finite indices"
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

test_that("predict_surface covers dispatch and validation branches", {
  full_draws <- .create_predict_draws(n = 30)
  grid <- .create_predict_grid()

  as.data.frame.fake_brmsfit_ps <<- function(x, ...) x$.draws_df
  registerS3method(
    "as.data.frame", "fake_brmsfit_ps",
    as.data.frame.fake_brmsfit_ps,
    envir = asNamespace("base")
  )

  fake_brmsfit <- structure(
    list(.draws_df = full_draws),
    class = c("fake_brmsfit_ps", "brmsfit")
  )

  fake_brsm_fit <- structure(
    list(
      fit = fake_brmsfit,
      factor_names = c("x1", "x2"),
      model_terms = "second_order"
    ),
    class = "brsm_fit"
  )

  out_fit <- predict_surface(fake_brsm_fit, newdata = grid, summary = TRUE)
  out_brms <- predict_surface(fake_brmsfit,
    factor_names = c("x1", "x2"),
    newdata = grid,
    summary = TRUE
  )
  expect_equal(nrow(out_fit), nrow(grid))
  expect_equal(nrow(out_brms), nrow(grid))

  expect_error(
    predict_surface(full_draws, c("x1", "x2"), as.matrix(grid)),
    "newdata must be a data.frame"
  )
  expect_error(
    predict_surface(full_draws, c("x1", "x2"), grid[0, , drop = FALSE]),
    "newdata must contain at least one row"
  )
  expect_error(
    predict_surface(full_draws, c("x1", "x2"), grid, draw_subset = "bad"),
    "draw_subset must be NULL, a logical vector, or numeric indices"
  )

  none_selected <- rep(FALSE, nrow(full_draws))
  expect_error(
    predict_surface(full_draws, c("x1", "x2"), grid, draw_subset = none_selected),
    "No draws remain after applying draw_subset/max_draws"
  )

  draws_no_intercept <- full_draws
  draws_no_intercept$b_Intercept <- NULL
  expect_error(
    predict_surface(draws_no_intercept, c("x1", "x2"), grid),
    "draws must contain 'b_Intercept'"
  )

  grid_char <- grid
  grid_char$x1 <- as.character(grid_char$x1)
  expect_warning(
    predict_surface(full_draws, c("x1", "x2"), grid_char, summary = FALSE),
    "coerced to numeric"
  )

  expect_warning(
    predict_surface(full_draws, c("x1", "x2"), grid, summary = TRUE, return_matrix = TRUE),
    "return_matrix is ignored when summary = TRUE"
  )
})

test_that("predict_surface enters chunked long-output path", {
  set.seed(903)
  n_draws <- 2300
  n_points <- 2301

  draws <- data.frame(
    b_Intercept = rnorm(n_draws),
    b_x1 = rnorm(n_draws),
    check.names = FALSE
  )
  grid <- data.frame(x1 = seq(-1, 1, length.out = n_points))

  out <- predict_surface(
    draws = draws,
    factor_names = c("x1"),
    newdata = grid,
    summary = FALSE,
    output_chunk_size = 128
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), n_draws * n_points)
  expect_true(all(c("draw", "point_id", "x1", "estimate") %in% names(out)))
})
