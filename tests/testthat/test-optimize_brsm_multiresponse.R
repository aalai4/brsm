# Tests for optimize_brsm_multiresponse

.create_mro_draws_response1 <- function(n = 80, seed = 1001) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 8, 0.4),
    b_x1 = rnorm(n, 1.5, 0.1),
    b_x2 = rnorm(n, 0.3, 0.1),
    "b_I(x1^2)" = rnorm(n, -0.2, 0.03),
    "b_I(x2^2)" = rnorm(n, -0.1, 0.03),
    "b_x1:x2" = rnorm(n, 0.05, 0.02),
    sigma = abs(rnorm(n, 0.6, 0.05)),
    check.names = FALSE
  )
}

.create_mro_draws_response2 <- function(n = 80, seed = 1002) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n, 5, 0.4),
    b_x1 = rnorm(n, -1.0, 0.1),
    b_x2 = rnorm(n, 0.2, 0.1),
    "b_I(x1^2)" = rnorm(n, 0.15, 0.03),
    "b_I(x2^2)" = rnorm(n, 0.08, 0.03),
    "b_x1:x2" = rnorm(n, 0.00, 0.02),
    sigma = abs(rnorm(n, 0.5, 0.05)),
    check.names = FALSE
  )
}

test_that("optimize_brsm_multiresponse returns expected structure", {
  models <- list(
    y1 = .create_mro_draws_response1(),
    y2 = .create_mro_draws_response2()
  )

  specs <- list(
    y1 = list(goal = "max", low = 6, high = 12),
    y2 = list(goal = "min", low = 2, high = 8)
  )

  cand <- surface_grid(list(x1 = c(-2, 2), x2 = c(-2, 2)), n = 7)

  out <- optimize_brsm_multiresponse(
    models = models,
    desirability_specs = specs,
    factor_names = c("x1", "x2"),
    candidate_points = cand
  )

  expect_type(out, "list")
  expect_true("candidate_points" %in% names(out))
  expect_true("best_point" %in% names(out))
  expect_true("response_at_best" %in% names(out))
  expect_true("draw_count" %in% names(out))
  expect_true("desirability_specs" %in% names(out))
})

test_that("optimize_brsm_multiresponse builds grid from ranges", {
  models <- list(
    y1 = .create_mro_draws_response1(),
    y2 = .create_mro_draws_response2()
  )
  specs <- list(
    y1 = list(goal = "maximize", low = 6, high = 12),
    y2 = list(goal = "minimize", low = 2, high = 8)
  )

  out <- optimize_brsm_multiresponse(
    models = models,
    desirability_specs = specs,
    factor_names = c("x1", "x2"),
    ranges = list(x1 = c(-1, 1), x2 = c(-1, 1)),
    n_grid = 6
  )

  expect_equal(nrow(out$candidate_points), 36)
})

test_that("optimize_brsm_multiresponse can return draw-level desirability", {
  models <- list(
    y1 = .create_mro_draws_response1(n = 50),
    y2 = .create_mro_draws_response2(n = 50)
  )
  specs <- list(
    y1 = list(goal = "max", low = 6, high = 12, importance = 2),
    y2 = list(goal = "min", low = 2, high = 8, importance = 1)
  )
  cand <- surface_grid(list(x1 = c(-1, 1), x2 = c(-1, 1)), n = 5)

  out <- optimize_brsm_multiresponse(
    models = models,
    desirability_specs = specs,
    factor_names = c("x1", "x2"),
    candidate_points = cand,
    return_draws = TRUE,
    max_draws = 30
  )

  expect_true("combined_desirability_draws" %in% names(out))
  expect_true(is.matrix(out$combined_desirability_draws))
  expect_equal(dim(out$combined_desirability_draws), c(30, nrow(cand)))
})

test_that("optimize_brsm_multiresponse supports target desirability", {
  models <- list(
    y1 = .create_mro_draws_response1(),
    y2 = .create_mro_draws_response2()
  )
  specs <- list(
    y1 = list(goal = "target", low = 7, target = 9, high = 11),
    y2 = list(goal = "min", low = 2, high = 8)
  )
  cand <- surface_grid(list(x1 = c(-1, 1), x2 = c(-1, 1)), n = 5)

  out <- optimize_brsm_multiresponse(
    models = models,
    desirability_specs = specs,
    factor_names = c("x1", "x2"),
    candidate_points = cand
  )

  expect_true(nrow(out$best_point) == 1)
  expect_true(out$best_point$mean >= 0)
  expect_true(out$best_point$mean <= 1)
})

test_that("optimize_brsm_multiresponse validates key inputs", {
  models <- list(
    y1 = .create_mro_draws_response1(),
    y2 = .create_mro_draws_response2()
  )
  good_specs <- list(
    y1 = list(goal = "max", low = 6, high = 12),
    y2 = list(goal = "min", low = 2, high = 8)
  )
  cand <- surface_grid(list(x1 = c(-1, 1), x2 = c(-1, 1)), n = 4)

  expect_error(
    optimize_brsm_multiresponse(
      models = models,
      desirability_specs = good_specs,
      factor_names = c("x1", "x2")
    ),
    "Provide candidate_points or ranges"
  )

  expect_error(
    optimize_brsm_multiresponse(
      models = list(y1 = .create_mro_draws_response1()),
      desirability_specs = list(y1 = list(goal = "max", low = 6, high = 12)),
      factor_names = c("x1", "x2"),
      candidate_points = cand
    ),
    "at least two"
  )

  expect_error(
    optimize_brsm_multiresponse(
      models = models,
      desirability_specs = list(y1 = list(goal = "max", low = 6, high = 12)),
      factor_names = c("x1", "x2"),
      candidate_points = cand
    ),
    "Missing desirability specs"
  )

  expect_error(
    optimize_brsm_multiresponse(
      models = models,
      desirability_specs = list(
        y1 = list(goal = "target", low = 6, high = 12),
        y2 = list(goal = "min", low = 2, high = 8)
      ),
      factor_names = c("x1", "x2"),
      candidate_points = cand
    ),
    "Target desirability requires finite target"
  )
})