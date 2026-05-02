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

# ---------------------------------------------------------------------------
# Additional branch coverage
# ---------------------------------------------------------------------------

test_that("optimize_brsm_multiresponse errors on unnamed models list", {
  unnamed_models <- list(
    .create_mro_draws_response1(),
    .create_mro_draws_response2()
  )
  expect_error(
    optimize_brsm_multiresponse(
      models = unnamed_models,
      desirability_specs = list(
        y1 = list(goal = "max", low = 6, high = 12),
        y2 = list(goal = "min", low = 2, high = 8)
      ),
      factor_names = c("x1", "x2"),
      candidate_points = surface_grid(list(x1 = c(-1, 1), x2 = c(-1, 1)), n = 3)
    ),
    "named list"
  )
})

test_that("optimize_brsm_multiresponse errors on return_draws = NA", {
  models <- list(
    y1 = .create_mro_draws_response1(),
    y2 = .create_mro_draws_response2()
  )
  specs <- list(
    y1 = list(goal = "max", low = 6, high = 12),
    y2 = list(goal = "min", low = 2, high = 8)
  )
  cand <- surface_grid(list(x1 = c(-1, 1), x2 = c(-1, 1)), n = 3)
  expect_error(
    optimize_brsm_multiresponse(
      models = models,
      desirability_specs = specs,
      factor_names = c("x1", "x2"),
      candidate_points = cand,
      return_draws = NA
    ),
    "return_draws must be TRUE or FALSE"
  )
})

test_that("optimize_brsm_multiresponse errors on non-data.frame candidate_points", {
  models <- list(
    y1 = .create_mro_draws_response1(),
    y2 = .create_mro_draws_response2()
  )
  specs <- list(
    y1 = list(goal = "max", low = 6, high = 12),
    y2 = list(goal = "min", low = 2, high = 8)
  )
  expect_error(
    optimize_brsm_multiresponse(
      models = models,
      desirability_specs = specs,
      factor_names = c("x1", "x2"),
      candidate_points = "not_a_df"
    ),
    "candidate_points must be a data.frame"
  )
})

test_that("optimize_brsm_multiresponse errors on zero-row candidate_points", {
  models <- list(
    y1 = .create_mro_draws_response1(),
    y2 = .create_mro_draws_response2()
  )
  specs <- list(
    y1 = list(goal = "max", low = 6, high = 12),
    y2 = list(goal = "min", low = 2, high = 8)
  )
  empty_cand <- data.frame(x1 = numeric(0), x2 = numeric(0))
  expect_error(
    optimize_brsm_multiresponse(
      models = models,
      desirability_specs = specs,
      factor_names = c("x1", "x2"),
      candidate_points = empty_cand
    ),
    "candidate_points must contain at least one row"
  )
})

test_that("optimize_brsm_multiresponse warns and truncates differing draw counts", {
  # y1 has 80 draws, y2 has 50; with max_draws=NULL both keep full counts → mismatch
  models <- list(
    y1 = .create_mro_draws_response1(n = 80),
    y2 = .create_mro_draws_response2(n = 50)
  )
  specs <- list(
    y1 = list(goal = "max", low = 6, high = 12),
    y2 = list(goal = "min", low = 2, high = 8)
  )
  cand <- surface_grid(list(x1 = c(-1, 1), x2 = c(-1, 1)), n = 3)

  expect_warning(
    out <- optimize_brsm_multiresponse(
      models = models,
      desirability_specs = specs,
      factor_names = c("x1", "x2"),
      candidate_points = cand,
      max_draws = NULL
    ),
    "truncating all to"
  )
  expect_equal(out$draw_count, 50L)
})

test_that(".brsm_validate_desirability_specs covers all validation branches", {
  vds <- brsm:::.brsm_validate_desirability_specs

  # Empty list
  expect_error(vds(list(), "y1"), "non-empty named list")

  # Unnamed (NULL names)
  expect_error(
    vds(list(list(goal = "max", low = 1, high = 2), list(goal = "min", low = 1, high = 2)),
        c("y1", "y2")),
    "named list"
  )

  # Spec element not a list
  expect_error(
    vds(list(y1 = "bad", y2 = list(goal = "min", low = 1, high = 2)), c("y1", "y2")),
    "Each desirability spec must be a list"
  )

  # Missing goal
  expect_error(
    vds(list(y1 = list(low = 1, high = 2), y2 = list(goal = "min", low = 1, high = 2)),
        c("y1", "y2")),
    "Each desirability spec must include goal"
  )

  # Invalid goal value
  expect_error(
    vds(list(y1 = list(goal = "banana", low = 1, high = 2),
             y2 = list(goal = "min", low = 1, high = 2)),
        c("y1", "y2")),
    "goal must be one of"
  )

  # low >= high
  expect_error(
    vds(list(y1 = list(goal = "max", low = 5, high = 3),
             y2 = list(goal = "min", low = 1, high = 2)),
        c("y1", "y2")),
    "finite low < high"
  )

  # Non-finite low
  expect_error(
    vds(list(y1 = list(goal = "max", low = -Inf, high = 3),
             y2 = list(goal = "min", low = 1, high = 2)),
        c("y1", "y2")),
    "finite low < high"
  )

  # target outside [low, high]
  expect_error(
    vds(list(y1 = list(goal = "target", low = 1, high = 5, target = 6),
             y2 = list(goal = "min", low = 1, high = 2)),
        c("y1", "y2")),
    "low < target < high"
  )

  # bad weight
  expect_error(
    vds(list(y1 = list(goal = "max", low = 1, high = 5, weight = -0.5),
             y2 = list(goal = "min", low = 1, high = 2)),
        c("y1", "y2")),
    "weight must be a finite positive scalar"
  )

  # bad importance
  expect_error(
    vds(list(y1 = list(goal = "max", low = 1, high = 5, importance = 0),
             y2 = list(goal = "min", low = 1, high = 2)),
        c("y1", "y2")),
    "importance must be a finite positive scalar"
  )

  # Happy path — should return a named list without error
  result <- vds(
    list(y1 = list(goal = "max", low = 1, high = 5),
         y2 = list(goal = "min", low = 0, high = 3)),
    c("y1", "y2")
  )
  expect_type(result, "list")
  expect_equal(names(result), c("y1", "y2"))
})