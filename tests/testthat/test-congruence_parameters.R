# Tests for congruence_parameters, summarize_congruence, rope_congruence

.create_mock_draws <- function(n_draws = 100, seed = 123) {
  set.seed(seed)
  data.frame(
    b_Intercept = rnorm(n_draws, 50, 5),
    b_x1 = rnorm(n_draws, 0.5, 0.2),
    b_x2 = rnorm(n_draws, 0.5, 0.2),
    "b_I(x1^2)" = rnorm(n_draws, -0.3, 0.1),
    "b_I(x2^2)" = rnorm(n_draws, -0.25, 0.1),
    "b_x1:x2" = rnorm(n_draws, 0.1, 0.05),
    check.names = FALSE
  )
}

test_that("congruence_parameters computes a1-a5 correctly", {
  draws <- .create_mock_draws(n_draws = 10)

  result <- congruence_parameters(draws, factor_names = c("x1", "x2"))

  expect_data_frame(result, ncols = 5, nrows = 10)
  expect_named(result, c("a1", "a2", "a3", "a4", "a5"))

  expect_equal(result$a1[1], draws$b_x1[1] + draws$b_x2[1])
  expect_equal(result$a2[1], draws[["b_I(x1^2)"]][1] + draws[["b_I(x2^2)"]][1])
  expect_equal(result$a3[1], draws[["b_x1:x2"]][1])
  expect_equal(
    result$a4[1],
    -2 * draws[["b_I(x1^2)"]][1] * draws[["b_I(x2^2)"]][1]
  )
  expect_equal(
    result$a5[1],
    abs(draws[["b_I(x1^2)"]][1] - draws[["b_I(x2^2)"]][1])
  )
})

test_that("congruence_parameters preserves row count", {
  for (n in c(1, 10, 100)) {
    draws <- .create_mock_draws(n_draws = n)
    result <- congruence_parameters(draws, factor_names = c("x1", "x2"))
    expect_equal(nrow(result), n)
  }
})

test_that("summarize_congruence computes statistics and quantiles", {
  draws <- .create_mock_draws(n_draws = 200)
  summary_df <- summarize_congruence(
    draws,
    factor_names = c("x1", "x2"),
    probs = c(0.025, 0.5, 0.975)
  )

  expect_data_frame(summary_df, nrows = 5)
  expect_named(
    summary_df,
    c("parameter", "mean", "median", "q_002", "q_050", "q_097")
  )
  expect_equal(summary_df$parameter, c("a1", "a2", "a3", "a4", "a5"))
})

test_that("rope_congruence classifies practical equivalence when all inside", {
  draws <- data.frame(
    b_Intercept = rep(0, 50),
    b_x1 = rep(0, 50),
    b_x2 = rep(0, 50),
    "b_I(x1^2)" = rep(0, 50),
    "b_I(x2^2)" = rep(0, 50),
    "b_x1:x2" = rep(0, 50),
    check.names = FALSE
  )

  result <- rope_congruence(
    draws,
    factor_names = c("x1", "x2"),
    rope = c(-0.1, 0.1)
  )

  expect_data_frame(result, nrows = 5)
  expect_named(
    result,
    c("parameter", "inside_rope", "outside_rope", "decision")
  )
  expect_true(all(result$inside_rope == 1))
  expect_true(all(result$outside_rope == 0))
  expect_true(all(result$decision == "Practical Equivalence"))
})

test_that("rope_congruence classifies not equivalent when all outside", {
  draws <- data.frame(
    b_Intercept = rep(0, 50),
    b_x1 = rep(2, 50),
    b_x2 = rep(2, 50),
    "b_I(x1^2)" = rep(3, 50),
    "b_I(x2^2)" = rep(1, 50),
    "b_x1:x2" = rep(2, 50),
    check.names = FALSE
  )

  result <- rope_congruence(
    draws,
    factor_names = c("x1", "x2"),
    rope = c(-0.1, 0.1)
  )

  expect_true(all(result$inside_rope == 0))
  expect_true(all(result$outside_rope == 1))
  expect_true(all(result$decision == "Not Equivalent"))
})

test_that("rope_congruence can be inconclusive", {
  # a1 will be 50% inside and 50% outside ROPE
  draws <- data.frame(
    b_Intercept = rep(0, 100),
    b_x1 = c(rep(0.02, 50), rep(0.8, 50)),
    b_x2 = rep(0, 100),
    "b_I(x1^2)" = rep(0, 100),
    "b_I(x2^2)" = rep(0, 100),
    "b_x1:x2" = rep(0, 100),
    check.names = FALSE
  )

  result <- rope_congruence(
    draws,
    factor_names = c("x1", "x2"),
    rope = c(-0.1, 0.1)
  )
  a1_row <- result[result$parameter == "a1", , drop = FALSE]

  expect_equal(a1_row$inside_rope, 0.5)
  expect_equal(a1_row$outside_rope, 0.5)
  expect_equal(a1_row$decision, "Inconclusive")
})

test_that("rope_congruence validates rope bounds", {
  draws <- .create_mock_draws(n_draws = 20)
  expect_error(
    rope_congruence(draws, factor_names = c("x1", "x2"), rope = c(1, -1)),
    "must be less than"
  )
})

test_that(
  "congruence_parameters handles sanitized quadratic and interaction names",
  {
  set.seed(404)
  draws <- data.frame(
    b_Intercept = rnorm(30),
    b_x1 = rnorm(30),
    b_x2 = rnorm(30),
    b_Ix1E2 = rnorm(30, -0.2, 0.05),
    b_Ix2E2 = rnorm(30, -0.25, 0.05),
    b_x1.x2 = rnorm(30, 0.1, 0.02),
    check.names = FALSE
  )

  out <- congruence_parameters(draws, factor_names = c("x1", "x2"))

  expect_data_frame(out, ncols = 5, nrows = 30)
  expect_equal(out$a2[1], draws$b_Ix1E2[1] + draws$b_Ix2E2[1])
  expect_equal(out$a3[1], draws$b_x1.x2[1])
}
)

test_that("as_brsm_draws canonicalizes sanitized quadratic naming", {
  set.seed(405)
  draws_raw <- data.frame(
    b_Intercept = rnorm(20),
    b_x1 = rnorm(20),
    b_x2 = rnorm(20),
    b_Ix1E2 = rnorm(20, -0.2, 0.05),
    b_Ix2E2 = rnorm(20, -0.25, 0.05),
    b_x1.x2 = rnorm(20, 0.1, 0.02),
    check.names = FALSE
  )

  out <- as_brsm_draws(draws_raw, factor_names = c("x1", "x2"))

  expect_true("b_I(x1^2)" %in% names(out))
  expect_true("b_I(x2^2)" %in% names(out))
  expect_true("b_x1:x2" %in% names(out))
  expect_false("b_Ix1E2" %in% names(out))
  expect_false("b_Ix2E2" %in% names(out))
  expect_false("b_x1.x2" %in% names(out))
})
