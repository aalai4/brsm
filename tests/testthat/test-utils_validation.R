test_that("internal validators cover edge branches", {
  # .brsm_validate_factor_names
  expect_error(
    brsm:::.brsm_validate_factor_names(
      c("x1"),
      require_two = TRUE,
      two_factors_message = "need exactly two"
    ),
    "need exactly two"
  )
  expect_error(
    brsm:::.brsm_validate_factor_names(c("x1", "x2", "x3"), require_two = TRUE),
    "exactly two factors"
  )

  # .brsm_validate_draws
  expect_error(
    brsm:::.brsm_validate_draws(data.frame()[0, , drop = FALSE], empty_message = "empty draws"),
    "empty draws"
  )

  # .brsm_validate_probs finite message path
  expect_error(
    brsm:::.brsm_validate_probs(c(0.1, NA_real_), require_finite = TRUE),
    "finite values between 0 and 1"
  )
})

test_that(".brsm_validate_bayesian_input handles all accepted and rejected classes", {
  expect_error(
    brsm:::.brsm_validate_bayesian_input(structure(list(fit = NULL), class = "brsm_fit")),
    "valid brmsfit"
  )

  fake_brmsfit <- structure(list(), class = "brmsfit")
  valid_brsm_fit <- structure(list(fit = fake_brmsfit), class = "brsm_fit")
  expect_invisible(brsm:::.brsm_validate_bayesian_input(valid_brsm_fit))
  expect_invisible(brsm:::.brsm_validate_bayesian_input(fake_brmsfit))

  expect_error(
    brsm:::.brsm_validate_bayesian_input(data.frame(a = 1:3)),
    "lacks Bayesian posterior columns"
  )

  expect_error(
    brsm:::.brsm_validate_bayesian_input(list(a = 1)),
    "must be strictly Bayesian"
  )
})

test_that("coding policy and helper checks cover warning and stop branches", {
  dat <- data.frame(x1 = c(-1, 0, 1), x2 = c(-1, 0, 1))
  expect_warning(
    brsm:::.brsm_enforce_coding_policy(
      data = dat,
      factor_names = c("x1", "x2"),
      coding = NULL,
      coding_policy = "warn",
      caller = "fit_brsm"
    ),
    "approximately centered/scaled already"
  )

  expect_error(
    brsm:::.brsm_check_columns(c("x1", "x3"), dat, "missing columns"),
    "missing columns"
  )

  expect_true(is.na(brsm:::.brsm_find_quadratic_col("x9", names(dat))))

  # Cover regex-fallback match branch in .brsm_find_quadratic_col.
  odd_cols <- c("b_Intercept", "b_x1", "b_Ix1\\.2\\.")
  expect_equal(
    brsm:::.brsm_find_quadratic_col("x1", odd_cols),
    "b_Ix1\\.2\\."
  )
})

test_that(".brsm_hessian_array validates numeric quadratic and interaction coefficients", {
  bad_quad <- data.frame(
    b_x1 = c(1, 1),
    b_x2 = c(1, 1),
    `b_I(x1^2)` = c("bad", "bad"),
    `b_I(x2^2)` = c(-1, -1),
    `b_x1:x2` = c(0, 0),
    check.names = FALSE
  )
  expect_error(
    brsm:::.brsm_hessian_array(bad_quad, c("x1", "x2")),
    "numeric columns for quadratic terms"
  )

  bad_inter <- data.frame(
    b_x1 = c(1, 1),
    b_x2 = c(1, 1),
    `b_I(x1^2)` = c(-1, -1),
    `b_I(x2^2)` = c(-1, -1),
    `b_x1:x2` = c("bad", "bad"),
    check.names = FALSE
  )
  expect_error(
    brsm:::.brsm_hessian_array(bad_inter, c("x1", "x2")),
    "numeric columns for interaction terms"
  )
})
