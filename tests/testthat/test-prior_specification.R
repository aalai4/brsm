# Tests for prior specification utilities

test_that("specify_brsm_priors returns brmsprior object", {
  skip_if_not_installed("brms")

  priors <- specify_brsm_priors(
    factor_names = c("x1", "x2"),
    model_terms = "second_order"
  )

  expect_s3_class(priors, "brmsprior")
  p <- as.data.frame(priors)
  expect_true(all(c("prior", "class", "coef") %in% names(p)))
})

test_that("specify_brsm_priors includes expected terms by model_terms", {
  skip_if_not_installed("brms")

  p_second <- as.data.frame(specify_brsm_priors(
    factor_names = c("x1", "x2", "x3"),
    model_terms = "second_order"
  ))

  expect_true("Intercept" %in% p_second$class)
  expect_true("sigma" %in% p_second$class)
  expect_true("x1" %in% p_second$coef)
  expect_true("x1:x2" %in% p_second$coef)
  expect_true("I(x1^2)" %in% p_second$coef)

  p_first <- as.data.frame(specify_brsm_priors(
    factor_names = c("x1", "x2"),
    model_terms = "first_order"
  ))

  expect_false(any(grepl(":", p_first$coef, fixed = TRUE), na.rm = TRUE))
  expect_false(any(grepl("I\\(", p_first$coef), na.rm = TRUE))

  p_twi <- as.data.frame(specify_brsm_priors(
    factor_names = c("x1", "x2"),
    model_terms = "first_order_twi"
  ))

  expect_true("x1:x2" %in% p_twi$coef)
  expect_false(any(grepl("I\\(", p_twi$coef), na.rm = TRUE))

  p_quad <- as.data.frame(specify_brsm_priors(
    factor_names = c("x1", "x2"),
    model_terms = "pure_quadratic"
  ))

  expect_false(any(grepl(":", p_quad$coef, fixed = TRUE), na.rm = TRUE))
  expect_true("I(x1^2)" %in% p_quad$coef)
})

test_that("specify_brsm_priors supports student_t coefficient family", {
  skip_if_not_installed("brms")

  p <- as.data.frame(specify_brsm_priors(
    factor_names = c("x1", "x2"),
    coefficient_family = "student_t",
    student_df = 5
  ))

  slope_rows <- p[p$class %in% c("b", "Intercept"), , drop = FALSE]
  expect_true(all(grepl("student_t\\(", slope_rows$prior)))
})

test_that("specify_brsm_priors autoscale uses response sd", {
  skip_if_not_installed("brms")

  dat <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50),
    y = rnorm(50, sd = 4)
  )

  p <- as.data.frame(specify_brsm_priors(
    factor_names = c("x1", "x2"),
    linear_sd = 2,
    autoscale = TRUE,
    data = dat,
    response = "y"
  ))

  linear_row <- p[p$class == "b" & p$coef == "x1", , drop = FALSE]
  expect_true(grepl("normal\\(0,", linear_row$prior))
  expect_false(grepl("normal\\(0, 2\\)", linear_row$prior))
})

test_that("specify_brsm_priors validates autoscale inputs", {
  skip_if_not_installed("brms")

  dat <- data.frame(x1 = 1:5, x2 = 1:5, y = 1:5)

  expect_error(
    specify_brsm_priors(
      factor_names = c("x1", "x2"),
      autoscale = TRUE
    ),
    "data must be provided"
  )

  expect_error(
    specify_brsm_priors(
      factor_names = c("x1", "x2"),
      autoscale = TRUE,
      data = dat
    ),
    "response must be a single character"
  )

  expect_error(
    specify_brsm_priors(
      factor_names = c("x1", "x2"),
      autoscale = TRUE,
      data = dat,
      response = "missing"
    ),
    "response column not found"
  )
})

test_that("specify_brsm_priors validates numeric scales and flags", {
  skip_if_not_installed("brms")

  expect_error(
    specify_brsm_priors(
      factor_names = c("x1", "x2"),
      linear_sd = 0
    ),
    "linear_sd must be a finite numeric scalar > 0"
  )

  expect_error(
    specify_brsm_priors(
      factor_names = c("x1", "x2"),
      include_sigma = NA
    ),
    "include_sigma must be TRUE or FALSE"
  )
})