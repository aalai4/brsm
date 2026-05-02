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

test_that("prior term resolver maps sanitized quadratic and interaction names", {
  term_groups <- list(
    linear = c("x1", "x2"),
    interaction = c("x1:x2"),
    quadratic = c("I(x1^2)", "I(x2^2)")
  )

  available <- c("x1", "x2", "x1.x2", "Ix1E2", "Ix2E2")
  resolved <- .brsm_resolve_b_prior_targets(term_groups, available)

  expect_equal(resolved$matched$linear, c("x1", "x2"))
  expect_equal(resolved$matched$interaction, "x1.x2")
  expect_equal(resolved$matched$quadratic, c("Ix1E2", "Ix2E2"))
  expect_length(unlist(resolved$unmatched), 0)
})

test_that("prior term resolver records unmatched terms for class-b fallback", {
  term_groups <- list(
    linear = c("x1", "x2"),
    interaction = c("x1:x2"),
    quadratic = c("I(x1^2)", "I(x2^2)")
  )

  available <- c("x1", "x2", "x1:x2", "I(x1^2)")
  resolved <- .brsm_resolve_b_prior_targets(term_groups, available)

  expect_equal(resolved$matched$quadratic, "I(x1^2)")
  expect_equal(resolved$unmatched$quadratic, "I(x2^2)")
})

test_that("specify_brsm_priors validates additional flags and autoscale response", {
  skip_if_not_installed("brms")

  dat_non_numeric <- data.frame(x1 = 1:4, x2 = 2:5, y = letters[1:4])
  expect_error(
    specify_brsm_priors(
      factor_names = c("x1", "x2"),
      autoscale = TRUE,
      data = dat_non_numeric,
      response = "y"
    ),
    "response column must be numeric"
  )

  dat_zero_sd <- data.frame(x1 = 1:4, x2 = 2:5, y = c(1, 1, 1, 1))
  expect_error(
    specify_brsm_priors(
      factor_names = c("x1", "x2"),
      autoscale = TRUE,
      data = dat_zero_sd,
      response = "y"
    ),
    "response standard deviation must be finite"
  )

  expect_error(
    specify_brsm_priors(
      factor_names = c("x1", "x2"),
      include_intercept = NA
    ),
    "include_intercept must be TRUE or FALSE"
  )

  expect_error(
    specify_brsm_priors(
      factor_names = c("x1", "x2"),
      autoscale = NA
    ),
    "autoscale must be TRUE or FALSE"
  )
})

test_that("adaptive profile shrinks priors and can add global b fallback", {
  skip_if_not_installed("brms")

  tiny_dat <- data.frame(x1 = 0, x2 = 0, y = 1)
  p_adapt <- as.data.frame(specify_brsm_priors(
    factor_names = c("x1", "x2"),
    model_terms = "second_order",
    prior_profile = "adaptive",
    data = tiny_dat,
    response = "y"
  ))
  x1_row <- p_adapt[p_adapt$class == "b" & p_adapt$coef == "x1", , drop = FALSE]
  expect_equal(x1_row$prior[[1]], "normal(0, 0.8)")

  ns <- asNamespace("brsm")
  unlockBinding(".brsm_default_prior_b_coefs", ns)
  old_default <- get(".brsm_default_prior_b_coefs", envir = ns)
  assign(".brsm_default_prior_b_coefs", function(...) c("x1"), envir = ns)
  lockBinding(".brsm_default_prior_b_coefs", ns)
  on.exit({
    unlockBinding(".brsm_default_prior_b_coefs", ns)
    assign(".brsm_default_prior_b_coefs", old_default, envir = ns)
    lockBinding(".brsm_default_prior_b_coefs", ns)
  }, add = TRUE)

  p_fallback <- as.data.frame(specify_brsm_priors(
    factor_names = c("x1", "x2"),
    model_terms = "first_order",
    data = data.frame(x1 = 1:4, x2 = 2:5, y = 3:6),
    response = "y"
  ))
  expect_true(any(p_fallback$class == "b" & p_fallback$coef == ""))
})

test_that("prior helper edge branches return expected NA or fallback names", {
  expect_equal(
    .brsm_resolve_linear_coef("x 1", c("x.1")),
    "x.1"
  )

  expect_true(is.na(.brsm_resolve_interaction_coef("x1", c("x1:x2"))))

  expect_true(is.na(.brsm_resolve_quadratic_coef("I(x9^2)", c("x1", "x2"))))

  out <- .brsm_default_prior_b_coefs(
    response = "y",
    linear_terms = character(0),
    interaction_terms = character(0),
    quadratic_terms = character(0),
    include_interactions = FALSE,
    include_quadratic = FALSE,
    data = data.frame(y = 1:3)
  )
  expect_null(out)

  term_groups <- list(
    linear = c("x1", "x2"),
    interaction = c("x1:x2"),
    quadratic = character(0)
  )
  resolved <- .brsm_resolve_b_prior_targets(term_groups, c("x1"))
  expect_equal(resolved$unmatched$linear, "x2")
  expect_equal(resolved$unmatched$interaction, "x1:x2")
})