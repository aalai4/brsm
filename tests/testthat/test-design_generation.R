# Tests for generate_brsm_design

test_that("generate_brsm_design creates CCD with expected size", {
  out <- generate_brsm_design(
    factor_names = c("x1", "x2"),
    design = "ccd",
    n_center = 3,
    randomize = FALSE
  )

  # 2^k factorial + 2k axial + n_center
  expect_equal(nrow(out), 4 + 4 + 3)
  expect_true(all(c("x1", "x2", ".design_type", ".block", ".run_id") %in% names(out)))
  expect_true(all(out$.design_type == "ccd"))
})

test_that("generate_brsm_design creates BBD with expected size", {
  out <- generate_brsm_design(
    factor_names = c("x1", "x2", "x3"),
    design = "bbd",
    n_center = 2,
    randomize = FALSE
  )

  # 4 * choose(k, 2) + n_center
  expect_equal(nrow(out), 4 * choose(3, 2) + 2)
  expect_true(all(out$.design_type == "bbd"))
})

test_that("generate_brsm_design randomization can be seeded", {
  a <- generate_brsm_design(
    factor_names = c("x1", "x2"),
    design = "ccd",
    randomize = TRUE,
    seed = 101
  )
  b <- generate_brsm_design(
    factor_names = c("x1", "x2"),
    design = "ccd",
    randomize = TRUE,
    seed = 101
  )

  expect_equal(a$x1, b$x1)
  expect_equal(a$x2, b$x2)
  expect_equal(a$.block, b$.block)
})

test_that("generate_brsm_design supports natural output with ranges", {
  out <- generate_brsm_design(
    factor_names = c("x1", "x2"),
    design = "ccd",
    n_center = 1,
    randomize = FALSE,
    output = "natural",
    ranges = list(x1 = c(10, 20), x2 = c(100, 200))
  )

  expect_true(min(out$x1) < 10)
  expect_true(max(out$x1) > 20)
  expect_true(min(out$x2) < 100)
  expect_true(max(out$x2) > 200)
})

test_that("generate_brsm_design coded output can carry coding metadata", {
  out <- generate_brsm_design(
    factor_names = c("x1", "x2"),
    design = "ccd",
    randomize = FALSE,
    output = "coded",
    ranges = list(x1 = c(0, 10), x2 = c(-5, 5))
  )

  expect_true(inherits(out, "brsm_coded_data"))
  coding <- get_brsm_coding(out)
  expect_equal(coding$method, "range")
  expect_true(all(c("x1", "x2") %in% names(coding$factors)))
})

test_that("generate_brsm_design validates BBD factor count", {
  expect_error(
    generate_brsm_design(
      factor_names = c("x1", "x2"),
      design = "bbd"
    ),
    "BBD requires at least 3 factors"
  )
})

test_that("generate_brsm_design validates alpha", {
  expect_error(
    generate_brsm_design(
      factor_names = c("x1", "x2"),
      design = "ccd",
      alpha = "bad"
    ),
    "alpha must be"
  )
})

test_that("generate_brsm_design validates ranges for natural output", {
  expect_error(
    generate_brsm_design(
      factor_names = c("x1", "x2"),
      design = "ccd",
      output = "natural"
    ),
    "ranges must be provided"
  )

  expect_error(
    generate_brsm_design(
      factor_names = c("x1", "x2"),
      design = "ccd",
      output = "natural",
      ranges = list(x1 = c(0, 1))
    ),
    "include all factor_names"
  )
})