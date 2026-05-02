# Tests for prepare_brsm_data, get_brsm_coding, decode_brsm_data

test_that("prepare_brsm_data with zscore method centers and scales correctly", {
  df <- data.frame(
    x1 = c(1, 2, 3, 4, 5),
    x2 = c(10, 20, 30, 40, 50),
    y = c(5, 10, 15, 20, 25)
  )

  result <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "zscore"
  )

  # Check center (should be original means)
  expect_equal(mean(df$x1), 3)
  expect_equal(mean(df$x2), 30)

  # Check transformed data has mean near 0 and SD near 1
  expect_true(abs(mean(result$x1)) < 1e-10)
  expect_true(abs(mean(result$x2)) < 1e-10)
  expect_true(abs(sd(result$x1) - 1) < 1e-10)
  expect_true(abs(sd(result$x2) - 1) < 1e-10)

  # y should be unchanged
  expect_identical(result$y, df$y)
})

test_that("prepare_brsm_data with range method scales to [-1,1]", {
  df <- data.frame(
    x1 = c(1, 2, 3, 4, 5),
    x2 = c(10, 20, 30, 40, 50),
    y = c(5, 10, 15, 20, 25)
  )

  result <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "range"
  )

  # Check all scaled values in [-1, 1]
  expect_true(all(result$x1 >= -1 & result$x1 <= 1))
  expect_true(all(result$x2 >= -1 & result$x2 <= 1))

  # Check min/max are exactly -1 and 1
  expect_equal(min(result$x1), -1)
  expect_equal(max(result$x1), 1)
  expect_equal(min(result$x2), -1)
  expect_equal(max(result$x2), 1)

  # y should be unchanged
  expect_identical(result$y, df$y)
})

test_that(
  "prepare_brsm_data with identity method preserves pre-coded factors",
  {
  df <- data.frame(
    x1 = c(-1, -0.5, 0, 0.5, 1),
    x2 = c(-0.8, -0.2, 0, 0.3, 0.9),
    y = c(5, 10, 15, 20, 25)
  )

  result <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "identity"
  )

  expect_identical(result$x1, df$x1)
  expect_identical(result$x2, df$x2)

  coding <- get_brsm_coding(result)
  expect_equal(coding$method, "identity")
  expect_equal(coding$factors$x1$center, 0)
  expect_equal(coding$factors$x1$scale, 1)
  expect_equal(coding$factors$x2$center, 0)
  expect_equal(coding$factors$x2$scale, 1)
})

test_that("prepare_brsm_data stores coding metadata in attributes", {
  df <- data.frame(x1 = c(1, 2, 3), x2 = c(10, 20, 30), y = c(5, 10, 15))

  result <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "zscore"
  )

  # Check brsm_coded_data attribute
  expect_true(inherits(result, "brsm_coded_data"))

  # Check brsm_coding attribute exists
  coding <- attr(result, "brsm_coding", exact = TRUE)
  expect_false(is.null(coding))

  # Check structure: list with method and factors
  expect_true(all(c("method", "factors") %in% names(coding)))
  expect_equal(coding$method, "zscore")
  expect_named(coding$factors, c("x1", "x2"))
  expect_true(all(c("center", "scale", "method") %in% names(coding$factors$x1)))
  expect_true(all(c("center", "scale", "method") %in% names(coding$factors$x2)))
})

test_that("get_brsm_coding extracts metadata from prepared data", {
  df <- data.frame(x1 = c(1, 2, 3), x2 = c(10, 20, 30), y = c(5, 10, 15))
  result <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "range"
  )

  coding <- get_brsm_coding(result)

  expect_false(is.null(coding))
  expect_equal(coding$method, "range")
  expect_length(coding$factors, 2)
})

test_that("get_brsm_coding errors when no coding present", {
  df <- data.frame(x1 = c(1, 2, 3), y = c(5, 10, 15))

  expect_error(get_brsm_coding(df), "No brsm coding metadata found")
})

test_that("decode_brsm_data reverses zscore transformation correctly", {
  df <- data.frame(
    x1 = c(1, 2, 3, 4, 5),
    x2 = c(10, 20, 30, 40, 50),
    y = c(5, 10, 15, 20, 25)
  )

  prepared <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "zscore"
  )

  # Create a subset of prepared data for decoding (simulating predictions)
  pred_data <- prepared[, c("x1", "x2")]

  decoded <- decode_brsm_data(pred_data, coding = get_brsm_coding(prepared))

  # Check that decoding reverses the transformation
  expect_true(all(abs(decoded$x1 - df$x1) < 1e-10))
  expect_true(all(abs(decoded$x2 - df$x2) < 1e-10))
})

test_that("decode_brsm_data reverses range transformation correctly", {
  df <- data.frame(
    x1 = c(1, 2, 3, 4, 5),
    x2 = c(10, 20, 30, 40, 50),
    y = c(5, 10, 15, 20, 25)
  )

  prepared <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "range"
  )

  # Create a subset of prepared data for decoding
  pred_data <- prepared[, c("x1", "x2")]

  decoded <- decode_brsm_data(pred_data, coding = get_brsm_coding(prepared))

  # Check that decoding reverses the transformation
  expect_true(all(abs(decoded$x1 - df$x1) < 1e-10))
  expect_true(all(abs(decoded$x2 - df$x2) < 1e-10))
})

test_that("decode_brsm_data is a no-op for identity coding", {
  df <- data.frame(
    x1 = c(-1, -0.5, 0, 0.5, 1),
    x2 = c(-0.8, -0.2, 0, 0.3, 0.9)
  )

  prepared <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "identity"
  )

  decoded <- decode_brsm_data(
    prepared[, c("x1", "x2")],
    coding = get_brsm_coding(prepared)
  )

  expect_identical(decoded$x1, df$x1)
  expect_identical(decoded$x2, df$x2)
}
)

test_that("decode_brsm_data with suffix appends correctly", {
  df <- data.frame(x1 = c(1, 2, 3), x2 = c(10, 20, 30))
  prepared <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "zscore"
  )

  pred_data <- prepared[, c("x1", "x2")]

  decoded <- decode_brsm_data(
    pred_data,
    coding = get_brsm_coding(prepared),
    overwrite = FALSE,
    suffix = "_decoded"
  )

  expect_named(decoded, c("x1", "x2", "x1_decoded", "x2_decoded"))
})

test_that("decode_brsm_data works with brsm_fit object coding", {
  skip_if_no_brms_tests()

  df <- generate_simulation_data(n = 25, seed = 701)
  prepared <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "zscore"
  )

  fit <- fit_brsm(
    data = prepared,
    response = "y",
    factor_names = c("x1", "x2"),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 701,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  decoded <- decode_brsm_data(
    prepared[, c("x1", "x2")],
    coding = get_brsm_coding(fit)
  )
  expect_true(all(abs(decoded$x1 - df$x1) < 1e-8))
  expect_true(all(abs(decoded$x2 - df$x2) < 1e-8))
})

test_that("prepare_brsm_data handles single factor", {
  df <- data.frame(x1 = c(1, 2, 3, 4, 5), y = c(5, 10, 15, 20, 25))

  result <- prepare_brsm_data(df, factor_names = "x1", method = "zscore")

  expect_true(inherits(result, "brsm_coded_data"))
  expect_true(abs(mean(result$x1)) < 1e-10)
  expect_true(abs(sd(result$x1) - 1) < 1e-10)
})

test_that("prepare_brsm_data handles edge case: zero standard deviation", {
  df <- data.frame(x1 = c(5, 5, 5, 5), y = c(1, 2, 3, 4))

  expect_error(
    prepare_brsm_data(df, factor_names = "x1", method = "zscore"),
    "zero/invalid scale"
  )
})

test_that("prepare_brsm_data handles edge case: single value after range", {
  df <- data.frame(x1 = c(5, 5, 5, 5), y = c(1, 2, 3, 4))

  expect_error(
    prepare_brsm_data(df, factor_names = "x1", method = "range"),
    "zero/invalid scale"
  )
})

test_that("prepare_brsm_data preserves non-factor columns", {
  df <- data.frame(
    x1 = c(1, 2, 3),
    x2 = c(10, 20, 30),
    y = c(5, 10, 15),
    group = c("A", "B", "C")
  )

  result <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "zscore"
  )

  expect_identical(result$y, df$y)
  expect_identical(result$group, df$group)
})

test_that("decode_brsm_data with missing factor columns errors gracefully", {
  df <- data.frame(x1 = c(1, 2, 3), x2 = c(10, 20, 30), y = c(5, 10, 15))
  prepared <- prepare_brsm_data(
    df,
    factor_names = c("x1", "x2"),
    method = "zscore"
  )

  # Try to decode with missing x2
  pred_data <- prepared[, "x1", drop = FALSE]

  expect_error(decode_brsm_data(pred_data, prepared))
})

test_that("coded_data_helpers cover remaining validation and metadata branches", {
  expect_error(
    prepare_brsm_data(1, factor_names = c("x1", "x2")),
    "data must be a data.frame"
  )

  df <- data.frame(x1 = c(1, 2, 3), y = c(2, 3, 4))
  expect_error(
    prepare_brsm_data(df, factor_names = c("x1", "x2")),
    "Missing factor columns"
  )

  bad_num <- data.frame(x1 = c("a", "b"), x2 = c(1, 2), y = c(1, 2))
  expect_error(
    prepare_brsm_data(bad_num, factor_names = c("x1", "x2")),
    "All factor columns must be numeric"
  )

  all_na <- data.frame(x1 = c(NA_real_, NA_real_), x2 = c(1, 2), y = c(1, 2))
  expect_error(
    prepare_brsm_data(all_na, factor_names = c("x1", "x2")),
    "contains only NA values"
  )

  fake_fit <- structure(
    list(coding = list(method = "zscore", factors = list(x1 = list(center = 0, scale = 1)))),
    class = "brsm_fit"
  )
  coding_fit <- get_brsm_coding(fake_fit)
  expect_equal(coding_fit$method, "zscore")

  bad_coding_attr <- structure(
    data.frame(x1 = 1:2),
    brsm_coding = list(method = "zscore")
  )
  expect_error(
    get_brsm_coding(bad_coding_attr),
    "Invalid brsm coding metadata format"
  )

  expect_error(
    decode_brsm_data(1, coding = list(factors = list(x1 = list(center = 0, scale = 1)))),
    "data must be a data.frame"
  )

  plain <- data.frame(x1 = c(0, 1), x2 = c(0, 1))
  expect_error(
    decode_brsm_data(plain),
    "coding is NULL and no 'brsm_coding' attribute found"
  )

  expect_error(
    decode_brsm_data(plain, coding = list(method = "zscore")),
    "coding must be a list with a 'factors' list entry"
  )

  coding <- list(
    method = "zscore",
    factors = list(
      x1 = list(center = 1, scale = 2),
      x2 = list(center = 2, scale = 3)
    )
  )
  expect_error(
    decode_brsm_data(plain, coding = coding, factor_names = c("x1", "x3")),
    "No coding metadata found for factors"
  )

  expect_error(
    decode_brsm_data(plain[, "x1", drop = FALSE], coding = coding, factor_names = c("x1", "x2")),
    "data is missing coded factor column"
  )
})
