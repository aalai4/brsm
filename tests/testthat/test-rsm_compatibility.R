# Tests for rsm coded.data compatibility

test_that(".brsm_parse_rsm_formula parses standard formulas correctly", {
  # Standard subtraction + division
  f1 <- x1 ~ (Temp - 150) / 10
  p1 <- brsm:::.brsm_parse_rsm_formula(f1)
  expect_equal(p1$coded, "x1")
  expect_equal(p1$original, "Temp")
  expect_equal(p1$center, 150)
  expect_equal(p1$scale, 10)

  # Addition + division
  f2 <- x2 ~ (Concentration + 0.5) / 0.1
  p2 <- brsm:::.brsm_parse_rsm_formula(f2)
  expect_equal(p2$coded, "x2")
  expect_equal(p2$original, "Concentration")
  expect_equal(p2$center, -0.5)
  expect_equal(p2$scale, 0.1)

  # Simple subtraction without division
  f3 <- x3 ~ Time - 5
  p3 <- brsm:::.brsm_parse_rsm_formula(f3)
  expect_equal(p3$coded, "x3")
  expect_equal(p3$original, "Time")
  expect_equal(p3$center, 5)
  expect_equal(p3$scale, 1)

  # Identity expression
  f4 <- x4 ~ Speed
  p4 <- brsm:::.brsm_parse_rsm_formula(f4)
  expect_equal(p4$coded, "x4")
  expect_equal(p4$original, "Speed")
  expect_equal(p4$center, 0)
  expect_equal(p4$scale, 1)
})

test_that("prepare_brsm_data imports coded.data objects and decode_brsm_data decodes them", {
  # Build a mock coded.data frame
  mock_coded <- data.frame(
    x1 = c(-1, 0, 1),
    x2 = c(-1, 0, 1)
  )
  attr(mock_coded, "codings") <- list(
    x1 = x1 ~ (Temp - 150) / 10,
    x2 = x2 ~ (Concentration - 0.5) / 0.1
  )
  class(mock_coded) <- c("coded.data", class(mock_coded))

  # Ingestion without specifying factor names
  coded_brsm <- prepare_brsm_data(mock_coded)
  expect_s3_class(coded_brsm, "brsm_coded_data")
  expect_s3_class(coded_brsm, "coded.data")

  coding <- get_brsm_coding(coded_brsm)
  expect_equal(coding$method, "rsm")
  expect_equal(names(coding$factors), c("x1", "x2"))

  expect_equal(coding$factors$x1$center, 150)
  expect_equal(coding$factors$x1$scale, 10)
  expect_equal(coding$factors$x1$original_range, c(140, 160))

  expect_equal(coding$factors$x2$center, 0.5)
  expect_equal(coding$factors$x2$scale, 0.1)
  expect_equal(coding$factors$x2$original_range, c(0.4, 0.6))

  # Decoding verification
  decoded <- decode_brsm_data(coded_brsm)
  expect_equal(decoded$x1, c(140, 150, 160))
  expect_equal(decoded$x2, c(0.4, 0.5, 0.6))
  expect_false(inherits(decoded, "brsm_coded_data"))
})
