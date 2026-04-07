# Tests for coding policy enforcement in fitting functions

test_that("coding policy helper warns/errors/ignores as expected", {
  dat <- data.frame(
    x1 = c(-2, -1, 0, 1, 2),
    x2 = c(10, 20, 30, 40, 50),
    y = c(1, 2, 3, 4, 5)
  )

  expect_warning(
    .brsm_enforce_coding_policy(
      data = dat,
      factor_names = c("x1", "x2"),
      coding = NULL,
      coding_policy = "warn",
      caller = "fit_brsm"
    ),
    "no brsm coding metadata"
  )

  expect_error(
    .brsm_enforce_coding_policy(
      data = dat,
      factor_names = c("x1", "x2"),
      coding = NULL,
      coding_policy = "error",
      caller = "fit_brsm"
    ),
    "no brsm coding metadata"
  )

  expect_no_warning(
    .brsm_enforce_coding_policy(
      data = dat,
      factor_names = c("x1", "x2"),
      coding = NULL,
      coding_policy = "ignore",
      caller = "fit_brsm"
    )
  )
})

test_that("coding policy helper is silent when coding metadata exists", {
  dat <- data.frame(x1 = 1:5, x2 = 6:10, y = 11:15)
  coding <- list(
    method = "zscore",
    factors = list(
      x1 = list(center = 3, scale = 1, method = "zscore"),
      x2 = list(center = 8, scale = 1, method = "zscore")
    )
  )

  expect_no_warning(
    .brsm_enforce_coding_policy(
      data = dat,
      factor_names = c("x1", "x2"),
      coding = coding,
      coding_policy = "warn",
      caller = "fit_brsm"
    )
  )
})
