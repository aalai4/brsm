# Focused tests for previously uncovered exported helpers.

test_that(".brsm_require_ggplot2 errors when ggplot2 is unavailable", {
  base_ns <- asNamespace("base")
  old_req <- get("requireNamespace", envir = base_ns)

  unlockBinding("requireNamespace", base_ns)
  assign("requireNamespace", function(package, quietly = TRUE) FALSE, envir = base_ns)
  lockBinding("requireNamespace", base_ns)

  on.exit({
    unlockBinding("requireNamespace", base_ns)
    assign("requireNamespace", old_req, envir = base_ns)
    lockBinding("requireNamespace", base_ns)
  }, add = TRUE)

  expect_error(
    brsm:::.brsm_require_ggplot2(),
    "package 'ggplot2' is required"
  )
})

test_that("fit_brsm validates inputs before fitting", {
  dat <- data.frame(x1 = c(-1, 0, 1), x2 = c(-2, 0, 2), y = c(1, 2, 4))

  expect_error(brsm::fit_brsm(1, "y", c("x1", "x2")), "data.frame")
  expect_error(brsm::fit_brsm(dat, c("y", "z"), c("x1", "x2")), "single character")

  expect_error(
    suppressWarnings(brsm::fit_brsm(dat, "missing", c("x1", "x2"))),
    "Missing required columns"
  )

  bad_factor <- dat
  bad_factor$x1 <- letters[1:3]
  expect_error(
    suppressWarnings(brsm::fit_brsm(bad_factor, "y", c("x1", "x2"))),
    "All factor columns must be numeric"
  )

  bad_response <- dat
  bad_response$y <- letters[1:3]
  expect_error(
    suppressWarnings(brsm::fit_brsm(bad_response, "y", c("x1", "x2"))),
    "response column must be numeric"
  )
})

test_that("check_brsm_ppc validates inputs before execution", {
  fake_fit <- list(
    X = matrix(1:6, ncol = 2),
    V_n = diag(2),
    y = rnorm(3),
    draws = data.frame(sigma = rnorm(3))
  )
  class(fake_fit) <- "brsm_conjugate_fit"
  fit <- structure(list(fit = fake_fit), class = "brsm_fit")

  expect_error(brsm::check_brsm_ppc(fit, ndraws = 0), "ndraws")
  expect_error(brsm::check_brsm_ppc(fit, ndraws = -1), "ndraws")
  expect_error(brsm::check_brsm_ppc(fit, probs = c(-0.1, 0.9)), "between 0 and 1")
  expect_error(brsm::check_brsm_ppc(fit, probs = c(0.1, 0.5, 0.9)), "between 0 and 1|length 2")
  expect_error(brsm::check_brsm_ppc(list(a = 1)), "brsm_fit or brsm_conjugate_fit")
})
