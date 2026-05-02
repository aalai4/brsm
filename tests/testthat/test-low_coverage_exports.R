# Focused tests for previously uncovered exported helpers.

.make_local_mock_brsm_fit <- function(include_coding = TRUE) {
  obj <- list(
    formula = y ~ x1 + x2 + I(x1^2) + I(x2^2) + x1:x2,
    response = "y",
    factor_names = c("x1", "x2"),
    coding = if (include_coding) list(method = "coded") else NULL,
    ranges = list(x1 = c(-1, 1), x2 = c(-2, 2)),
    sampling = list(chains = 4, iter = 2000, warmup = 1000,
      sampling_preset = "standard"
    ),
    fit = NULL,
    call = quote(fit_brsm(dat, response = "y", factor_names = c("x1", "x2")))
  )
  class(obj) <- "brsm_fit"
  obj
}

.create_fake_brmsfit <- function(data = data.frame(y = c(1, 2, 4)), summary_obj = NULL) {
  structure(
    list(
      data = data,
      .summary = if (is.null(summary_obj)) list(fixed = data.frame()) else summary_obj,
      score = sum(data$y)
    ),
    class = c("fake_brmsfit", "brmsfit")
  )
}

summary.fake_brmsfit <- function(object, ...) {
  object$.summary
}

registerS3method(
  "summary",
  "fake_brmsfit",
  summary.fake_brmsfit,
  envir = asNamespace("base")
)

get_y.fake_brmsfit <- function(object, ...) {
  object$data$y
}
pp_check.fake_brmsfit <- function(object, ndraws, type, ...) {
  list(ndraws = ndraws, type = type)
}
posterior_predict.fake_brmsfit <- function(object, ndraws, ...) {
  if (isTRUE(object$.bad_pp)) {
    return(1)
  }
  matrix(c(1, 2, 3, 2, 3, 4), nrow = ndraws, byrow = TRUE)
}
loo.fake_brmsfit <- function(object, ...) {
  ll <- matrix(log(c(0.6, 0.7, 0.8, 0.65, 0.75, 0.85)), nrow = 2, byrow = TRUE)
  loo::loo(ll)
}
waic.fake_brmsfit <- function(object, ...) {
  ll <- matrix(log(c(0.6, 0.7, 0.8, 0.65, 0.75, 0.85)), nrow = 2, byrow = TRUE)
  loo::waic(ll)
}
fixef.fake_brmsfit <- function(object, ...) {
  matrix(c(1, 2), ncol = 1, dimnames = list(c("b_x1", "b_x2"), "Estimate"))
}
rhat.fake_brmsfit <- function(object, ...) {
  c(b_x1 = 1.00, b_x2 = 1.01)
}

if (requireNamespace("brms", quietly = TRUE)) {
  for (.mn in c("get_y", "posterior_predict", "pp_check", "loo", "waic",
                "fixef", "rhat")) {
    .fn <- get(paste0(.mn, ".fake_brmsfit"))
    registerS3method(.mn, "fake_brmsfit", .fn, envir = asNamespace("brms"))
  }
}

test_that("check_brsm_ppc validates inputs before execution", {
  skip_if_not_installed("brms")

  fit <- .create_fake_brmsfit(data = data.frame(y = c(1, 2, 4)))

  expect_error(brsm::check_brsm_ppc(fit, ndraws = 0), "ndraws")
  expect_error(brsm::check_brsm_ppc(fit, ndraws = -1), "ndraws")
  expect_error(brsm::check_brsm_ppc(fit, probs = c(-0.1, 0.9)), "between 0 and 1")
  expect_error(brsm::check_brsm_ppc(fit, probs = c(0.1, 0.5, 0.9)), "between 0 and 1|length 2")
  expect_error(brsm::check_brsm_ppc(list(a = 1)), "brsm_fit or brmsfit")
})

test_that("check_brsm_ppc returns summary, optional plot, and matrix validation", {
  skip_if_not_installed("brms")

  fit <- .create_fake_brmsfit(data = data.frame(y = c(1, 2, 4)))

  ns <- asNamespace("brms")
  old_get_y <- get("get_y", envir = ns)

  unlockBinding("get_y", ns)
  assign("get_y", function(object, ...) object$data$y, envir = ns)
  lockBinding("get_y", ns)

  on.exit({
    unlockBinding("get_y", ns)
    assign("get_y", old_get_y, envir = ns)
    lockBinding("get_y", ns)
  }, add = TRUE)

  out <- brsm::check_brsm_ppc(
    fit,
    ndraws = 2,
    probs = c(0.1, 0.9),
    seed = 123,
    include_plot = TRUE
  )

  expect_true(all(c("summary", "observed", "predicted_mean", "plot") %in% names(out)))
  expect_s3_class(out$summary, "data.frame")
  expect_equal(out$summary$n_obs, 3)
  expect_equal(out$summary$ndraws, 2)
  expect_equal(out$summary$interval_lower_prob, 0.1)
  expect_equal(out$summary$interval_upper_prob, 0.9)
  expect_equal(out$plot$type, "hist")

  bad_pp <- .create_fake_brmsfit(data = data.frame(y = c(1, 2, 4)))
  bad_pp$.bad_pp <- TRUE
  expect_error(
    brsm::check_brsm_ppc(bad_pp),
    "posterior_predict did not return a matrix-like object"
  )
})

.replace_brms_binding <- function(name, new_fn) {
  ns <- asNamespace("brms")
  old <- get(name, envir = ns)
  unlockBinding(name, ns)
  assign(name, new_fn, envir = ns)
  lockBinding(name, ns)
  old
}

.restore_brms_binding <- function(name, old_fn) {
  ns <- asNamespace("brms")
  unlockBinding(name, ns)
  assign(name, old_fn, envir = ns)
  lockBinding(name, ns)
}

test_that("compare_brsm_models validates input before dispatch", {
  skip_if_not_installed("brms")

  fit1 <- .create_fake_brmsfit(data = data.frame(y = c(1, 2, 3)))
  fit2 <- .create_fake_brmsfit(data = data.frame(y = c(2, 3, 4)))

  expect_error(
    brsm::compare_brsm_models(list(fit1), criterion = "loo"),
    "at least two"
  )

})

test_that("compare_brsm_models handles unnamed and named model lists", {
  skip_if_not_installed("brms")
  skip_if_not_installed("loo")

  fit1 <- .create_fake_brmsfit(data = data.frame(y = c(1, 2, 3)))
  fit2 <- .create_fake_brmsfit(data = data.frame(y = c(2, 3, 4)))

  res1 <- suppressWarnings(
    brsm::compare_brsm_models(list(fit1, fit2), criterion = "loo")
  )
  expect_identical(res1$criterion, "loo")
  expect_equal(nrow(res1$comparison), 2)

  fit1_brsm <- structure(list(fit = fit1), class = "brsm_fit")
  res2 <- suppressWarnings(
    brsm::compare_brsm_models(list(a = fit1_brsm, b = fit2), criterion = "waic")
  )
  expect_identical(res2$criterion, "waic")
  expect_equal(nrow(res2$comparison), 2)
})

test_that("compare_brsm_models returns structured output for loo and waic", {
  skip_if_not_installed("brms")
  skip_if_not_installed("loo")

  fit1 <- .create_fake_brmsfit(data = data.frame(y = c(1, 2, 3)))
  fit2 <- .create_fake_brmsfit(data = data.frame(y = c(2, 3, 5)))

  res_loo <- suppressWarnings(
    brsm::compare_brsm_models(list(fit1, fit2), criterion = "loo")
  )
  expect_identical(res_loo$criterion, "loo")
  expect_true(all(c("model", "elpd_diff", "se_diff") %in% names(res_loo$comparison)))
  expect_equal(nrow(res_loo$comparison), 2)

  fit1_brsm <- structure(list(fit = fit1), class = "brsm_fit")
  res_waic <- suppressWarnings(
    brsm::compare_brsm_models(list(a = fit1_brsm, b = fit2), criterion = "waic")
  )
  expect_identical(res_waic$criterion, "waic")
  expect_equal(nrow(res_waic$comparison), 2)

  expect_error(
    brsm::compare_brsm_models(list(a = fit1, b = fit2), criterion = "bad"),
    "'arg' should be one of"
  )
})

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
  skip_if_not_installed("brms")

  dat <- data.frame(x1 = c(-1, 0, 1), x2 = c(-2, 0, 2), y = c(1, 2, 4))

  expect_error(brsm::fit_brsm(1, "y", c("x1", "x2")), "data.frame")
  expect_error(brsm::fit_brsm(dat, c("y", "z"), c("x1", "x2")), "single character")
  expect_error(brsm::fit_brsm(dat, "y", c("x1", "x2"), control = 1), "control")

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

test_that("fit_brsm builds brsm_fit with mocked brm backend", {
  skip_if_not_installed("brms")

  dat <- data.frame(x1 = c(-1, 0, 1), x2 = c(-2, 0, 2), y = c(1, 2, 4))
  attr(dat, "brsm_coding") <- list(method = "zscore")

  fake_brm <- function(formula, data, prior, family, chains, iter, warmup, seed, control, ...) {
    extra <- list(...)
    structure(
      list(
        data = data, formula = formula, prior = prior, family = family,
        chains = chains, iter = iter, warmup = warmup, seed = seed,
        control = control, extra = extra
      ),
      class = "brmsfit"
    )
  }

  old_brm <- .replace_brms_binding("brm", fake_brm)
  on.exit(.restore_brms_binding("brm", old_brm), add = TRUE)

  fit <- brsm::fit_brsm(
    data = dat, response = "y", factor_names = c("x1", "x2"),
    prior_profile = "regularized", sampling_preset = "robust",
    model_terms = "pure_quadratic", control = list(adapt_delta = 0.95),
    backend = "cmdstanr", chains = 2, iter = 20, warmup = 10, seed = 1
  )

  expect_s3_class(fit, "brsm_fit")
  expect_identical(format(fit$formula), "y ~ x1 + x2 + I(x1^2) + I(x2^2)")
  expect_equal(fit$ranges$x1, c(-1, 1))
  expect_identical(fit$coding$method, "zscore")
  expect_equal(fit$sampling$control$adapt_delta, 0.95)
  expect_equal(fit$sampling$control$max_treedepth, 15)
  expect_identical(fit$sampling$backend, "cmdstanr")

  .restore_brms_binding("brm", old_brm)

  fake_brm2 <- function(formula, data, prior, family, chains, iter, warmup, seed, control, ...) {
    structure(list(data = data, control = control), class = "brmsfit")
  }
  old_brm2 <- .replace_brms_binding("brm", fake_brm2)
  on.exit(.restore_brms_binding("brm", old_brm2), add = TRUE)

  expect_warning(
    brsm::fit_brsm(
      data = data.frame(x1 = c(-1, 0, 1), x2 = c(-2, 0, 2), y = c(1, 2, 4)),
      response = "y", factor_names = c("x1", "x2"),
      prior = brms::prior("normal(0, 1)", class = "b"),
      coding_policy = "warn", chains = 1, iter = 10, warmup = 5
    ),
    "fit_brsm.*called with no brsm coding metadata"
  )
})

test_that("print and summary methods cover fit diagnostics branches", {
  skip_if_not_installed("brms")

  mock <- .make_local_mock_brsm_fit()
  mock$fit <- structure(
    list(data = data.frame(y = c(1, 2, 3))),
    class = c("fake_brmsfit", "brmsfit")
  )

  printed <- capture.output(brsm:::print.brsm_fit(mock))
  expect_true(any(grepl("Observations", printed)))
  expect_true(any(grepl("Parameters estimated", printed)))
  expect_true(any(grepl("Max Rhat", printed)))

  summary_obj <- list(fixed = data.frame(Estimate = 1, row.names = "b_x1"))
  mock2 <- .make_local_mock_brsm_fit()
  mock2$fit <- .create_fake_brmsfit(summary_obj = summary_obj)

  summarized <- brsm:::summary.brsm_fit(mock2)
  expect_s3_class(summarized, "summary.brsm_fit")
  expect_identical(summarized$response, "y")

  printed_summary <- capture.output(brsm:::print.summary.brsm_fit(summarized))
  expect_true(any(grepl("Bayesian Response Surface Model Summary", printed_summary)))
  expect_true(any(grepl("Coefficient Summary", printed_summary)))
})

# ---------------------------------------------------------------------------
# brms-absent guard branches
# ---------------------------------------------------------------------------

# Helper: patch base::requireNamespace so that `pkg` always returns FALSE,
# then restore on exit.  Returns the old function for use in on.exit.
.patch_requireNamespace_for <- function(pkg) {
  base_ns <- asNamespace("base")
  old_rns <- get("requireNamespace", envir = base_ns)
  unlockBinding("requireNamespace", base_ns)
  assign(
    "requireNamespace",
    function(package, quietly = TRUE) if (identical(package, pkg)) FALSE else old_rns(package, quietly = quietly),
    envir = base_ns
  )
  lockBinding("requireNamespace", base_ns)
  invisible(list(ns = base_ns, old = old_rns))
}

.restore_requireNamespace <- function(info) {
  unlockBinding("requireNamespace", info$ns)
  assign("requireNamespace", info$old, envir = info$ns)
  lockBinding("requireNamespace", info$ns)
}

test_that("check_brsm_ppc errors when brms is unavailable", {
  info <- .patch_requireNamespace_for("brms")
  on.exit(.restore_requireNamespace(info), add = TRUE)

  fake <- structure(list(), class = c("cbf_fake_brmsfit", "brmsfit"))
  expect_error(
    brsm::check_brsm_ppc(fake),
    "package 'brms' is required"
  )
})

test_that("fit_brsm errors when brms is unavailable", {
  info <- .patch_requireNamespace_for("brms")
  on.exit(.restore_requireNamespace(info), add = TRUE)

  dat <- data.frame(x1 = c(-1, 0, 1), x2 = c(-1, 0, 1), y = c(1, 2, 3))
  expect_error(
    brsm::fit_brsm(dat, "y", c("x1", "x2")),
    "package 'brms' is required"
  )
})
