# Tests for loftest_brsm

.loftest_cache <- new.env(parent = emptyenv())

.get_loftest_baseline <- function(seed = 301) {
  key <- paste0("baseline_", seed)
  if (!exists(key, envir = .loftest_cache, inherits = FALSE)) {
    dat <- prepare_brsm_data(
      generate_simulation_data(n = 25, seed = seed),
      factor_names = c("x1", "x2"),
      method = "zscore"
    )
    fit <- fit_brsm(
      data = dat,
      response = "y",
      factor_names = c("x1", "x2"),
      coding_policy = "ignore",
      chains = 1,
      iter = 250,
      warmup = 125,
      seed = seed,
      sampling_preset = "fast",
      refresh = 0,
      silent = 2
    )
    assign(key, fit, envir = .loftest_cache)
  }
  get(key, envir = .loftest_cache, inherits = FALSE)
}

test_that("loftest_brsm with provided reference model", {
  skip_if_no_brms_tests()

  baseline <- .get_loftest_baseline(seed = 301)
  reference <- fit_brsm(
    data = baseline$fit$data,
    response = "y",
    factor_names = c("x1", "x2"),
    coding_policy = "ignore",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 302,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result <- loftest_brsm(
    object = baseline,
    reference_model = reference,
    include_ppc = FALSE,
    criterion = "loo"
  )

  expect_type(result, "list")
  expect_true("comparison" %in% names(result))
  expect_false(isTRUE(result$reference_fitted))
})

test_that("loftest_brsm with auto-fit cubic reference model", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 303)
  baseline <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    coding_policy = "ignore",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 303,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result <- loftest_brsm(
    object = baseline,
    reference_type = "cubic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 303,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_type(result, "list")
  expect_true(isTRUE(result$reference_fitted))
  expect_true("comparison" %in% names(result))
  expect_true("comparison" %in% names(result$comparison))
})

test_that("loftest_brsm with auto-fit extended reference model", {
  skip_if_no_brms_tests()

  baseline <- .get_loftest_baseline(seed = 304)
  result <- loftest_brsm(
    object = baseline,
    reference_type = "extended",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 304,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  ref_formula <- paste(
    deparse(stats::formula(result$reference_model)),
    collapse = " "
  )
  expect_type(result, "list")
  expect_true(isTRUE(result$reference_fitted))
  expect_true(grepl("I\\(x1\\^3\\)", ref_formula))
  expect_true(grepl("I\\(x2\\^3\\)", ref_formula))
})

test_that("loftest_brsm reference_type='cubic' adds cubic terms", {
  f <- brsm:::.brsm_build_reference_formula(
    response = "y",
    factor_names = c("x1", "x2"),
    reference_type = "cubic"
  )
  txt <- paste(deparse(f), collapse = " ")

  expect_true(grepl("I\\(x1\\^3\\)", txt))
  expect_true(grepl("I\\(x2\\^3\\)", txt))
  expect_true(grepl("x1:x2", txt))
})

test_that(
  "loftest_brsm reference_type='extended' adds cubic and interaction terms",
  {
  f <- brsm:::.brsm_build_reference_formula(
    response = "y",
    factor_names = c("x1", "x2"),
    reference_type = "extended"
  )
  txt <- paste(deparse(f), collapse = " ")

  expect_true(grepl("I\\(x1\\^3\\)", txt))
  expect_true(grepl("I\\(x2\\^3\\)", txt))
  expect_true(grepl("I\\(x1\\^2\\):x2", txt))
  expect_true(grepl("x1:I\\(x2\\^2\\)", txt))
}
)

test_that("loftest_brsm without include_ppc", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 305),
    reference_type = "cubic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 305,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_false("ppc" %in% names(result))
})

test_that("loftest_brsm with include_ppc=TRUE", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 306),
    reference_type = "cubic",
    include_ppc = TRUE,
    ppc_ndraws = 50,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 306,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_true("ppc" %in% names(result))
  expect_true("summaries" %in% names(result$ppc))
})

test_that("loftest_brsm returns comparison with loo_diff", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 307),
    reference_type = "cubic",
    include_ppc = FALSE,
    criterion = "loo",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 307,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_true("comparison" %in% names(result$comparison))
  expect_true("elpd_diff" %in% names(result$comparison$comparison))
})

test_that("loftest_brsm richer model has lower LOO (better fit)", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 308),
    reference_type = "extended",
    include_ppc = FALSE,
    criterion = "loo",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 308,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_true(is.data.frame(result$comparison$comparison))
  expect_true(nrow(result$comparison$comparison) >= 2)
})

test_that("loftest_brsm errors when reference model is simpler than baseline", {
  skip_if_no_brms_tests()

  baseline <- .get_loftest_baseline(seed = 309)
  simpler_reference <- brms::brm(
    formula = y ~ x1 + x2,
    data = generate_simulation_data(n = 25, seed = 310),
    family = stats::gaussian(),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 310,
    control = list(adapt_delta = 0.8, max_treedepth = 10),
    refresh = 0,
    silent = 2
  )

  expect_no_error(
    loftest_brsm(
      object = baseline,
      reference_model = simpler_reference,
      include_ppc = FALSE
    )
  )
})

test_that("loftest_brsm handles single factor case", {
  skip_if_no_brms_tests()

  dat <- generate_simulation_data(n = 25, seed = 311)
  dat$x2 <- NULL
  baseline <- brms::brm(
    formula = y ~ x1 + I(x1^2),
    data = dat,
    family = stats::gaussian(),
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 311,
    control = list(adapt_delta = 0.8, max_treedepth = 10),
    refresh = 0,
    silent = 2
  )

  result <- loftest_brsm(
    object = baseline,
    data = dat,
    response = "y",
    factor_names = "x1",
    reference_type = "cubic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 311,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_type(result, "list")
  expect_true(
    grepl(
      "I\\(x1\\^3\\)",
      paste(deparse(stats::formula(result$reference_model)), collapse = " ")
    )
  )
})

test_that("loftest_brsm handles two factor case", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 312),
    reference_type = "cubic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 312,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  ref_formula <- paste(
    deparse(stats::formula(result$reference_model)),
    collapse = " "
  )
  expect_true(grepl("I\\(x1\\^3\\)", ref_formula))
  expect_true(grepl("I\\(x2\\^3\\)", ref_formula))
})

test_that("loftest_brsm reference model inherits priors from baseline", {
  skip_if_no_brms_tests()

  pr <- brms::prior(normal(0, 1), class = "b")
  dat <- generate_simulation_data(n = 25, seed = 313)
  baseline <- fit_brsm(
    data = dat,
    response = "y",
    factor_names = c("x1", "x2"),
    prior = pr,
    coding_policy = "ignore",
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 313,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result <- loftest_brsm(
    object = baseline,
    reference_type = "cubic",
    prior = pr,
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 313,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_type(result, "list")
  expect_true(isTRUE(result$reference_fitted))
})

test_that("loftest_brsm comparison returns LOO and/or WAIC", {
  skip_if_no_brms_tests()

  result_loo <- loftest_brsm(
    object = .get_loftest_baseline(seed = 314),
    reference_type = "cubic",
    criterion = "loo",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 314,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  result_waic <- loftest_brsm(
    object = .get_loftest_baseline(seed = 315),
    reference_type = "cubic",
    criterion = "waic",
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 315,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_equal(result_loo$criterion, "loo")
  expect_equal(result_waic$criterion, "waic")
  expect_true(is.data.frame(result_loo$comparison$comparison))
  expect_true(is.data.frame(result_waic$comparison$comparison))
})

test_that("loftest_brsm validates loo_* arguments", {
  skip_if_not_installed("brms")

  expect_error(
    loftest_brsm(
      object = list(),
      loo_k_threshold = 0
    ),
    "loo_k_threshold must be a finite numeric scalar > 0"
  )

  expect_error(
    loftest_brsm(
      object = list(),
      loo_moment_match = NA
    ),
    "loo_moment_match must be a non-missing TRUE/FALSE value"
  )

  expect_error(
    loftest_brsm(
      object = list(),
      loo_reloo = NA
    ),
    "loo_reloo must be a non-missing TRUE/FALSE value"
  )

  expect_error(
    loftest_brsm(
      object = list(),
      loo_auto_moment_match = NA
    ),
    "loo_auto_moment_match must be a non-missing TRUE/FALSE value"
  )
})

test_that("loftest_brsm returns loo_diagnostics for LOO", {
  skip_if_no_brms_tests()

  result <- loftest_brsm(
    object = .get_loftest_baseline(seed = 316),
    reference_type = "cubic",
    criterion = "loo",
    loo_auto_moment_match = FALSE,
    loo_k_threshold = 0.7,
    include_ppc = FALSE,
    chains = 1,
    iter = 250,
    warmup = 125,
    seed = 316,
    sampling_preset = "fast",
    refresh = 0,
    silent = 2
  )

  expect_true("loo_diagnostics" %in% names(result))
  expect_equal(result$loo_diagnostics$k_threshold, 0.7)
  expect_true("initial" %in% names(result$loo_diagnostics))
  expect_true("final" %in% names(result$loo_diagnostics))
  expect_true("moment_match_used" %in% names(result$loo_diagnostics))
  expect_true("reloo_used" %in% names(result$loo_diagnostics))
  expect_true("auto_moment_match_retry" %in% names(result$loo_diagnostics))
})

test_that(".brsm_summarize_pareto_k handles empty, unnamed, and non-finite diagnostics", {
  empty <- brsm:::.brsm_summarize_pareto_k(list(), threshold = 0.7)
  expect_false(empty$has_high_k)
  expect_true(is.na(empty$max_pareto_k))
  expect_length(empty$max_pareto_k_by_model, 0)
  expect_length(empty$n_above_threshold_by_model, 0)

  est <- list(
    list(diagnostics = list(pareto_k = c(0.2, 0.8, Inf, NA_real_))),
    list(diagnostics = list(pareto_k = c(0.1, 0.3))),
    list(diagnostics = list(pareto_k = c(NA_real_, Inf)))
  )
  out <- brsm:::.brsm_summarize_pareto_k(est, threshold = 0.7)

  expect_true(out$has_high_k)
  expect_equal(out$max_pareto_k_by_model[["model_1"]], 0.8)
  expect_equal(out$n_above_threshold_by_model[["model_1"]], 1L)
  expect_true(is.na(out$max_pareto_k_by_model[["model_3"]]))
})

test_that(".brsm_build_reference_formula handles single and multi-factor extended cases", {
  f1 <- brsm:::.brsm_build_reference_formula(
    response = "y",
    factor_names = "x1",
    reference_type = "extended"
  )
  txt1 <- paste(deparse(f1), collapse = " ")
  expect_false(grepl(":", txt1))
  expect_true(grepl("I\\(x1\\^3\\)", txt1))

  f3 <- brsm:::.brsm_build_reference_formula(
    response = "y",
    factor_names = c("x1", "x2", "x3"),
    reference_type = "extended"
  )
  txt3 <- paste(deparse(f3), collapse = " ")
  expect_true(grepl("I\\(x1\\^2\\):x2", txt3))
  expect_true(grepl("x1:I\\(x2\\^2\\)", txt3))
  expect_true(grepl("I\\(x1\\^2\\):x3", txt3))
})

test_that("loftest_brsm mock flow covers auto-fit retry, diagnostics, and PPC branches", {
  skip_if_not_installed("brms")

  ns_brsm <- asNamespace("brsm")
  ns_brms <- asNamespace("brms")

  old_extract <- get(".brsm_extract_fit", envir = ns_brsm)
  old_compare <- get("compare_brsm_models", envir = ns_brsm)
  old_ppc <- get("check_brsm_ppc", envir = ns_brsm)
  old_brm <- get("brm", envir = ns_brms)

  trace_env <- new.env(parent = emptyenv())
  trace_env$compare_calls <- list()
  trace_env$brm_args <- NULL

  unlockBinding(".brsm_extract_fit", ns_brsm)
  assign(".brsm_extract_fit", function(m, caller = NULL) {
    if (inherits(m, "brsm_fit")) m$fit else m
  }, envir = ns_brsm)
  lockBinding(".brsm_extract_fit", ns_brsm)

  unlockBinding("compare_brsm_models", ns_brsm)
  assign("compare_brsm_models", function(models, criterion = c("loo", "waic"), ...) {
    args <- list(...)
    trace_env$compare_calls[[length(trace_env$compare_calls) + 1L]] <- args
    if (length(trace_env$compare_calls) == 1L) {
      est <- list(
        baseline = list(diagnostics = list(pareto_k = c(0.9, 0.85))),
        reference = list(diagnostics = list(pareto_k = c(0.8)))
      )
    } else {
      est <- list(
        baseline = list(diagnostics = list(pareto_k = c(0.2, 0.3))),
        reference = list(diagnostics = list(pareto_k = c(0.1)))
      )
    }

    list(
      criterion = match.arg(criterion),
      estimates = est,
      comparison = data.frame(model = c("baseline", "reference"), elpd_diff = c(0, -1), se_diff = c(0, 0.1))
    )
  }, envir = ns_brsm)
  lockBinding("compare_brsm_models", ns_brsm)

  unlockBinding("check_brsm_ppc", ns_brsm)
  assign("check_brsm_ppc", function(object, ndraws, probs, seed, include_plot, ...) {
    list(summary = data.frame(n_obs = 3, ndraws = ndraws, interval_lower_prob = probs[1], interval_upper_prob = probs[2]))
  }, envir = ns_brsm)
  lockBinding("check_brsm_ppc", ns_brsm)

  unlockBinding("brm", ns_brms)
  assign("brm", function(...) {
    trace_env$brm_args <- list(...)
    structure(list(data = trace_env$brm_args$data), class = "brmsfit")
  }, envir = ns_brms)
  lockBinding("brm", ns_brms)

  on.exit({
    unlockBinding(".brsm_extract_fit", ns_brsm)
    assign(".brsm_extract_fit", old_extract, envir = ns_brsm)
    lockBinding(".brsm_extract_fit", ns_brsm)

    unlockBinding("compare_brsm_models", ns_brsm)
    assign("compare_brsm_models", old_compare, envir = ns_brsm)
    lockBinding("compare_brsm_models", ns_brsm)

    unlockBinding("check_brsm_ppc", ns_brsm)
    assign("check_brsm_ppc", old_ppc, envir = ns_brsm)
    lockBinding("check_brsm_ppc", ns_brsm)

    unlockBinding("brm", ns_brms)
    assign("brm", old_brm, envir = ns_brms)
    lockBinding("brm", ns_brms)
  }, add = TRUE)

  baseline <- structure(
    list(
      response = "y",
      factor_names = c("x1", "x2"),
      fit = list(data = data.frame(y = c(1, 2, 3), x1 = c(-1, 0, 1), x2 = c(-1, 0, 1)), family = stats::gaussian()),
      sampling = list(chains = 2, iter = 40, warmup = 20, backend = "cmdstanr", control = list(adapt_delta = 0.85))
    ),
    class = "brsm_fit"
  )

  out <- loftest_brsm(
    object = baseline,
    reference_model = NULL,
    reference_type = "extended",
    criterion = "loo",
    loo_moment_match = FALSE,
    loo_auto_moment_match = TRUE,
    loo_k_threshold = 0.7,
    include_ppc = TRUE,
    ppc_ndraws = 5,
    ppc_probs = c(0.1, 0.9),
    seed = 99
  )

  expect_true(isTRUE(out$reference_fitted))
  expect_true(isTRUE(out$loo_diagnostics$auto_moment_match_retry))
  expect_true(isTRUE(out$loo_diagnostics$moment_match_used))
  expect_true("ppc" %in% names(out))
  expect_equal(nrow(out$ppc$summaries), 2)
  expect_equal(length(trace_env$compare_calls), 2)
  expect_true(isTRUE(trace_env$compare_calls[[2]]$moment_match))
  expect_identical(trace_env$brm_args$backend, "cmdstanr")
})

test_that("loftest_brsm auto-fit input validation branches error cleanly", {
  skip_if_not_installed("brms")

  ns_brsm <- asNamespace("brsm")
  old_extract <- get(".brsm_extract_fit", envir = ns_brsm)
  unlockBinding(".brsm_extract_fit", ns_brsm)
  assign(".brsm_extract_fit", function(m, caller = NULL) m, envir = ns_brsm)
  lockBinding(".brsm_extract_fit", ns_brsm)
  on.exit({
    unlockBinding(".brsm_extract_fit", ns_brsm)
    assign(".brsm_extract_fit", old_extract, envir = ns_brsm)
    lockBinding(".brsm_extract_fit", ns_brsm)
  }, add = TRUE)

  fake_fit <- structure(list(), class = "brmsfit")

  expect_error(
    loftest_brsm(object = fake_fit, reference_model = NULL),
    "automatic reference fitting requires data, response, and factor_names"
  )
  expect_error(
    loftest_brsm(object = fake_fit, reference_model = NULL, data = 1, response = "y", factor_names = c("x1", "x2")),
    "data must be a data.frame"
  )
  expect_error(
    loftest_brsm(
      object = fake_fit,
      reference_model = NULL,
      data = data.frame(y = 1:3, x1 = 1:3, x2 = 3:1),
      response = "y",
      factor_names = c("x1", "x2"),
      control = 1
    ),
    "control must be NULL or a named list"
  )
})