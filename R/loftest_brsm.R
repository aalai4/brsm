#' Bayesian Lack-of-Fit Convenience Wrapper for brsm
#'
#' Runs a Bayesian lack-of-fit style check by comparing a baseline model against
#' a richer reference model using LOO/WAIC, with optional side-by-side posterior
#' predictive summaries.
#'
#' The function can either:
#' 1. Accept an already fitted richer model via \code{reference_model}, or
#' 2. Fit a richer reference model automatically (for \code{brsm_fit} baseline
#'    objects, or when fitting inputs are supplied).
#'
#' @param object Baseline model as a \code{brsm_fit} or \code{brmsfit} object.
#' @param reference_model Optional richer reference model as \code{brsm_fit} or
#'   \code{brmsfit}. If \code{NULL}, the function attempts to fit one.
#' @param reference_type Type of richer reference formula used when
#'   \code{reference_model = NULL}. One of \code{"cubic"} (adds cubic terms)
#'   or \code{"extended"} (adds cubic and quadratic-by-linear terms).
#' @param criterion Information criterion for model comparison: \code{"loo"}
#'   (default) or \code{"waic"}.
#' @param loo_moment_match Logical; when \code{criterion = "loo"}, pass
#'   \code{moment_match} to [brms::loo()] via [compare_brsm_models()].
#' @param loo_reloo Logical; when \code{criterion = "loo"}, pass
#'   \code{reloo} to [brms::loo()] via [compare_brsm_models()].
#' @param loo_k_threshold Pareto-k threshold used to flag unstable LOO
#'   diagnostics when \code{criterion = "loo"}. Defaults to \code{0.7}.
#' @param loo_auto_moment_match Logical; when \code{TRUE} and
#'   \code{criterion = "loo"}, automatically re-run model comparison with
#'   \code{moment_match = TRUE} if any Pareto-k values exceed
#'   \code{loo_k_threshold}.
#' @param include_ppc Logical; if \code{TRUE}, compute
#'   [check_brsm_ppc()] for both baseline and reference models.
#' @param ppc_ndraws Number of posterior predictive draws for PPC.
#' @param ppc_probs Length-2 numeric vector of predictive interval probabilities
#'   for PPC.
#' @param seed Optional random seed for reference fitting and PPC.
#' @param data Optional data frame used to fit the reference model when
#'   \code{reference_model = NULL} and \code{object} is not \code{brsm_fit}.
#' @param response Optional response variable name for automatic reference
#'   fitting. Inferred from \code{object} when possible.
#' @param factor_names Optional factor names for automatic reference fitting.
#'   Inferred from \code{object} when possible.
#' @param prior Optional prior passed to \code{brms::brm()} when fitting the
#'   reference model.
#' @param family Optional family passed to \code{brms::brm()} when fitting the
#'   reference model.
#' @param chains Number of chains for automatic reference fitting.
#' @param iter Number of iterations per chain for automatic reference fitting.
#' @param warmup Number of warmup iterations per chain for automatic reference
#'   fitting.
#' @param sampling_preset Sampling profile for auto-fitted reference model.
#'   One of \code{"fast"}, \code{"balanced"}, or \code{"robust"}.
#' @param backend Optional backend passed to \code{brms::brm()}.
#' @param control Optional named list of NUTS control arguments.
#' @param ... Additional arguments passed to \code{brms::brm()} when fitting
#'   the reference model.
#'
#' @return A list with components:
#'   \code{criterion}, \code{comparison} (from [compare_brsm_models()]),
#'   \code{baseline_model}, \code{reference_model}, \code{reference_type},
#'   \code{reference_fitted} (logical), optional \code{loo_diagnostics} (for
#'   \code{criterion = "loo"}), and optional \code{ppc} list with side-by-side
#'   summaries.
#'
#' @examples
#' \dontrun{
#' fit_q <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
#'
#' # Automatic cubic reference fit + LOO comparison
#' lof <- loftest_brsm(
#'   object = fit_q,
#'   reference_type = "cubic",
#'   include_ppc = TRUE
#' )
#' lof$comparison$comparison
#' lof$ppc$summaries
#' }
#' @export
loftest_brsm <- function(object,
                         reference_model = NULL,
                         reference_type = c("cubic", "extended"),
                         criterion = c("loo", "waic"),
                         loo_moment_match = FALSE,
                         loo_reloo = FALSE,
                         loo_k_threshold = 0.7,
                         loo_auto_moment_match = TRUE,
                         include_ppc = FALSE,
                         ppc_ndraws = 200,
                         ppc_probs = c(0.025, 0.975),
                         seed = NULL,
                         data = NULL,
                         response = NULL,
                         factor_names = NULL,
                         prior = NULL,
                         family = NULL,
                         chains = 4,
                         iter = 2000,
                         warmup = floor(iter / 2),
                         sampling_preset = c("fast", "balanced", "robust"),
                         backend = NULL,
                         control = NULL,
                         ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("package 'brms' is required for loftest_brsm().")
  }

  reference_type <- match.arg(reference_type)
  criterion <- match.arg(criterion)
  sampling_preset <- match.arg(sampling_preset)

  if (!is.logical(loo_moment_match) || length(loo_moment_match) != 1L ||
      is.na(loo_moment_match)) {
    stop("loo_moment_match must be a non-missing TRUE/FALSE value.")
  }
  if (!is.logical(loo_reloo) || length(loo_reloo) != 1L ||
      is.na(loo_reloo)) {
    stop("loo_reloo must be a non-missing TRUE/FALSE value.")
  }
  if (!is.numeric(loo_k_threshold) || length(loo_k_threshold) != 1L ||
      !is.finite(loo_k_threshold) || loo_k_threshold <= 0) {
    stop("loo_k_threshold must be a finite numeric scalar > 0.")
  }
  if (!is.logical(loo_auto_moment_match) ||
      length(loo_auto_moment_match) != 1L ||
      is.na(loo_auto_moment_match)) {
    stop("loo_auto_moment_match must be a non-missing TRUE/FALSE value.")
  }

  baseline_fit <- .brsm_extract_fit(object, caller = "loftest_brsm")
  reference_fitted <- FALSE

  if (is.null(reference_model)) {
    # Pull defaults from brsm_fit metadata when available
    if (inherits(object, "brsm_fit")) {
      if (is.null(data)) {
        data <- baseline_fit$data
      }
      if (is.null(response)) {
        response <- object$response
      }
      if (is.null(factor_names)) {
        factor_names <- object$factor_names
      }

      if (is.null(family) && !is.null(object$fit$family)) {
        family <- object$fit$family
      }

      if (!is.null(object$sampling)) {
        if (missing(chains)) {
          chains <- object$sampling$chains
        }
        if (missing(iter)) {
          iter <- object$sampling$iter
        }
        if (missing(warmup)) {
          warmup <- object$sampling$warmup
        }
        if (missing(backend) && !is.null(object$sampling$backend)) {
          backend <- object$sampling$backend
        }
        if (is.null(control) && !is.null(object$sampling$control)) {
          control <- object$sampling$control
        }
      }
    }

    if (is.null(data) || is.null(response) || is.null(factor_names)) {
      stop(
        "When reference_model is NULL, automatic reference",
        " fitting requires data, ",
        "response, and factor_names (or a brsm_fit object containing them)."
      )
    }

    if (!is.data.frame(data)) {
      stop("data must be a data.frame when fitting a reference model.")
    }

    factor_names <- .brsm_validate_factor_names(factor_names)

    if (is.null(family)) {
      family <- stats::gaussian()
    }

    if (!is.null(control) && !is.list(control)) {
      stop("control must be NULL or a named list.")
    }

    preset_control <- switch(sampling_preset,
      fast = list(adapt_delta = 0.8, max_treedepth = 10),
      balanced = list(adapt_delta = 0.9, max_treedepth = 12),
      robust = list(adapt_delta = 0.99, max_treedepth = 15)
    )

    if (is.null(control)) {
      control_use <- preset_control
    } else {
      control_use <- utils::modifyList(preset_control, control)
    }

    ref_formula <- .brsm_build_reference_formula(
      response = response,
      factor_names = factor_names,
      reference_type = reference_type
    )

    if (!is.null(seed)) {
      set.seed(seed)
    }

    brm_args <- list(
      formula = ref_formula,
      data = data,
      prior = prior,
      family = family,
      chains = chains,
      iter = iter,
      warmup = warmup,
      seed = seed,
      control = control_use,
      ...
    )
    if (!is.null(backend)) {
      brm_args$backend <- backend
    }

    reference_model <- do.call(brms::brm, brm_args)
    reference_fitted <- TRUE
  }

  reference_fit <- .brsm_extract_fit(reference_model, caller = "loftest_brsm")

  compare_args <- list(
    models = list(baseline = baseline_fit, reference = reference_fit),
    criterion = criterion
  )
  if (criterion == "loo") {
    compare_args$moment_match <- loo_moment_match
    compare_args$reloo <- loo_reloo
  }

  comparison <- do.call(compare_brsm_models, compare_args)

  loo_diagnostics <- NULL
  auto_retry_used <- FALSE
  if (criterion == "loo") {
    initial_diag <- .brsm_summarize_pareto_k(
      estimates = comparison$estimates,
      threshold = loo_k_threshold
    )

    if (isTRUE(loo_auto_moment_match) &&
        !isTRUE(loo_moment_match) &&
        isTRUE(initial_diag$has_high_k)) {
      message(
        "Found Pareto-k values above ", loo_k_threshold,
        "; retrying LOO with moment_match = TRUE."
      )
      compare_args$moment_match <- TRUE
      comparison <- do.call(compare_brsm_models, compare_args)
      auto_retry_used <- TRUE
      loo_moment_match <- TRUE
    }

    final_diag <- .brsm_summarize_pareto_k(
      estimates = comparison$estimates,
      threshold = loo_k_threshold
    )

    loo_diagnostics <- list(
      k_threshold = loo_k_threshold,
      initial = initial_diag,
      final = final_diag,
      moment_match_used = isTRUE(loo_moment_match),
      reloo_used = isTRUE(loo_reloo),
      auto_moment_match_retry = auto_retry_used
    )
  }

  out <- list(
    criterion = criterion,
    comparison = comparison,
    baseline_model = baseline_fit,
    reference_model = reference_fit,
    reference_type = reference_type,
    reference_fitted = reference_fitted
  )

  if (!is.null(loo_diagnostics)) {
    out$loo_diagnostics <- loo_diagnostics
  }

  if (isTRUE(include_ppc)) {
    ppc_base <- check_brsm_ppc(
      baseline_fit,
      ndraws = ppc_ndraws,
      probs = ppc_probs,
      seed = seed,
      include_plot = FALSE
    )
    ppc_ref <- check_brsm_ppc(
      reference_fit,
      ndraws = ppc_ndraws,
      probs = ppc_probs,
      seed = seed,
      include_plot = FALSE
    )

    ppc_summaries <- rbind(
      cbind(model = "baseline", ppc_base$summary, row.names = NULL),
      cbind(model = "reference", ppc_ref$summary, row.names = NULL)
    )
    rownames(ppc_summaries) <- NULL

    out$ppc <- list(
      summaries = ppc_summaries,
      baseline = ppc_base,
      reference = ppc_ref
    )
  }

  out
}


#' @keywords internal
.brsm_summarize_pareto_k <- function(estimates, threshold = 0.7) {
  if (!is.list(estimates) || length(estimates) == 0) {
    return(list(
      has_high_k = FALSE,
      max_pareto_k = NA_real_,
      max_pareto_k_by_model = numeric(0),
      n_above_threshold_by_model = integer(0)
    ))
  }

  model_names <- names(estimates)
  if (is.null(model_names) || any(model_names == "")) {
    model_names <- paste0("model_", seq_along(estimates))
  }

  max_by_model <- stats::setNames(rep(NA_real_, length(estimates)), model_names)
  n_above_by_model <- stats::setNames(rep(0L, length(estimates)), model_names)

  for (i in seq_along(estimates)) {
    est <- estimates[[i]]
    pareto_k <- tryCatch(est$diagnostics$pareto_k, error = function(e) NULL)
    if (is.null(pareto_k)) {
      next
    }

    pareto_k <- as.numeric(pareto_k)
    pareto_k <- pareto_k[is.finite(pareto_k)]
    if (length(pareto_k) == 0L) {
      next
    }

    max_by_model[i] <- max(pareto_k)
    n_above_by_model[i] <- sum(pareto_k > threshold)
  }

  max_overall <- if (all(is.na(max_by_model))) {
    NA_real_
  } else {
    max(max_by_model, na.rm = TRUE)
  }

  list(
    has_high_k = any(n_above_by_model > 0L),
    max_pareto_k = max_overall,
    max_pareto_k_by_model = max_by_model,
    n_above_threshold_by_model = n_above_by_model
  )
}


#' @keywords internal
.brsm_build_reference_formula <- function(
    response,
    factor_names,
    reference_type) {
  linear_terms <- factor_names

  interaction_terms <- character(0)
  if (length(factor_names) > 1) {
    combos <- utils::combn(factor_names, 2, simplify = FALSE)
    interaction_terms <- vapply(combos, function(pair) {
      paste0(pair[[1]], ":", pair[[2]])
    }, character(1))
  }

  quadratic_terms <- paste0("I(", factor_names, "^2)")
  cubic_terms <- paste0("I(", factor_names, "^3)")

  extra_terms <- character(0)
  if (reference_type == "extended" && length(factor_names) > 1) {
    combos <- utils::combn(factor_names, 2, simplify = FALSE)
    extra_terms <- unlist(lapply(combos, function(pair) {
      c(
        paste0("I(", pair[[1]], "^2):", pair[[2]]),
        paste0(pair[[1]], ":I(", pair[[2]], "^2)")
      )
    }), use.names = FALSE)
  }

  rhs_terms <- c(
    linear_terms, interaction_terms,
    quadratic_terms, cubic_terms, extra_terms
  )
  formula_text <- paste(response, "~", paste(rhs_terms, collapse = " + "))
  stats::as.formula(formula_text)
}