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
#'   \code{reference_fitted} (logical), and optional \code{ppc} list with
#'   side-by-side summaries.
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

  comparison <- compare_brsm_models(
    models = list(baseline = baseline_fit, reference = reference_fit),
    criterion = criterion
  )

  out <- list(
    criterion = criterion,
    comparison = comparison,
    baseline_model = baseline_fit,
    reference_model = reference_fit,
    reference_type = reference_type,
    reference_fitted = reference_fitted
  )

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