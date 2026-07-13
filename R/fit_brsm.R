#' Fit a Bayesian Quadratic Response Surface Model
#'
#' Fits a second-order (quadratic) response surface model using
#' \code{brms::brm()} and returns a \code{brsm_fit} object that stores
#' the fitted model and metadata for downstream brsm analysis.
#'
#' @param data Data frame containing response and factor columns.
#' @param response Name of the response variable.
#' @param factor_names Character vector of factor names.
#' @param ranges Optional named list of factor ranges. If \code{NULL},
#'   ranges are inferred from \code{data}.
#' @param prior Optional prior specification passed to \code{brms::brm()}.
#' @param prior_profile Prior profile to use when \code{prior = NULL}. One of
#'   \code{"legacy"} (flat, wide normal priors), \code{"regularized"}
#'   (moderate shrinkage), or \code{"adaptive"} (data-scaled priors).
#'   Default is \code{"legacy"}.
#' @param family Model family passed to \code{brms::brm()}.
#' @param chains Number of MCMC chains.
#' @param iter Number of iterations per chain.
#' @param warmup Number of warmup iterations per chain.
#' @param seed Optional random seed.
#' @param sampling_preset Sampling profile for NUTS tuning. One of
#'   \code{"fast"}, \code{"balanced"}, or \code{"robust"}.
#' @param backend Optional backend passed to \code{brms::brm()} (e.g.,
#'   \code{"rstan"}, \code{"cmdstanr"}). If \code{NULL}, uses brms defaults.
#' @param control Optional named list of NUTS control arguments. Values here
#'   override defaults from \code{sampling_preset}.
#' @param model_terms Polynomial term specification. One of
#'   \code{"second_order"} (default; linear + two-way interactions +
#'   pure quadratic terms), \code{"first_order"} (linear only),
#'   \code{"first_order_twi"} (linear + two-way interactions), or
#'   \code{"pure_quadratic"} (linear + pure quadratic terms).
#' @param coding_policy How to handle missing coding metadata from
#'   [prepare_brsm_data()]. One of \code{"warn"} (default),
#'   \code{"error"}, or \code{"ignore"}.
#' @param ... Additional arguments passed to \code{brms::brm()}.
#'
#' @return An object of class \code{brsm_fit} with elements:
#'   \code{fit}, \code{formula}, \code{response}, \code{factor_names},
#'   \code{ranges}, \code{model_terms}, and \code{call}.
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(
#'   x1 = runif(100, -2, 2),
#'   x2 = runif(100, -2, 2)
#' )
#' dat$y <- 5 + 2 * dat$x1 + 4 * dat$x2 - dat$x1^2 - 2 * dat$x2^2 +
#'   0.5 * dat$x1 * dat$x2 + rnorm(100, sd = 0.5)
#'
#' fit <- fit_brsm(
#'   data = dat,
#'   response = "y",
#'   factor_names = c("x1", "x2"),
#'   sampling_preset = "balanced",
#'   control = list(adapt_delta = 0.95),
#'   chains = 2,
#'   iter = 1000,
#'   warmup = 500,
#'   seed = 123
#' )
#'
#' out <- brsm_workflow(
#'   object = fit,
#'   factor_names = fit$factor_names,
#'   ranges = fit$ranges
#' )
#' }
#' @export
fit_brsm <- function(data,
                     response,
                     factor_names,
                     ranges = NULL,
                     prior = NULL,
                     prior_profile = c("legacy", "regularized", "adaptive"),
                     draws = 2000,
                     seed = NULL,
                     model_terms = c(
                       "second_order", "first_order",
                       "first_order_twi", "pure_quadratic"
                     ),
                     coding_policy = c("warn", "error", "ignore"),
                     ...) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }

  coding <- attr(data, "brsm_coding", exact = TRUE)
  coding_policy <- match.arg(coding_policy)

  if (!is.character(response) || length(response) != 1) {
    stop("response must be a single character string.")
  }

  prior_profile <- match.arg(prior_profile)
  model_terms <- match.arg(model_terms)

  factor_names <- .brsm_validate_factor_names(factor_names)

  .brsm_enforce_coding_policy(
    data = data,
    factor_names = factor_names,
    coding = coding,
    coding_policy = coding_policy,
    caller = "fit_brsm"
  )

  model_vars <- c(response, factor_names)
  missing_vars <- setdiff(model_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "Missing required columns in data: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  if (!all(vapply(data[factor_names], is.numeric, logical(1)))) {
    stop("All factor columns must be numeric.")
  }

  if (!is.numeric(data[[response]])) {
    stop("response column must be numeric.")
  }

  if (is.null(ranges)) {
    ranges <- stats::setNames(
      lapply(factor_names, function(f) range(data[[f]], na.rm = TRUE)),
      factor_names
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  fit <- .brsm_fit_conjugate(
    data = data,
    response = response,
    factor_names = factor_names,
    ranges = ranges,
    prior = prior,
    prior_profile = prior_profile,
    draws = draws,
    model_terms = model_terms
  )

  out <- list(
    fit = fit,
    formula = fit$formula,
    response = response,
    factor_names = factor_names,
    ranges = ranges,
    coding = coding,
    model_terms = model_terms,
    sampling = list(
      prior_profile = prior_profile,
      draws = draws,
      seed = seed
    ),
    call = match.call()
  )

  class(out) <- c("brsm_fit", "list")
  out
}