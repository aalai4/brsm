#' Fit Constrained Congruence Models for RSA Hypothesis Testing
#'
#' Fits Bayesian quadratic response surface models with congruence constraints,
#' enabling direct testing of RSA hypothesis ordering: unconstrained > broad
#' congruence > strict congruence.
#'
#' Congruence hypotheses impose coefficient constraints:
#'   - \code{"unconstrained"}: Full quadratic model (y ~ x1 + x2 + I(x1^2) +
#'     I(x2^2) + x1:x2)
#'   - \code{"broad"}: b1 = b2 and b11 = b22, but interaction free
#'     (y ~ I(x1 + x2) + I(I(x1^2) + I(x2^2)) + x1:x2)
#'   - \code{"strict"}: All broad constraints PLUS b12 = 0 (no interaction)
#'     (y ~ I(x1 + x2) + I(I(x1^2) + I(x2^2)))
#'
#' The broad and strict specifications enforce congruence hypotheses: outcomes
#' depend only on X-Y congruence/fit, not on the direction of discrepancy.
#'
#' @param data Data frame containing response and factor columns.
#' @param response Name of the response variable.
#' @param factor_names Character vector of exactly two factor names.
#' @param congruence_type One of \code{"unconstrained"}, \code{"broad"}, or
#'   \code{"strict"}. Default is \code{"unconstrained"}. Determines which
#'   coefficient constraints are applied.
#' @param ranges Optional named list of factor ranges.
#' @param prior Optional prior specification passed to \code{brms::brm()}.
#' @param family Model family passed to \code{brms::brm()}.
#' @param chains Number of MCMC chains.
#' @param iter Number of iterations per chain.
#' @param warmup Number of warmup iterations per chain.
#' @param seed Optional random seed.
#' @param sampling_preset Sampling profile for NUTS tuning. One of
#'   \code{"fast"}, \code{"balanced"}, or \code{"robust"}.
#' @param backend Optional backend passed to \code{brms::brm()}.
#' @param control Optional named list of NUTS control arguments.
#' @param coding_policy How to handle missing coding metadata from
#'   [prepare_brsm_data()]. One of \code{"warn"} (default),
#'   \code{"error"}, or \code{"ignore"}.
#' @param ... Additional arguments passed to \code{brms::brm()}.
#'
#' @return An object of class \code{brsm_fit} with additional metadata:
#'   \code{congruence_type} field records the constraint specification.
#'   All other elements (\code{fit}, \code{formula}, \code{response},
#'   \code{factor_names}, \code{ranges}) follow [fit_brsm()] convention.
#'
#' @details
#' Broad and strict specifications are used to test RSA congruence hypotheses
#' by model comparison (LOO/WAIC). Typical workflow:
#'
#' 1. Fit all three models (unconstrained, broad, strict)
#' 2. Compare using [compare_brsm_models()]
#' 3. Compute congruence parameters [congruence_parameters()] for each
#' 4. Apply ROPE tests [rope_congruence()] to interpret coefficients
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(
#'   x1 = runif(100, -2, 2),
#'   x2 = runif(100, -2, 2)
#' )
#' dat$y <- 5 + 2 * (dat$x1 + dat$x2) + (-1) * (dat$x1^2 + dat$x2^2) +
#'   0.3 * dat$x1 * dat$x2 + rnorm(100, sd = 0.5)
#'
#' # Fit all three models
#' fit_unconstrained <- fit_brsm_congruence(
#'   data = dat,
#'   response = "y",
#'   factor_names = c("x1", "x2"),
#'   congruence_type = "unconstrained",
#'   chains = 2, iter = 1000
#' )
#'
#' fit_broad <- fit_brsm_congruence(
#'   data = dat,
#'   response = "y",
#'   factor_names = c("x1", "x2"),
#'   congruence_type = "broad",
#'   chains = 2, iter = 1000
#' )
#'
#' fit_strict <- fit_brsm_congruence(
#'   data = dat,
#'   response = "y",
#'   factor_names = c("x1", "x2"),
#'   congruence_type = "strict",
#'   chains = 2, iter = 1000
#' )
#'
#' # Compare models
#' comparisons <- compare_brsm_models(
#'   list(
#'     unconstrained = fit_unconstrained,
#'     broad = fit_broad,
#'     strict = fit_strict
#'   ),
#'   criterion = "loo"
#' )
#' print(comparisons)
#' }
#'
#' @export
fit_brsm_congruence <- function(data,
                                response,
                                factor_names,
                                congruence_type = c(
                                  "unconstrained", "broad", "strict"
                                ),
                                ranges = NULL,
                                prior = NULL,
                                family = stats::gaussian(),
                                chains = 4,
                                iter = 2000,
                                warmup = floor(iter / 2),
                                seed = NULL,
                                sampling_preset = c(
                                  "fast", "balanced", "robust"
                                ),
                                backend = NULL,
                                control = NULL,
                                coding_policy = c("warn", "error", "ignore"),
                                ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("package 'brms' is required for fit_brsm_congruence().")
  }

  congruence_type <- match.arg(congruence_type)
  coding_policy <- match.arg(coding_policy)

  # Basic validation
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }

  factor_names <- .brsm_validate_factor_names(
    factor_names,
    require_two = TRUE,
    two_factors_message = "fit_brsm_congruence requires exactly two factors."
  )

  coding <- attr(data, "brsm_coding", exact = TRUE)

  .brsm_enforce_coding_policy(
    data = data,
    factor_names = factor_names,
    coding = coding,
    coding_policy = coding_policy,
    caller = "fit_brsm_congruence"
  )

  if (!is.character(response) || length(response) != 1) {
    stop("response must be a character scalar.")
  }

  if (!(response %in% colnames(data))) {
    stop("response '", response, "' not found in data columns.")
  }

  # Build constrained formula based on type
  formula_str <- .brsm_build_congruence_formula(
    response = response,
    factor_names = factor_names,
    congruence_type = congruence_type
  )

  model_formula <- stats::as.formula(formula_str)

  # Set up brms arguments
  sampling_preset <- match.arg(sampling_preset)

  # Sampling preset defaults
  preset_controls <- list(
    fast = list(adapt_delta = 0.8, max_treedepth = 10),
    balanced = list(adapt_delta = 0.9, max_treedepth = 12),
    robust = list(adapt_delta = 0.99, max_treedepth = 15)
  )

  control_default <- preset_controls[[sampling_preset]]
  if (!is.null(control)) {
    control_default <- utils::modifyList(control_default, control)
  }

  brm_args <- list(
    formula = model_formula,
    data = data,
    prior = prior,
    family = family,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    control = control_default,
    ...
  )
  if (!is.null(backend)) {
    brm_args$backend <- backend
  }

  fit <- do.call(brms::brm, brm_args)

  # Preserve coding metadata from input data if present

  out <- list(
    fit = fit,
    formula = model_formula,
    response = response,
    factor_names = factor_names,
    ranges = ranges,
    coding = coding,
    congruence_type = congruence_type,
    sampling = list(
      sampling_preset = sampling_preset,
      control = control_default,
      chains = chains,
      iter = iter,
      warmup = warmup,
      seed = seed,
      backend = backend
    ),
    call = match.call()
  )

  class(out) <- c("brsm_fit", "list")
  out
}


#' Build Congruence-Constrained Formula
#'
#' Internal function to construct a model formula with appropriate
#' congruence constraints for RSA hypothesis testing.
#'
#' @param response Character scalar, response variable name.
#' @param factor_names Character vector of two factor names.
#' @param congruence_type Type of constraint: "unconstrained",
#'   "broad", "strict".
#'
#' @return Character string representing the constrained formula.
#'
#' @keywords internal
.brsm_build_congruence_formula <- function(
    response,
    factor_names,
    congruence_type) {
  f1 <- factor_names[1]
  f2 <- factor_names[2]

  if (congruence_type == "unconstrained") {
    # Full quadratic: y ~ x1 + x2 + I(x1^2) + I(x2^2) + x1:x2
    formula_str <- paste(
      response,
      "~",
      f1, "+", f2, "+",
      "I(", f1, "^2) + I(", f2, "^2) +",
      f1, ":", f2,
      sep = ""
    )
  } else if (congruence_type == "broad") {
    # Broad congruence: y ~ I(x1 + x2) + I(I(x1^2) + I(x2^2)) + x1:x2
    # Enforces b1=b2 and b11=b22
    formula_str <- paste(
      response,
      "~",
      "I(", f1, "+", f2, ") +",
      "I(I(", f1, "^2) + I(", f2, "^2)) +",
      f1, ":", f2,
      sep = ""
    )
  } else if (congruence_type == "strict") {
    # Strict congruence: y ~ I(x1 + x2) + I(I(x1^2) + I(x2^2))
    # Enforces b1=b2, b11=b22, AND b12=0
    formula_str <- paste(
      response,
      "~",
      "I(", f1, "+", f2, ") +",
      "I(I(", f1, "^2) + I(", f2, "^2))",
      sep = ""
    )
  } else {
    stop("congruence_type must be one of: unconstrained, broad, strict")
  }

  formula_str
}