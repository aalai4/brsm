#' Run a Complete BRSM Analysis Workflow
#'
#' Executes a full response-surface workflow from model/draw conversion through
#' prediction, optimization diagnostics, and optional plots.
#'
#' By default (\code{fit_mode = "none"}), input must be Bayesian:
#' \code{brsm_fit}, \code{brmsfit}, or posterior-draw data frame with required
#' coefficient columns (\code{b_Intercept}, \code{b_x1}, \code{b_I(x1^2)},
#' interactions, etc.).
#'
#' If \code{fit_mode} is \code{"fit_brsm"} or \code{"fit_brsm_congruence"},
#' then \code{object} is treated as raw input data and model fitting is done
#' internally before running selected workflow steps.
#'
#' @param object Either (1) a \code{brsm_fit} object from [fit_brsm()],
#'   (2) a \code{brmsfit} object from \code{brms::brm()},
#'   (3) a posterior-draw data frame with Bayesian coefficient columns, or
#'   (4) raw analysis data when \code{fit_mode != "none"}.
#' @param factor_names Character vector of factor names.
#' @param ranges Named list with factor ranges.
#' @param probs Posterior probabilities for prediction summaries.
#' @param grid_n Number of points per factor dimension in the prediction grid.
#' @param radii Numeric vector of radii for ridge analysis.
#' @param start Optional named numeric vector of starting values for steepest
#'   ascent. Defaults to zeros for each factor.
#' @param step_size Step size for steepest ascent.
#' @param n_steps Number of steps for steepest ascent.
#' @param include_plots Logical; if `TRUE`, generate contour/optimum/ridge plots
#'   when exactly two factors are provided.
#' @param contour_type Surface type for contour plotting. Passed to
#'   [plot_posterior_contours()].
#' @param seed Optional random seed used for stochastic plotting operations.
#' @param fit_mode One of `"none"`, `"fit_brsm"`, or
#'   `"fit_brsm_congruence"`. Use `"none"` (default) to analyze an existing
#'   Bayesian fit/draw object. Use a fit mode to fit from raw data in `object`
#'   and then run the workflow.
#' @param response Required when `fit_mode != "none"`. Name of the response
#'   column in raw input data.
#' @param fit_args Named list of additional arguments forwarded to
#'   [fit_brsm()] or [fit_brsm_congruence()] when `fit_mode != "none"`.
#'   This enables custom priors and sampling controls without re-writing the
#'   pipeline manually.
#' @param steps Character vector selecting which workflow stages to run. Any of
#'   `"predictions"`, `"stationary"`, `"classification"`,
#'   `"credible_region"`, `"steepest_ascent"`, `"ridge"`, `"plots"`.
#'   Default runs all stages.
#'
#' @return A named list with components:
#'   always `draws`, and conditionally `fit`, `ranges`, `grid`, `predictions`,
#'   `surface`, `stationary_points`, `classification`, `credible_region`,
#'   `steepest_ascent`, `ridge_analysis`, and `plots` depending on `steps`.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   x1 = runif(50, -2, 2),
#'   x2 = runif(50, -2, 2)
#' )
#' dat$y <- 5 + 2 * dat$x1 + 4 * dat$x2 - dat$x1^2 - 2 * dat$x2^2 +
#'   0.5 * dat$x1 * dat$x2 + rnorm(50, sd = 0.5)
#' draws_raw <- data.frame(
#'   b_Intercept = rnorm(200, 5, 0.2),
#'   b_x1 = rnorm(200, 2, 0.1),
#'   b_x2 = rnorm(200, 4, 0.1),
#'   "b_I(x1^2)" = rnorm(200, -1, 0.05),
#'   "b_I(x2^2)" = rnorm(200, -2, 0.05),
#'   "b_x1:x2" = rnorm(200, 0.5, 0.05),
#'   check.names = FALSE
#' )
#'
#' out <- brsm_workflow(
#'   object = draws_raw,
#'   factor_names = c("x1", "x2"),
#'   ranges = list(x1 = c(-2, 2), x2 = c(-2, 2))
#' )
#'
#' names(out)
#' @export
brsm_workflow <- function(object,
                          factor_names,
                          ranges = NULL,
                          probs = c(0.025, 0.5, 0.975),
                          grid_n = 40,
                          radii = seq(0, 1, length.out = 5),
                          start = NULL,
                          step_size = 0.1,
                          n_steps = 10,
                          include_plots = TRUE,
                          contour_type = c("mean", "uncertainty", "quantile"),
                          seed = NULL,
                          fit_mode = c(
                            "none", "fit_brsm", "fit_brsm_congruence"
                          ),
                          response = NULL,
                          fit_args = list(),
                          steps = c(
                            "predictions", "stationary", "classification",
                            "credible_region", "steepest_ascent",
                            "ridge", "plots"
                          )) {
  factor_names <- .brsm_validate_factor_names(factor_names)
  contour_type <- match.arg(contour_type)
  fit_mode <- match.arg(fit_mode)

  if (!is.list(fit_args)) {
    stop("fit_args must be a named list.")
  }

  if (!is.character(steps) || length(steps) == 0L) {
    stop("steps must be a non-empty character vector.")
  }

  allowed_steps <- c(
    "predictions", "stationary", "classification",
    "credible_region", "steepest_ascent", "ridge", "plots"
  )
  bad_steps <- setdiff(steps, allowed_steps)
  if (length(bad_steps) > 0L) {
    stop(
      "Unknown step(s): ", paste(bad_steps, collapse = ", "),
      ". Allowed: ", paste(allowed_steps, collapse = ", "), "."
    )
  }
  steps <- unique(steps)

  if (isTRUE(include_plots) && !("plots" %in% steps)) {
    steps <- c(steps, "plots")
  }
  if (!isTRUE(include_plots)) {
    steps <- setdiff(steps, "plots")
  }

  # Resolve workflow input: either fit from raw data or validate Bayesian input.
  workflow_object <- object
  fitted_model <- NULL

  if (fit_mode == "none") {
    .brsm_validate_bayesian_input(workflow_object)
  } else {
    if (!is.data.frame(workflow_object)) {
      stop("When fit_mode != 'none', object must be a raw data.frame.")
    }
    if (is.null(response) || !is.character(response) ||
          length(response) != 1L) {
      stop(
        "response must be supplied as a character scalar",
        " when fit_mode != 'none'."
      )
    }

    base_fit_args <- list(
      data = workflow_object,
      response = response,
      factor_names = factor_names,
      ranges = ranges
    )
    merged_fit_args <- utils::modifyList(base_fit_args, fit_args)

    if (fit_mode == "fit_brsm") {
      fitted_model <- do.call(fit_brsm, merged_fit_args)
    } else {
      fitted_model <- do.call(fit_brsm_congruence, merged_fit_args)
    }
    workflow_object <- fitted_model
  }

  if (!is.numeric(grid_n) || length(grid_n) != 1 ||
        !is.finite(grid_n) || grid_n < 2) {
    stop("grid_n must be a finite integer >= 2.")
  }
  grid_n <- as.integer(grid_n)

  probs <- .brsm_validate_probs(probs)

  if (is.null(start)) {
    start <- stats::setNames(rep(0, length(factor_names)), factor_names)
  }

  draws <- as_brsm_draws(workflow_object, factor_names = factor_names)

  draws_complete <- stats::complete.cases(draws)
  if (!all(draws_complete)) {
    n_dropped <- sum(!draws_complete)
    warning(
      "Dropped ", n_dropped,
      " draw row(s) with missing/non-finite coefficients",
      " before analysis."
    )
    draws <- draws[draws_complete, , drop = FALSE]
  }
  if (nrow(draws) == 0) {
    stop(
      "No valid draw rows remain after filtering",
      " missing/non-finite coefficients. ",
      "Check model specification for rank deficiency or missing coefficients."
    )
  }

  if (is.null(ranges) &&
    inherits(workflow_object, "brsm_fit") &&
    !is.null(workflow_object$ranges)) {
    ranges <- workflow_object$ranges
  }

  # Determine whether ranges/grid are needed.
  needs_ranges <- any(c("predictions", "plots") %in% steps)
  if (needs_ranges && is.null(ranges)) {
    stop(
      "ranges must be supplied (or available in brsm_fit$ranges)",
      " when steps include predictions or plots."
    )
  }

  out <- list(draws = draws)
  if (!is.null(fitted_model)) {
    out$fit <- fitted_model
  }
  if (!is.null(ranges)) {
    out$ranges <- ranges
  }

  grid <- NULL
  predictions <- NULL
  stationary <- NULL
  ridge_full <- NULL

  dispatch_object <- if (inherits(workflow_object, "brsm_fit")) {
    workflow_object
  } else {
    draws
  }

  if ("predictions" %in% steps) {
    grid <- surface_grid(ranges = ranges, n = grid_n)
    predictions <- predict_surface(
      draws = dispatch_object,
      factor_names = factor_names,
      newdata = grid,
      summary = TRUE,
      probs = probs
    )
    out$grid <- grid
    out$predictions <- predictions
    # Keep `surface` for backward compatibility.
    out$surface <- predictions
  }

  needs_stationary <- "stationary" %in% steps ||
    "classification" %in% steps ||
    "credible_region" %in% steps ||
    "plots" %in% steps
  if (needs_stationary) {
    stationary <- stationary_point(
      object = dispatch_object,
      factor_names = factor_names
    )
  }

  if ("stationary" %in% steps) {
    out$stationary_points <- stationary
  }

  if ("classification" %in% steps) {
    out$classification <- classify_stationarity_point(
      object = dispatch_object,
      factor_names = factor_names
    )
  }

  if ("credible_region" %in% steps) {
    out$credible_region <- credible_optimum_region(
      object = dispatch_object,
      factor_names = factor_names,
      probs = probs[c(1, length(probs))],
      summary = TRUE
    )
  }

  if ("steepest_ascent" %in% steps) {
    out$steepest_ascent <- steepest_ascent(
      object = dispatch_object,
      factor_names = factor_names,
      start = start,
      step_size = step_size,
      n_steps = n_steps
    )
  }

  if ("ridge" %in% steps || "plots" %in% steps) {
    ridge_full <- posterior_ridge_analysis(
      object = dispatch_object,
      factor_names = factor_names,
      radii = radii,
      summary = FALSE
    )
  }

  if ("ridge" %in% steps) {
    out$ridge_analysis <- posterior_ridge_analysis(
      object = dispatch_object,
      factor_names = factor_names,
      radii = radii,
      summary = TRUE
    )
  }

  if ("plots" %in% steps) {
    if (length(factor_names) == 2 && !is.null(ranges)) {
      if (is.null(stationary)) {
        stationary <- stationary_point(
          object = dispatch_object,
          factor_names = factor_names
        )
      }
      plots <- list(
        contours = plot_posterior_contours(
          draws = draws,
          factor_names = factor_names,
          ranges = ranges,
          n = grid_n,
          type = contour_type,
          overlay_stationary = TRUE,
          stationary_draws = stationary,
          seed = seed
        ),
        optimum = plot_optimum_posterior(
          draws = draws,
          factor_names = factor_names,
          stationary_draws = stationary,
          seed = seed
        )
      )

      if (!is.null(ridge_full)) {
        ridge_df <- ridge_full[, c(factor_names, "radius"), drop = FALSE]
        ridge_df <- ridge_df[stats::complete.cases(ridge_df), , drop = FALSE]

        if (nrow(ridge_df) > 0) {
          plots$ridge_path <- plot_ridge_path(
            ridge_draws = ridge_df,
            factor_names = factor_names,
            response = "radius",
            show_draws = TRUE
          )
        }
      }

      out$plots <- plots
    } else {
      warning(
        "plots step ignored because plotting currently",
        " requires exactly two factors and valid ranges."
      )
    }
  }

  out
}