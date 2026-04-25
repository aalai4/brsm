#' Run a Complete BRSM Analysis Workflow
#'
#' Runs a response-surface workflow from model/draw conversion through
#' prediction, diagnostics, and optional plots.
#'
#' By default (\code{fit_mode = "none"}), input must be Bayesian:
#' \code{brsm_fit}, \code{brmsfit}, or posterior-draw data frame with required
#' coefficient columns (\code{b_Intercept}, \code{b_x1}, \code{b_I(x1^2)},
#' interactions, etc.).
#'
#' If \code{fit_mode} is \code{"fit_brsm"}, \code{object} is treated as raw
#' input data and fitting is performed before selected workflow steps run.
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
#'   when exactly two factors are provided. For non-quadratic draw objects,
#'   contour plots are returned with a steepest-ascent direction overlay and
#'   optimization-geometry plots are omitted.
#' @param contour_type Surface type for contour plotting. Passed to
#'   [plot_posterior_contours()].
#' @param plot_vary_factors Optional character vector of length 2 specifying
#'   which factors vary in contour plots when `length(factor_names) > 2` and
#'   `plot_pairwise = FALSE`. Defaults to the first two factors.
#' @param plot_fixed Optional named numeric vector of fixed values used when
#'   `plot_conditioning = "user"`.
#' @param plot_conditioning Conditioning strategy for non-varied factors in
#'   contour plots. One of `"center"`, `"optimum_mean"`, or `"user"`.
#' @param plot_slice Optional named list with one element specifying sliced
#'   values for a conditioning factor, e.g. `list(x3 = c(-1, 0, 1))`.
#' @param plot_pairwise Logical; if `TRUE`, generate pairwise contour panels for
#'   all factor pairs.
#' @param seed Optional random seed used for stochastic plotting operations.
#' @param fit_mode One of `"none"` or `"fit_brsm"`. Use `"none"` (default) to analyze an existing
#'   Bayesian fit/draw object. Use `"fit_brsm"` to fit from raw data in
#'   \code{object} and then run the workflow.
#' @param response Required when `fit_mode != "none"`. Name of the response
#'   column in raw input data.
#' @param fit_args Named list of additional arguments forwarded to
#'   [fit_brsm()] when `fit_mode != "none"`.
#'   Useful for custom priors and sampling controls.
#' @param steps Character vector selecting which workflow stages to run. Any of
#'   `"predictions"`, `"stationary"`, `"classification"`,
#'   `"credible_region"`, `"steepest_ascent"`, `"ridge"`, `"plots"`.
#'   Default runs all stages.
#'
#' @return A named list with components:
#'   always `draws`, and conditionally `fit`, `ranges`, `grid`, `predictions`,
#'   `surface`, `stationary_points`, `classification`, `credible_region`,
#'   `steepest_ascent`, `ridge_analysis`, and `plots` depending on `steps`.
#'   When plots are produced, `plots` also includes
#'   `plot_mode` (`"optimization_geometry"` or `"direction_only"`) and
#'   `geometry_available` (logical).
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
                          plot_vary_factors = NULL,
                          plot_fixed = NULL,
                          plot_conditioning = c(
                            "center", "optimum_mean", "user"
                          ),
                          plot_slice = NULL,
                          plot_pairwise = FALSE,
                          seed = NULL,
                          fit_mode = c("none", "fit_brsm"),
                          response = NULL,
                          fit_args = list(),
                          steps = c(
                            "predictions", "stationary", "classification",
                            "credible_region", "steepest_ascent",
                            "ridge", "plots"
                          )) {
  factor_names <- .brsm_validate_factor_names(factor_names)
  contour_type <- match.arg(contour_type)
  plot_conditioning <- match.arg(plot_conditioning)
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

    fitted_model <- do.call(fit_brsm, merged_fit_args)
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

  draws <- if (is.data.frame(workflow_object)) {
    as_brsm_draws.data.frame(
      workflow_object,
      factor_names = factor_names,
      require_quadratic = FALSE,
      require_interactions = FALSE
    )
  } else {
    as_brsm_draws(workflow_object, factor_names = factor_names)
  }

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

  has_full_quadratic <- all(vapply(
    factor_names,
    function(f) {
      qcol <- .brsm_find_quadratic_col(f, colnames(draws))
      !is.na(qcol) && qcol != ""
    },
    logical(1)
  ))

  geometry_steps <- c("stationary", "classification", "credible_region", "ridge")
  unsupported_steps <- intersect(steps, geometry_steps)
  if (!has_full_quadratic && length(unsupported_steps) > 0) {
    warning(
      "Skipping geometry-dependent step(s) because quadratic terms are ",
      "missing in draws: ",
      paste(unsupported_steps, collapse = ", "),
      "."
    )
    steps <- setdiff(steps, unsupported_steps)
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

  needs_stationary <- has_full_quadratic && (
    "stationary" %in% steps ||
      "classification" %in% steps ||
      "credible_region" %in% steps ||
      "plots" %in% steps
  )
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

  if (has_full_quadratic && ("ridge" %in% steps || "plots" %in% steps)) {
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
    if (length(factor_names) >= 2 && !is.null(ranges)) {
      overlay_stationary <- has_full_quadratic
      contour_overlay_stationary <- overlay_stationary &&
        !isTRUE(plot_pairwise) &&
        is.null(plot_slice)

      selected_plot_pair <- if (!is.null(plot_vary_factors)) {
        .brsm_validate_factor_names(
          plot_vary_factors,
          require_two = TRUE,
          two_factors_message =
            "plot_vary_factors must contain exactly two factors."
        )
      } else {
        factor_names[1:2]
      }

      if (!all(selected_plot_pair %in% factor_names)) {
        stop("plot_vary_factors must be a subset of factor_names.")
      }

      optimization_pairs <- if (length(factor_names) == 2L) {
        list(factor_names)
      } else if (isTRUE(plot_pairwise)) {
        utils::combn(factor_names, 2, simplify = FALSE)
      } else {
        list(selected_plot_pair)
      }

      if (overlay_stationary && is.null(stationary)) {
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
          vary_factors = plot_vary_factors,
          fixed = plot_fixed,
          conditioning = plot_conditioning,
          slice = plot_slice,
          pairwise = plot_pairwise,
          overlay_stationary = contour_overlay_stationary,
          stationary_draws = stationary,
          seed = seed
        )
      )

      plots$plot_mode <- if (overlay_stationary) {
        "optimization_geometry"
      } else {
        "direction_only"
      }
      plots$geometry_available <- overlay_stationary

      if (overlay_stationary && length(factor_names) == 2L) {
        plots$optimum <- plot_optimum_posterior(
          draws = draws,
          factor_names = factor_names,
          stationary_draws = stationary,
          seed = seed
        )
      } else if (overlay_stationary && length(factor_names) > 2L) {
        optimum_plots <- lapply(optimization_pairs, function(pair) {
          stationary_pair <- stationary[, pair, drop = FALSE]
          p <- plot_optimum_posterior(
            draws = draws,
            factor_names = pair,
            stationary_draws = stationary_pair,
            seed = seed
          )
          p + ggplot2::labs(
            subtitle = paste0("Projection: ", pair[1], " vs ", pair[2])
          )
        })
        names(optimum_plots) <- vapply(
          optimization_pairs,
          function(pair) paste(pair, collapse = "_"),
          character(1)
        )

        if (length(optimum_plots) == 1L) {
          plots$optimum <- optimum_plots[[1L]]
        } else {
          plots$optimum_pairwise <- optimum_plots
        }
      } else {
        warning(
          "Skipping optimum/ridge plot overlays because stationary points ",
          "are undefined without quadratic terms."
        )

        sa_plot <- if (!is.null(out$steepest_ascent)) {
          out$steepest_ascent
        } else {
          steepest_ascent(
            object = dispatch_object,
            factor_names = factor_names,
            start = start,
            step_size = step_size,
            n_steps = n_steps
          )
        }

        mean_path <- sa_plot$mean_path
        if (!is.null(mean_path) &&
          nrow(mean_path) > 1 &&
          all(factor_names %in% names(mean_path))) {
          plots$contours <- plots$contours +
            ggplot2::geom_path(
              data = mean_path,
              ggplot2::aes(
                x = .data[[factor_names[1]]],
                y = .data[[factor_names[2]]]
              ),
              inherit.aes = FALSE,
              color = "black",
              linewidth = 1
            ) +
            ggplot2::geom_point(
              data = mean_path,
              ggplot2::aes(
                x = .data[[factor_names[1]]],
                y = .data[[factor_names[2]]]
              ),
              inherit.aes = FALSE,
              color = "red",
              size = 2
            )

          plots$direction <- ggplot2::ggplot(
            mean_path,
            ggplot2::aes(
              x = .data[[factor_names[1]]],
              y = .data[[factor_names[2]]]
            )
          ) +
            ggplot2::geom_path(color = "black", linewidth = 1) +
            ggplot2::geom_point(color = "red", size = 2) +
            ggplot2::labs(
              x = factor_names[1],
              y = factor_names[2],
              title = "Posterior Mean Direction Path",
              subtitle = "Steepest-ascent mean path for non-quadratic model"
            ) +
            ggplot2::theme_minimal()
        } else {
          warning(
            "Direction overlay could not be drawn because a mean ",
            "steepest-ascent path was unavailable."
          )
        }
      }

      if (!is.null(ridge_full) && length(factor_names) == 2L) {
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
      } else if (!is.null(ridge_full) && length(factor_names) > 2L) {
        ridge_plots <- lapply(optimization_pairs, function(pair) {
          ridge_df <- ridge_full[, c(pair, "radius"), drop = FALSE]
          ridge_df <- ridge_df[stats::complete.cases(ridge_df), , drop = FALSE]
          if (nrow(ridge_df) == 0) {
            return(NULL)
          }
          p <- plot_ridge_path(
            ridge_draws = ridge_df,
            factor_names = pair,
            response = "radius",
            show_draws = TRUE
          )
          p + ggplot2::labs(
            subtitle = paste0("Projection: ", pair[1], " vs ", pair[2])
          )
        })
        names(ridge_plots) <- vapply(
          optimization_pairs,
          function(pair) paste(pair, collapse = "_"),
          character(1)
        )
        ridge_plots <- Filter(Negate(is.null), ridge_plots)

        if (length(ridge_plots) == 1L) {
          plots$ridge_path <- ridge_plots[[1L]]
        } else if (length(ridge_plots) > 1L) {
          plots$ridge_path_pairwise <- ridge_plots
        }
      }

      out$plots <- plots
    } else {
      warning(
        "plots step ignored because plotting currently",
        " requires at least two factors and valid ranges."
      )
    }
  }

  out
}