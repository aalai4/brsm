#' Plot Posterior Distribution of the Optimum
#'
#' @param draws Validated draw data frame.
#' @param factor_names Character vector of factor names (exactly 2).
#' @param stationary_draws Optional precomputed stationary draws.
#' @param bins Number of density bins.
#' @param levels Number of density contour levels.
#' @param alpha Transparency for density fill.
#' @param show_points Logical; if TRUE, overlay posterior sample points.
#' @param point_alpha Transparency for point overlay.
#' @param seed Random seed.
#'
#' @return A ggplot2 object.
#'
#' @details
#' Density fills use the viridis "turbo" palette, which is perceptually uniform
#' and colorblind-safe. Density contour lines are drawn in black for contrast.
#' The posterior mean optimum is marked with a red cross symbol. This plot is only
#' available for optimization_geometry mode (when quadratic terms are present).
#' @keywords internal
plot_optimum_posterior <- function(draws,
                                   factor_names,
                                   stationary_draws = NULL,
                                   bins = 60,
                                   levels = 10,
                                   alpha = 0.9,
                                   show_points = FALSE,
                                   point_alpha = 0.2,
                                   seed = NULL) {
  draws <- .brsm_validate_draws(draws)
  factor_names <- .brsm_validate_factor_names(
    factor_names,
    require_two = TRUE,
    two_factors_message =
      "plot_optimum_posterior currently supports exactly two factors."
  )

  .brsm_require_ggplot2()

  if (!is.numeric(bins) || bins <= 0) {
    stop("bins must be a positive integer.")
  }
  bins <- as.integer(bins)

  if (!is.numeric(levels) || levels <= 0) {
    stop("levels must be a positive integer.")
  }
  levels <- as.integer(levels)

  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
  }

  if (!is.numeric(point_alpha) || point_alpha < 0 || point_alpha > 1) {
    stop("point_alpha must be between 0 and 1.")
  }

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1)) {
    stop("seed must be a single numeric value.")
  }

  # Obtain stationary point draws
  if (is.null(stationary_draws)) {
    stationary <- stationary_point(draws, factor_names)
    stationary <- as.data.frame(stationary)
    colnames(stationary) <- factor_names
  } else {
    stationary <- as.data.frame(stationary_draws)
    .brsm_check_columns(
      factor_names, stationary,
      "stationary_draws must contain columns for factor_names."
    )
  }
  stationary <- stationary[, factor_names, drop = FALSE]

  # Remove invalid rows and warn if any removed
  removed <- sum(!stats::complete.cases(stationary))
  stationary <- stationary[stats::complete.cases(stationary), , drop = FALSE]
  if (removed > 0) {
    warning(removed, " stationary points removed due to missing values.")
  }

  if (nrow(stationary) == 0) {
    stop("no valid stationary points available for plotting.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Posterior mean optimum (index by name for safety)
  opt_mean <- colMeans(stationary, na.rm = TRUE)

  # Create plot
  p <- ggplot2::ggplot(
    stationary,
    ggplot2::aes(
      x = .data[[factor_names[1]]],
      y = .data[[factor_names[2]]]
    )
  ) +
    ggplot2::stat_density_2d_filled(
      bins = bins,
      alpha = alpha,
      contour_var = "ndensity"
    ) +
    viridis::scale_fill_viridis(discrete = TRUE, option = "turbo") +
    ggplot2::stat_density_2d(
      color = "black",
      bins = levels,
      contour_var = "ndensity"
    ) +
    ggplot2::labs(
      x = factor_names[1],
      y = factor_names[2],
      title = "Posterior Distribution of the Optimum",
      fill = "Density"
    ) +
    ggplot2::theme_minimal()

  # Add sample points (if requested) before mean optimum marker
  if (show_points) {
    p <- p + ggplot2::geom_point(
      alpha = point_alpha,
      size = 1
    )
  }

  # Add mean optimum marker last so it's always visible, use [[ ]] for safety
  p <- p + ggplot2::geom_point(
    data = data.frame(
      x = opt_mean[[factor_names[1]]],
      y = opt_mean[[factor_names[2]]]
    ),
    ggplot2::aes(x = .data$x, y = .data$y),
    color = "red",
    size = 3,
    shape = 4,
    stroke = 1.5
  )

  return(p)
}