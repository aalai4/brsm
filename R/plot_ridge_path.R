plot_ridge_path <- function(ridge_draws,
                            factor_names,
                            response = NULL,
                            show_draws = FALSE,
                            alpha = 0.2,
                            mean_color = "red",
                            path_color = "black") {
  .brsm_require_ggplot2()

  factor_names <- .brsm_validate_factor_names(
    factor_names,
    require_two = TRUE,
    two_factors_message =
      "plot_ridge_path currently supports exactly two factors."
  )
  ridge_draws <- .brsm_validate_draws(ridge_draws)

  .brsm_check_columns(
    factor_names, ridge_draws,
    "factor_names must correspond to columns in ridge_draws."
  )

  if (is.null(response) || !response %in% colnames(ridge_draws)) {
    stop("response must be supplied and correspond to a column in ridge_draws.")
  }
  if (!is.numeric(ridge_draws[[response]])) {
    stop("response column must be numeric.")
  }

  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
  }

  # Remove missing values
  removed <- sum(!stats::complete.cases(ridge_draws))
  ridge_draws <- ridge_draws[stats::complete.cases(ridge_draws), , drop = FALSE]
  if (removed > 0) {
    warning(removed, " rows removed from ridge_draws due to missing values.")
  }

  if (nrow(ridge_draws) == 0) {
    stop("ridge_draws contains no valid rows.")
  }

  # Posterior mean ridge path
  mean_path <- stats::aggregate(
    ridge_draws[, factor_names],
    by = list(level = ridge_draws[[response]]),
    FUN = mean
  )
  colnames(mean_path) <- c("level", factor_names)
  mean_path <- mean_path[order(mean_path$level), ] # Ensure correct path order

  # Base plot
  p <- ggplot2::ggplot(
    ridge_draws,
    ggplot2::aes(
      x = .data[[factor_names[1]]],
      y = .data[[factor_names[2]]]
    )
  ) +
    ggplot2::labs(
      x = factor_names[1],
      y = factor_names[2],
      title = "Posterior Mean Ridge Optimization Path",
      subtitle = paste("Indexed by", response)
    ) +
    ggplot2::theme_minimal()

  # Optional posterior draws
  if (show_draws) {
    p <- p +
      ggplot2::geom_point(
        alpha = alpha,
        size = 0.7
      )
  }

  # Mean ridge path
  p <- p +
    ggplot2::geom_path(
      data = mean_path,
      ggplot2::aes(
        x = .data[[factor_names[1]]],
        y = .data[[factor_names[2]]]
      ),
      color = path_color,
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = mean_path,
      ggplot2::aes(
        x = .data[[factor_names[1]]],
        y = .data[[factor_names[2]]]
      ),
      color = mean_color,
      size = 2
    )

  return(p)
}