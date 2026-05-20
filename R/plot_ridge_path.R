plot_ridge_path <- function(ridge_draws,
                            factor_names,
                            response = NULL,
                            show_draws = FALSE,
                            alpha = 0.2,
                            show_interval_crossbars = TRUE,
                            interval_prob = 0.9,
                            interval_alpha = 0.35,
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

  if (!is.logical(show_interval_crossbars) ||
      length(show_interval_crossbars) != 1L) {
    stop("show_interval_crossbars must be a single logical value.")
  }

  if (!is.numeric(interval_prob) || length(interval_prob) != 1L ||
      !is.finite(interval_prob) || interval_prob <= 0 || interval_prob >= 1) {
    stop("interval_prob must be a finite numeric scalar in (0, 1).")
  }

  if (!is.numeric(interval_alpha) || length(interval_alpha) != 1L ||
      interval_alpha < 0 || interval_alpha > 1) {
    stop("interval_alpha must be between 0 and 1.")
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

  interval_probs <- c((1 - interval_prob) / 2, 1 - (1 - interval_prob) / 2)
  interval_path <- do.call(
    rbind,
    lapply(split(ridge_draws, ridge_draws[[response]]), function(df_level) {
      data.frame(
        level = df_level[[response]][1],
        x = mean(df_level[[factor_names[1]]], na.rm = TRUE),
        y = mean(df_level[[factor_names[2]]], na.rm = TRUE),
        x_low = stats::quantile(df_level[[factor_names[1]]], interval_probs[1], na.rm = TRUE),
        x_high = stats::quantile(df_level[[factor_names[1]]], interval_probs[2], na.rm = TRUE),
        y_low = stats::quantile(df_level[[factor_names[2]]], interval_probs[1], na.rm = TRUE),
        y_high = stats::quantile(df_level[[factor_names[2]]], interval_probs[2], na.rm = TRUE),
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    })
  )
  interval_path <- interval_path[order(interval_path$level), , drop = FALSE]

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
      subtitle = if (isTRUE(show_interval_crossbars)) {
        paste0(
          "Indexed by ", response,
          "; crossbars show ",
          formatC(100 * interval_prob, format = "f", digits = 0),
          "% pointwise posterior intervals"
        )
      } else {
        paste("Indexed by", response)
      }
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

  if (isTRUE(show_interval_crossbars)) {
    p <- p +
      ggplot2::geom_segment(
        data = interval_path,
        ggplot2::aes(x = .data$x_low, xend = .data$x_high, y = .data$y, yend = .data$y),
        color = path_color,
        alpha = interval_alpha,
        linewidth = 0.45
      ) +
      ggplot2::geom_segment(
        data = interval_path,
        ggplot2::aes(x = .data$x, xend = .data$x, y = .data$y_low, yend = .data$y_high),
        color = path_color,
        alpha = interval_alpha,
        linewidth = 0.45
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