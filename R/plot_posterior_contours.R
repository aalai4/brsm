plot_posterior_contours <- function(draws,
                                    factor_names,
                                    ranges,
                                    n = 60,
                                    type = c("mean", "uncertainty", "quantile"),
                                    probs = c(0.025, 0.975),
                                    quantile = 0.5,
                                    bins = 10,
                                    overlay_stationary = FALSE,
                                    stationary_draws = NULL,
                                    overlay_type = c("mean", "posterior"),
                                    overlay_alpha = 0.1,
                                    overlay_max_draws = 2000,
                                    seed = NULL) {
  draws <- .brsm_validate_draws(draws)
  factor_names <- .brsm_validate_factor_names(
    factor_names,
    require_two = TRUE,
    two_factors_message =
      "plot_posterior_contours currently supports exactly two factors."
  )

  if (is.null(ranges) || !is.list(ranges)) {
    stop("ranges must be a named list with ranges for each factor.")
  }
  .brsm_check_columns(
    factor_names, ranges,
    "ranges must contain entries for each factor in factor_names."
  )
  for (nm in factor_names) {
    r <- ranges[[nm]]
    if (!is.numeric(r) || length(r) != 2) {
      stop("each entry in ranges must be a numeric vector of length 2.")
    }
  }

  if (!is.numeric(n) || n <= 0) {
    stop("n must be a positive integer.")
  }
  n <- as.integer(n)

  if (!is.numeric(quantile) || quantile < 0 || quantile > 1) {
    stop("quantile must be between 0 and 1.")
  }

  if (!is.numeric(bins) || bins <= 0) {
    stop("bins must be a positive integer.")
  }
  bins <- as.integer(bins)

  .brsm_require_ggplot2()

  if (!is.null(stationary_draws)) {
    stationary_draws <- as.data.frame(stationary_draws)
    .brsm_check_columns(
      factor_names, stationary_draws,
      "stationary_draws must contain columns for factor_names."
    )
  }

  type <- match.arg(type)

  if (type == "uncertainty") {
    probs <- .brsm_validate_probs(
      probs,
      require_length = 2,
      message = paste0(
        "probs must be a numeric vector of length 2 ",
        "with values between 0 and 1."
      )
    )
    pred_probs <- probs
  } else if (type == "quantile") {
    probs <- .brsm_validate_probs(probs)
    pred_probs <- sort(unique(c(probs, quantile)))
  } else {
    probs <- .brsm_validate_probs(probs)
    pred_probs <- probs
  }

  # Generate full surface grid
  grid <- surface_grid(ranges = ranges, n = n)

  # Compute predictions
  pred <- predict_surface(
    draws = draws,
    factor_names = factor_names,
    newdata = grid,
    summary = TRUE,
    probs = pred_probs
  )

  # Determine which variable to plot
  if (type == "mean") {
    zvar <- "mean"
    pred$z <- pred$mean
    plot_title <- "Posterior Mean Surface"
  } else if (type == "uncertainty") {
    lower_name <- paste0("q", formatC(probs[1] * 100, format = "f", digits = 1))
    upper_name <- paste0("q", formatC(probs[2] * 100, format = "f", digits = 1))
    if (!all(c(lower_name, upper_name) %in% names(pred))) {
      stop("predict_surface() did not return expected quantile columns.")
    }
    pred$z <- pred[[upper_name]] - pred[[lower_name]]
    zvar <- "uncertainty"
    plot_title <- "Posterior Uncertainty Surface"
  } else if (type == "quantile") {
    qname <- paste0("q", formatC(quantile * 100, format = "f", digits = 1))
    if (!(qname %in% names(pred))) {
      stop("requested quantile not found in predictions.")
    }
    pred$z <- pred[[qname]]
    zvar <- qname
    plot_title <- paste0("Posterior Quantile Surface (", quantile, ")")
  }

  p <- ggplot2::ggplot(
    pred,
    ggplot2::aes(
      x = .data[[factor_names[1]]],
      y = .data[[factor_names[2]]]
    )
  ) +
    ggplot2::geom_contour_filled(
      ggplot2::aes(z = .data$z),
      bins = bins
    ) +
    ggplot2::labs(
      x = factor_names[1],
      y = factor_names[2],
      title = plot_title
    ) +
    ggplot2::theme_minimal()

  if (overlay_stationary) {
    if (!is.numeric(overlay_alpha) || overlay_alpha < 0 || overlay_alpha > 1) {
      stop("overlay_alpha must be between 0 and 1.")
    }
    if (!is.numeric(overlay_max_draws) || overlay_max_draws <= 0) {
      stop("overlay_max_draws must be a positive integer.")
    }
    overlay_max_draws <- as.integer(overlay_max_draws)
    overlay_type <- match.arg(overlay_type)
    if (is.null(stationary_draws)) {
      stationary <- stationary_point(draws, factor_names)
      stationary <- as.data.frame(stationary)
      colnames(stationary) <- factor_names
    } else {
      stationary <- as.data.frame(stationary_draws)
    }
    stationary <- stationary[, factor_names, drop = FALSE]
    # Remove all-NA rows before overlay
    stationary <- stationary[stats::complete.cases(stationary), , drop = FALSE]
    if (nrow(stationary) == 0) {
      warning("No valid stationary points available for overlay.")
    } else if (overlay_type == "mean") {
      stationary_mean <- as.data.frame(t(colMeans(stationary, na.rm = TRUE)))
      colnames(stationary_mean) <- factor_names
      p <- p + ggplot2::geom_point(
        data = stationary_mean,
        ggplot2::aes(
          x = .data[[factor_names[1]]],
          y = .data[[factor_names[2]]]
        ),
        color = "red", size = 3, shape = 4, stroke = 2
      )
    } else if (overlay_type == "posterior") {
      n_stat <- nrow(stationary)
      if (n_stat > overlay_max_draws) {
        if (!is.null(seed)) set.seed(seed)
        stationary <- stationary[
          sample.int(n_stat, overlay_max_draws), ,
          drop = FALSE
        ]
      }
      p <- p + ggplot2::geom_point(
        data = stationary,
        ggplot2::aes(
          x = .data[[factor_names[1]]],
          y = .data[[factor_names[2]]]
        ),
        color = "red", alpha = overlay_alpha, size = 2, shape = 16
      )
    }
  }

  return(p)
}