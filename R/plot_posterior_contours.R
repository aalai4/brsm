#' Plot Posterior Contours
#'
#' @param draws Validated draw data frame.
#' @param factor_names Character vector of factor names.
#' @param ranges Named list of factor ranges.
#' @param n Grid resolution.
#' @param type Surface type ("mean", "uncertainty", "quantile").
#' @param probs Probability quantiles for uncertainty.
#' @param quantile Quantile value if type="quantile".
#' @param bins Number of contour bins.
#' @param vary_factors Optional 2-factor pair for conditional plotting.
#' @param fixed Optional named numeric vector for fixed factor values.
#' @param conditioning Strategy for fixed factors.
#' @param slice Optional named list for sliced factor values.
#' @param pairwise Logical; if TRUE, generate all factor-pair panels.
#' @param overlay_stationary Logical; if TRUE, overlay stationary points.
#' @param stationary_draws Optional precomputed stationary draws.
#' @param overlay_type Type of stationary overlay ("mean" or "posterior").
#' @param overlay_alpha Transparency for posterior draw overlay.
#' @param overlay_max_draws Cap on posterior draws for overlay.
#' @param seed Random seed.
#'
#' @return A ggplot2 object.
#'
#' @details
#' Contour fills use the viridis "turbo" palette, which is perceptually uniform
#' and colorblind-safe. Stationary points are overlaid in red (single-pair, non-sliced
#' contours only). For 3+ factors, use `pairwise=TRUE` for all pairs or `slice`
#' for conditional slices.
#' @keywords internal
plot_posterior_contours <- function(draws,
                                    factor_names,
                                    ranges,
                                    n = 60,
                                    type = c("mean", "uncertainty", "quantile"),
                                    probs = c(0.025, 0.975),
                                    quantile = 0.5,
                                    bins = 10,
                                    vary_factors = NULL,
                                    fixed = NULL,
                                    conditioning = c(
                                      "center", "optimum_mean", "user"
                                    ),
                                    slice = NULL,
                                    pairwise = FALSE,
                                    overlay_stationary = FALSE,
                                    stationary_draws = NULL,
                                    overlay_type = c("mean", "posterior"),
                                    overlay_alpha = 0.1,
                                    overlay_max_draws = 2000,
                                    seed = NULL) {
  draws <- .brsm_validate_draws(draws)
  factor_names <- .brsm_validate_factor_names(factor_names)

  if (length(factor_names) < 2) {
    stop("plot_posterior_contours requires at least two factors.")
  }

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

  conditioning <- match.arg(conditioning)

  if (!is.null(vary_factors)) {
    vary_factors <- .brsm_validate_factor_names(
      vary_factors,
      require_two = TRUE,
      two_factors_message = "vary_factors must contain exactly two factors."
    )
    if (!all(vary_factors %in% factor_names)) {
      stop("vary_factors must be a subset of factor_names.")
    }
  }

  if (!is.logical(pairwise) || length(pairwise) != 1L) {
    stop("pairwise must be a single logical value.")
  }

  if (!is.null(slice)) {
    if (!is.list(slice) || length(slice) != 1L || is.null(names(slice))) {
      stop("slice must be a named list with exactly one sliced factor.")
    }
    slice_var <- names(slice)[1]
    slice_vals <- slice[[1]]
    if (!slice_var %in% factor_names) {
      stop("slice factor must be present in factor_names.")
    }
    if (!is.numeric(slice_vals) || length(slice_vals) == 0L ||
      any(!is.finite(slice_vals))) {
      stop("slice values must be a non-empty finite numeric vector.")
    }
  } else {
    slice_var <- NULL
    slice_vals <- NULL
  }

  resolve_fixed <- function(vary_pair) {
    non_vary <- setdiff(factor_names, vary_pair)
    if (length(non_vary) == 0L) {
      return(stats::setNames(numeric(0), character(0)))
    }

    if (conditioning == "user") {
      if (is.null(fixed) || !is.numeric(fixed) || is.null(names(fixed))) {
        stop(
          "When conditioning='user', fixed must be a named numeric ",
          "vector covering non-varied factors."
        )
      }
      if (!all(non_vary %in% names(fixed))) {
        stop(
          "fixed is missing values for non-varied factors: ",
          paste(setdiff(non_vary, names(fixed)), collapse = ", "),
          "."
        )
      }
      return(as.numeric(fixed[non_vary]))
    }

    center_vals <- vapply(non_vary, function(f) {
      mean(sort(as.numeric(ranges[[f]])))
    }, numeric(1))
    names(center_vals) <- non_vary

    if (conditioning == "center") {
      return(center_vals)
    }

    stationary_all <- tryCatch(
      {
        sp <- stationary_point(draws, factor_names = factor_names)
        sp <- as.data.frame(sp)
        colnames(sp) <- factor_names
        sp
      },
      error = function(e) NULL
    )
    if (is.null(stationary_all) || nrow(stationary_all) == 0) {
      warning(
        "Could not compute stationary points for conditioning='optimum_mean'; ",
        "falling back to center conditioning."
      )
      return(center_vals)
    }

    opt_vals <- vapply(non_vary, function(f) {
      x <- stationary_all[[f]]
      if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    }, numeric(1))
    names(opt_vals) <- non_vary

    bad <- is.na(opt_vals)
    if (any(bad)) {
      warning(
        "Optimum mean unavailable for some fixed factors; using center for: ",
        paste(names(opt_vals)[bad], collapse = ", "),
        "."
      )
      opt_vals[bad] <- center_vals[bad]
    }

    opt_vals
  }

  pair_list <- if (isTRUE(pairwise)) {
    utils::combn(factor_names, 2, simplify = FALSE)
  } else {
    list(if (is.null(vary_factors)) factor_names[1:2] else vary_factors)
  }

  if (!is.null(slice_var) && !isTRUE(pairwise)) {
    current_pair <- pair_list[[1]]
    if (slice_var %in% current_pair) {
      stop("slice factor must differ from vary_factors when pairwise=FALSE.")
    }
  }

  if (!is.null(stationary_draws)) {
    stationary_draws <- as.data.frame(stationary_draws)
    .brsm_check_columns(factor_names, stationary_draws,
      "stationary_draws must contain columns for factor_names.")
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

  # Generate conditional plotting grid(s)
  grid_list <- list()
  for (pair in pair_list) {
    pair <- as.character(pair)
    pair_label <- paste(pair, collapse = " vs ")

    fixed_vals <- resolve_fixed(pair)

    pair_ranges <- stats::setNames(lapply(pair, function(f) ranges[[f]]), pair)
    pair_grid <- surface_grid(ranges = pair_ranges, n = n)

    if (length(fixed_vals) > 0) {
      for (nm in names(fixed_vals)) {
        pair_grid[[nm]] <- fixed_vals[[nm]]
      }
    }

    pair_grid <- pair_grid[, factor_names, drop = FALSE]

    if (is.null(slice_var)) {
      pair_grid$.slice_label <- "base"
      pair_grid$.slice_value <- NA_real_
      pair_grid$.panel <- pair_label
      pair_grid$.pair1 <- pair[1]
      pair_grid$.pair2 <- pair[2]
      grid_list[[length(grid_list) + 1L]] <- pair_grid
    } else {
      for (sv in slice_vals) {
        slice_grid <- pair_grid
        slice_grid[[slice_var]] <- sv
        slice_grid$.slice_label <- paste0(slice_var, "=", signif(sv, 6))
        slice_grid$.slice_value <- sv
        slice_grid$.panel <- paste(pair_label, "|", slice_grid$.slice_label[1])
        slice_grid$.pair1 <- pair[1]
        slice_grid$.pair2 <- pair[2]
        grid_list[[length(grid_list) + 1L]] <- slice_grid
      }
    }
  }

  grid <- do.call(rbind, grid_list)
  rownames(grid) <- NULL

  # Compute predictions
  pred <- predict_surface(
    draws = draws,
    factor_names = factor_names,
    newdata = grid,
    summary = TRUE,
    probs = pred_probs
  )

  pred$.panel <- grid$.panel
  pred$.pair1 <- grid$.pair1
  pred$.pair2 <- grid$.pair2

  pred$.x <- vapply(seq_len(nrow(pred)), function(i) {
    as.numeric(pred[[pred$.pair1[i]]][i])
  }, numeric(1))
  pred$.y <- vapply(seq_len(nrow(pred)), function(i) {
    as.numeric(pred[[pred$.pair2[i]]][i])
  }, numeric(1))

  # Determine which variable to plot
  if (type == "mean") {
    pred$z <- pred$mean
    plot_title <- "Posterior Mean Surface"
  } else if (type == "uncertainty") {
    lower_name <- paste0("q", formatC(probs[1] * 100, format = "f", digits = 1))
    upper_name <- paste0("q", formatC(probs[2] * 100, format = "f", digits = 1))
    if (!all(c(lower_name, upper_name) %in% names(pred))) {
      stop("predict_surface() did not return expected quantile columns.")
    }
    pred$z <- pred[[upper_name]] - pred[[lower_name]]
    plot_title <- "Posterior Uncertainty Surface"
  } else if (type == "quantile") {
    qname <- paste0("q", formatC(quantile * 100, format = "f", digits = 1))
    if (!(qname %in% names(pred))) {
      stop("requested quantile not found in predictions.")
    }
    pred$z <- pred[[qname]]
    plot_title <- paste0("Posterior Quantile Surface (", quantile, ")")
  }

  p <- ggplot2::ggplot(
    pred,
    ggplot2::aes(
      x = .data$.x,
      y = .data$.y
    )
  ) +
    ggplot2::geom_contour_filled(
      ggplot2::aes(z = .data$z),
      bins = bins
    ) +
    viridis::scale_fill_viridis(discrete = TRUE, option = "turbo") +
    ggplot2::labs(
      x = "factor 1",
      y = "factor 2",
      title = plot_title
    ) +
    ggplot2::theme_minimal()

  if (length(unique(pred$.panel)) > 1L) {
    p <- p + ggplot2::facet_wrap(~.panel, scales = "free")
  }

  can_overlay_stationary <- overlay_stationary &&
    !isTRUE(pairwise) &&
    is.null(slice_var)

  if (can_overlay_stationary) {
    if (!is.numeric(overlay_alpha) || overlay_alpha < 0 || overlay_alpha > 1) {
      stop("overlay_alpha must be between 0 and 1.")
    }
    if (!is.numeric(overlay_max_draws) || overlay_max_draws <= 0) {
      stop("overlay_max_draws must be a positive integer.")
    }
    overlay_max_draws <- as.integer(overlay_max_draws)
    overlay_type <- match.arg(overlay_type)

    pair <- pair_list[[1]]
    if (is.null(stationary_draws)) {
      stationary <- stationary_point(draws, factor_names)
      stationary <- as.data.frame(stationary)
      colnames(stationary) <- factor_names
    } else {
      stationary <- as.data.frame(stationary_draws)
    }
    stationary <- stationary[, pair, drop = FALSE]
    # Remove all-NA rows before overlay
    stationary <- stationary[stats::complete.cases(stationary), , drop = FALSE]
    if (nrow(stationary) == 0) {
      warning("No valid stationary points available for overlay.")
    } else if (overlay_type == "mean") {
      stationary_mean <- as.data.frame(t(colMeans(stationary, na.rm = TRUE)))
      colnames(stationary_mean) <- pair
      p <- p + ggplot2::geom_point(
        data = stationary_mean,
        ggplot2::aes(
          x = .data[[pair[1]]],
          y = .data[[pair[2]]]
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
          x = .data[[pair[1]]],
          y = .data[[pair[2]]]
        ),
        color = "red", alpha = overlay_alpha, size = 2, shape = 16
      )
    }
  } else if (isTRUE(overlay_stationary)) {
    warning(
      "Stationary overlays are only available for single-pair, non-sliced ",
      "contours."
    )
  }

  return(p)
}