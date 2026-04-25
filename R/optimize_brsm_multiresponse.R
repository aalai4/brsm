#' Multi-Response Optimization via Posterior Desirability
#'
#' Combines response-specific desirability functions over posterior predictions
#' at candidate points. Predictions are obtained with
#' [posterior_predict_brsm()], mapped to \eqn{[0, 1]}, and aggregated with a
#' weighted geometric mean.
#'
#' @param models Named list of response models. Each element can be a
#'   \code{brsm_fit}, \code{brmsfit}, or posterior-draw data frame accepted by
#'   [posterior_predict_brsm()].
#' @param desirability_specs Named list of desirability specifications, one per
#'   model. Each element must include \code{goal} and bounds:
#'   \itemize{
#'   \item \code{goal = "max"} or \code{"maximize"}: requires \code{low},
#'   \code{high}
#'   \item \code{goal = "min"} or \code{"minimize"}: requires \code{low},
#'   \code{high}
#'   \item \code{goal = "target"}: requires \code{low}, \code{target},
#'   \code{high}
#'   }
#'   Optional fields: \code{weight} (shape, default 1),
#'   \code{importance} (combination weight, default 1).
#' @param factor_names Character vector of factor names.
#' @param candidate_points Optional data frame of candidate points.
#' @param ranges Optional named list of factor ranges used to generate a grid
#'   when \code{candidate_points} is \code{NULL}.
#' @param n_grid Grid resolution per factor when generating candidates from
#'   \code{ranges}.
#' @param include_residual Logical; forwarded to [posterior_predict_brsm()].
#' @param sigma Optional residual SD override forwarded to
#'   [posterior_predict_brsm()].
#' @param draw_subset Optional draw subset forwarded to
#'   [posterior_predict_brsm()].
#' @param max_draws Optional draw cap forwarded to [posterior_predict_brsm()].
#' @param probs Probabilities for desirability summaries across draws.
#' @param optimize_metric Which summary metric to maximize when selecting the
#'   best point. One of \code{"mean"} or \code{"median"}.
#' @param return_draws Logical; if \code{TRUE}, include per-draw combined
#'   desirability matrix in output.
#' @param seed Optional random seed (used when \code{include_residual = TRUE}).
#'
#' @return A list with components:
#'   \code{candidate_points} (with desirability summaries),
#'   \code{best_point} (single-row data frame), \code{response_at_best}
#'   (per-response summary at the selected point), \code{draw_count},
#'   \code{desirability_specs}, and optional
#'   \code{combined_desirability_draws}.
#' @export
optimize_brsm_multiresponse <- function(models,
                                        desirability_specs,
                                        factor_names,
                                        candidate_points = NULL,
                                        ranges = NULL,
                                        n_grid = 25,
                                        include_residual = FALSE,
                                        sigma = NULL,
                                        draw_subset = NULL,
                                        max_draws = NULL,
                                        probs = c(0.025, 0.5, 0.975),
                                        optimize_metric = c("mean", "median"),
                                        return_draws = FALSE,
                                        seed = NULL) {
  if (!is.list(models) || length(models) < 2L) {
    stop("models must be a list with at least two response models.")
  }
  if (is.null(names(models)) || any(names(models) == "")) {
    stop("models must be a named list.")
  }

  factor_names <- .brsm_validate_factor_names(factor_names)
  optimize_metric <- match.arg(optimize_metric)

  if (!is.logical(return_draws) || length(return_draws) != 1L ||
      is.na(return_draws)) {
    stop("return_draws must be TRUE or FALSE.")
  }

  probs <- .brsm_validate_probs(
    probs,
    require_finite = TRUE,
    message = "probs must be finite values between 0 and 1."
  )

  specs <- .brsm_validate_desirability_specs(
    desirability_specs = desirability_specs,
    model_names = names(models)
  )

  if (is.null(candidate_points)) {
    if (is.null(ranges)) {
      stop("Provide candidate_points or ranges.")
    }
    candidate_points <- surface_grid(ranges = ranges, n = n_grid)
  }

  if (!is.data.frame(candidate_points)) {
    stop("candidate_points must be a data.frame.")
  }

  .brsm_check_columns(
    factor_names,
    candidate_points,
    "candidate_points must contain columns for all factor_names."
  )

  if (nrow(candidate_points) == 0L) {
    stop("candidate_points must contain at least one row.")
  }

  pred_mats <- vector("list", length(models))
  names(pred_mats) <- names(models)

  for (i in seq_along(models)) {
    model_name <- names(models)[i]
    seed_i <- if (is.null(seed)) NULL else (as.integer(seed) + i - 1L)

    pred_mats[[i]] <- posterior_predict_brsm(
      object = models[[i]],
      factor_names = factor_names,
      newdata = candidate_points,
      include_residual = include_residual,
      sigma = sigma,
      summary = FALSE,
      draw_subset = draw_subset,
      max_draws = max_draws,
      return_matrix = TRUE,
      seed = seed_i
    )

    if (!is.matrix(pred_mats[[i]])) {
      stop("Internal error: posterior predictions must return a matrix.")
    }

    if (!identical(ncol(pred_mats[[i]]), nrow(candidate_points))) {
      stop("Internal error: prediction matrix column count mismatch.")
    }

    spec <- specs[[model_name]]
    pred_mats[[i]] <- .brsm_desirability_transform(
      y = pred_mats[[i]],
      goal = spec$goal,
      low = spec$low,
      high = spec$high,
      target = spec$target,
      weight = spec$weight
    )
  }

  draw_counts <- vapply(pred_mats, nrow, integer(1))
  n_draws <- min(draw_counts)
  if (any(draw_counts != n_draws)) {
    warning(
      "Responses had differing draw counts; truncating all to ",
      n_draws,
      " draw(s)."
    )
    pred_mats <- lapply(pred_mats, function(m) m[seq_len(n_draws), , drop = FALSE])
  }

  importances <- vapply(specs[names(models)], function(s) s$importance, numeric(1))
  wsum <- sum(importances)

  combined_log <- matrix(0, nrow = n_draws, ncol = nrow(candidate_points))
  any_zero <- matrix(FALSE, nrow = n_draws, ncol = nrow(candidate_points))
  for (i in seq_along(pred_mats)) {
    d <- pred_mats[[i]]
    zeros <- d <= 0
    any_zero <- any_zero | zeros
    d_safe <- pmax(d, .Machine$double.xmin)
    combined_log <- combined_log + importances[i] * log(d_safe)
  }
  combined <- exp(combined_log / wsum)
  combined[any_zero] <- 0

  desirability_summary <- apply(combined, 2, function(x) {
    qs <- stats::quantile(x, probs = probs, names = FALSE, na.rm = TRUE)
    names(qs) <- paste0("q", formatC(probs * 100, format = "f", digits = 1))
    c(mean = mean(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE), qs)
  })
  desirability_summary <- t(desirability_summary)

  candidates_out <- cbind(candidate_points, as.data.frame(desirability_summary))
  candidates_out$.point_id <- seq_len(nrow(candidates_out))

  metric_col <- optimize_metric
  best_idx <- which.max(candidates_out[[metric_col]])
  best_point <- candidates_out[best_idx, , drop = FALSE]

  response_at_best <- do.call(rbind, lapply(seq_along(models), function(i) {
    y <- pred_mats[[i]][, best_idx]
    qs <- stats::quantile(y, probs = probs, names = FALSE, na.rm = TRUE)
    data.frame(
      response = names(models)[i],
      mean = mean(y, na.rm = TRUE),
      median = stats::median(y, na.rm = TRUE),
      t(as.data.frame(qs)),
      stringsAsFactors = FALSE
    )
  }))

  q_names <- paste0("q", formatC(probs * 100, format = "f", digits = 1))
  names(response_at_best)[seq.int(from = 4, length.out = length(q_names))] <- q_names

  out <- list(
    candidate_points = candidates_out,
    best_point = best_point,
    response_at_best = response_at_best,
    draw_count = n_draws,
    desirability_specs = specs
  )

  if (isTRUE(return_draws)) {
    out$combined_desirability_draws <- combined
  }

  out
}


.brsm_validate_desirability_specs <- function(desirability_specs, model_names) {
  if (!is.list(desirability_specs) || length(desirability_specs) == 0L) {
    stop("desirability_specs must be a non-empty named list.")
  }
  if (is.null(names(desirability_specs)) || any(names(desirability_specs) == "")) {
    stop("desirability_specs must be a named list.")
  }

  missing_specs <- setdiff(model_names, names(desirability_specs))
  if (length(missing_specs) > 0L) {
    stop(
      "Missing desirability specs for model(s): ",
      paste(missing_specs, collapse = ", "),
      "."
    )
  }

  specs <- desirability_specs[model_names]

  for (nm in names(specs)) {
    s <- specs[[nm]]
    if (!is.list(s)) {
      stop("Each desirability spec must be a list.")
    }

    if (is.null(s$goal)) {
      stop("Each desirability spec must include goal.")
    }

    goal <- tolower(as.character(s$goal[[1L]]))
    if (goal %in% c("max", "maximize")) goal <- "maximize"
    if (goal %in% c("min", "minimize")) goal <- "minimize"

    if (!goal %in% c("maximize", "minimize", "target")) {
      stop("goal must be one of 'max'/'maximize', 'min'/'minimize', or 'target'.")
    }

    if (is.null(s$low) || is.null(s$high) ||
          !is.numeric(s$low) || !is.numeric(s$high) ||
          length(s$low) != 1L || length(s$high) != 1L ||
          !is.finite(s$low) || !is.finite(s$high) ||
          s$low >= s$high) {
      stop("Each desirability spec must include finite low < high.")
    }

    target <- NA_real_
    if (goal == "target") {
      if (is.null(s$target) || !is.numeric(s$target) || length(s$target) != 1L ||
            !is.finite(s$target)) {
        stop("Target desirability requires finite target.")
      }
      target <- as.numeric(s$target)
      if (target <= s$low || target >= s$high) {
        stop("For goal='target', target must satisfy low < target < high.")
      }
    }

    weight <- if (is.null(s$weight)) 1 else as.numeric(s$weight)
    importance <- if (is.null(s$importance)) 1 else as.numeric(s$importance)

    if (!is.numeric(weight) || length(weight) != 1L || !is.finite(weight) ||
          weight <= 0) {
      stop("weight must be a finite positive scalar.")
    }
    if (!is.numeric(importance) || length(importance) != 1L ||
          !is.finite(importance) || importance <= 0) {
      stop("importance must be a finite positive scalar.")
    }

    specs[[nm]] <- list(
      goal = goal,
      low = as.numeric(s$low),
      high = as.numeric(s$high),
      target = target,
      weight = weight,
      importance = importance
    )
  }

  specs
}


.brsm_desirability_transform <- function(y,
                                         goal,
                                         low,
                                         high,
                                         target = NA_real_,
                                         weight = 1) {
  d <- matrix(0, nrow = nrow(y), ncol = ncol(y))

  if (goal == "maximize") {
    idx_hi <- y >= high
    idx_mid <- y > low & y < high
    d[idx_hi] <- 1
    d[idx_mid] <- ((y[idx_mid] - low) / (high - low))^weight
    return(d)
  }

  if (goal == "minimize") {
    idx_lo <- y <= low
    idx_mid <- y > low & y < high
    d[idx_lo] <- 1
    d[idx_mid] <- ((high - y[idx_mid]) / (high - low))^weight
    return(d)
  }

  # target
  idx_left <- y >= low & y < target
  idx_right <- y > target & y <= high
  idx_target <- y == target

  d[idx_target] <- 1
  d[idx_left] <- ((y[idx_left] - low) / (target - low))^weight
  d[idx_right] <- ((high - y[idx_right]) / (high - target))^weight

  d
}