#' Compute Steepest Ascent Path on Response Surface
#'
#' Computes the path of steepest ascent (or descent) on the response surface
#' starting from a specified point. Useful for sequential design of experiments
#' to locate the vicinity of the optimum.
#'
#' @param object A `brsm_fit` object, `brmsfit` object, or
#'   data frame of posterior draws with Bayesian coefficient columns.
#' @param factor_names Character vector of factor names
#'   (if object is data frame).
#' @param start Starting point for ascent path
#'   (default: center of coded region).
#' @param step_size Step size along gradient direction.
#' @param n_steps Number of steps to compute.
#' @param tol Numerical tolerance for gradient computations.
#' @param return_mean_path If TRUE, return mean path;
#'   if FALSE, return all draws.
#' @param drop_null_paths Drop paths that don't increase response.
#'
#' @return Data frame with steepest ascent path coordinates and response values.
#'
#' @export
steepest_ascent <- function(
    object,
    factor_names = NULL,
    start = NULL,
    step_size = 1,
    n_steps = 10,
    tol = 1e-8,
    return_mean_path = TRUE,
    drop_null_paths = TRUE) {
  UseMethod("steepest_ascent")
}

#' @rdname steepest_ascent
#' @export
steepest_ascent.brsm_fit <- function(
    object,
    factor_names = NULL,
    start = NULL,
    step_size = 1,
    n_steps = 10,
    tol = 1e-8,
    return_mean_path = TRUE,
    drop_null_paths = TRUE) {
  if (is.null(factor_names)) {
    factor_names <- object$factor_names
  }
  if (is.null(start)) {
    start <- stats::setNames(rep(0, length(factor_names)), factor_names)
  }
  steepest_ascent.default(
    as_brsm_draws(object),
    factor_names = factor_names,
    start = start,
    step_size = step_size,
    n_steps = n_steps,
    tol = tol,
    return_mean_path = return_mean_path,
    drop_null_paths = drop_null_paths
  )
}

#' @rdname steepest_ascent
#' @export
steepest_ascent.default <- function(
    object,
    factor_names = NULL,
    start = NULL,
    step_size = 1,
    n_steps = 10,
    tol = 1e-8,
    return_mean_path = TRUE,
    drop_null_paths = TRUE) {
  draws <- .brsm_validate_draws(object)
  if (is.null(factor_names)) {
    stop("factor_names must be supplied if object does not contain metadata.")
  }
  factor_names <- .brsm_validate_factor_names(factor_names)

  if (!is.numeric(start)) {
    stop("start must be a numeric vector.")
  }

  if (is.null(names(start))) {
    stop("start must be a named numeric vector.")
  }

  if (!all(factor_names %in% names(start))) {
    stop("start must contain values for all factor_names.")
  }

  if (!is.numeric(step_size) || length(step_size) != 1 ||
    !is.finite(step_size) || step_size <= 0) {
    stop("step_size must be a positive, finite numeric value.")
  }

  if (!is.numeric(n_steps) || length(n_steps) != 1 || n_steps < 1) {
    stop("n_steps must be a positive integer.")
  }
  n_steps <- as.integer(n_steps)

  if (!is.numeric(tol) || length(tol) != 1 || tol < 0) {
    stop("tol must be a non-negative numeric scalar.")
  }

  n_draws <- nrow(draws)
  n_factors <- length(factor_names)

  # Storage
  direction_matrix <- matrix(NA_real_, nrow = n_draws, ncol = n_factors)
  colnames(direction_matrix) <- factor_names

  path_list <- vector("list", n_draws)
  start_vec <- as.numeric(start[factor_names])

  for (d in seq_len(n_draws)) {
    grad <- gradient_quadratic(
      draws = draws[d, , drop = FALSE],
      factor_names = factor_names,
      x = start_vec
    )

    grad_vec <- as.numeric(grad[1, paste0("d_d", factor_names), drop = TRUE])
    norm <- sqrt(sum(grad_vec^2))

    if (norm < tol || is.na(norm)) {
      next
    }

    direction <- grad_vec / norm
    direction_matrix[d, ] <- direction

    # Build ascent path in one vectorized operation.
    step_seq <- seq.int(0, n_steps)
    path <- matrix(
      rep(start_vec, each = n_steps + 1L),
      nrow = n_steps + 1L,
      ncol = n_factors,
      byrow = FALSE
    ) + step_size * tcrossprod(step_seq, direction)
    colnames(path) <- factor_names

    path_list[[d]] <- as.data.frame(path)
  }

  # Warn if many gradients are skipped
  n_valid <- sum(rowSums(is.na(direction_matrix)) == 0)
  if (n_valid < n_draws) {
    warning(
      n_draws - n_valid,
      " draws had near-zero gradients and were skipped."
    )
  }

  # Posterior mean direction
  mean_direction <- stats::setNames(
    colMeans(direction_matrix, na.rm = TRUE),
    factor_names
  )

  # Compute mean path if requested
  mean_path <- NULL
  if (return_mean_path) {
    valid_paths <- path_list[!sapply(path_list, is.null)]
    if (length(valid_paths) > 0) {
      # Build numeric array: step x factor x draw
      path_arrays <- lapply(valid_paths, function(p) {
        as.matrix(p[, factor_names, drop = FALSE])
      })
      arr <- simplify2array(path_arrays)

      # If only one valid path, promote 2D to 3D
      if (length(dim(arr)) == 2) {
        arr <- array(arr, dim = c(nrow(arr), ncol(arr), 1))
      }

      mean_path <- apply(arr, c(1, 2), mean, na.rm = TRUE)
      mean_path <- as.data.frame(mean_path)
      colnames(mean_path) <- factor_names
      mean_path$step <- 0:n_steps
      mean_path <- mean_path[, c("step", factor_names), drop = FALSE]
      rownames(mean_path) <- NULL
    }
  }

  directions_df <- data.frame(
    draw = seq_len(n_draws),
    as.data.frame(direction_matrix),
    row.names = NULL
  )

  mean_direction_df <- data.frame(
    factor = factor_names,
    direction = as.numeric(mean_direction),
    row.names = NULL
  )

  path_indices <- if (drop_null_paths) {
    which(!sapply(path_list, is.null))
  } else {
    seq_along(path_list)
  }

  out_paths <- do.call(
    rbind,
    lapply(path_indices, function(idx) {
      path_df <- path_list[[idx]]
      if (is.null(path_df)) {
        return(NULL)
      }
      path_df$step <- seq_len(nrow(path_df)) - 1L
      path_df$draw <- idx
      path_df[, c("draw", "step", factor_names), drop = FALSE]
    })
  )

  if (is.null(out_paths)) {
    out_paths <- data.frame(
      draw = integer(0),
      step = integer(0)
    )
    for (factor_name in factor_names) {
      out_paths[[factor_name]] <- numeric(0)
    }
  }
  rownames(out_paths) <- NULL

  list(
    directions = directions_df,
    mean_direction = mean_direction_df,
    paths = out_paths,
    mean_path = mean_path
  )
}