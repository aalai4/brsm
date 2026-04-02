gradient_quadratic <- function(draws, x, factor_names, normalize = FALSE) {
  draws <- .brsm_validate_draws(draws)
  factor_names <- .brsm_validate_factor_names(factor_names)

  n_factors <- length(factor_names)
  n_draws <- nrow(draws)

  linear_cols <- paste0("b_", factor_names)
  quad_cols <- vapply(
    factor_names,
    function(f) .brsm_find_quadratic_col(f, colnames(draws)),
    character(1)
  )
  missing_quad <- which(is.na(quad_cols) | quad_cols == "")
  if (length(missing_quad) > 0) {
    stop(
      "draws must contain quadratic terms for: ",
      paste(factor_names[missing_quad], collapse = ", ")
    )
  }

  cols_to_check <- draws[, c(linear_cols, quad_cols), drop = FALSE]

  # Input validation
  if (!is.numeric(x)) stop("x must be numeric.")
  if (anyNA(x)) stop("x cannot contain NA values.")
  if (!all(sapply(cols_to_check, is.numeric))) {
    stop("draws must contain numeric columns for linear and quadratic terms.")
  }

  # Ensure x is matrix
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  n_points <- nrow(x)
  if (ncol(x) != n_factors) {
    stop("number of columns of x must match number of factors.")
  }
  colnames(x) <- factor_names

  # Precompute linear and quadratic matrices
  linear <- as.matrix(draws[, linear_cols, drop = FALSE])
  quad <- as.matrix(draws[, quad_cols, drop = FALSE])

  # Repeat x for all draws
  xmat <- matrix(rep(t(x), each = n_draws),
    nrow = n_draws * n_points,
    ncol = n_factors,
    byrow = TRUE
  )

  # Repeat linear and quad for all points
  linear_big <- linear[rep(1:n_draws, times = n_points), , drop = FALSE]
  quad_big <- quad[rep(1:n_draws, times = n_points), , drop = FALSE]

  # Compute linear + quadratic contribution
  grad <- linear_big + 2 * (quad_big * xmat)

  # Interaction terms
  inter_mat <- matrix(0, nrow = n_draws * n_points, ncol = n_factors)
  if (n_factors > 1) {
    for (i in seq_len(n_factors - 1)) {
      for (j in (i + 1):n_factors) {
        inter_col <- .brsm_find_interaction_col(
          factor_names[i], factor_names[j], colnames(draws)
        )
        coef <- if (!is.null(inter_col)) draws[[inter_col]] else rep(0, n_draws)
        # Repeat coefficients for all points
        coef_big <- rep(coef, times = n_points)
        inter_mat[, i] <- inter_mat[, i] + coef_big * xmat[, j]
        inter_mat[, j] <- inter_mat[, j] + coef_big * xmat[, i]
      }
    }
  }

  grad <- grad + inter_mat
  colnames(grad) <- paste0("d_d", factor_names)

  if (normalize) {
    norms <- sqrt(rowSums(grad^2))
    norms[norms == 0] <- 1
    grad <- grad / norms
  }

  grad_df <- data.frame(
    draw = rep(seq_len(n_draws), each = n_points),
    point_id = rep(seq_len(n_points), times = n_draws),
    x[rep(seq_len(n_points), times = n_draws), , drop = FALSE],
    as.data.frame(grad)
  )
  rownames(grad_df) <- NULL

  grad_df
}