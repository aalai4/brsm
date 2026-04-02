#' Predict Response Surface from Posterior Draws or brsm_fit Objects
#'
#' Supports S3 dispatch for `brsm_fit`, `brmsfit`, and draw data frames.
#' For `brsm_fit` objects, `factor_names` are inferred from metadata when
#' omitted.
#'
#' @param draws A `brsm_fit` object, `brmsfit` object, or data frame of draws.
#' @param factor_names Character vector of factor names.
#' @param newdata Data frame of predictor values.
#' @param summary Logical; if `TRUE`, return posterior summaries per row.
#' @param probs Probabilities used for quantile summaries when `summary=TRUE`.
#' @param draw_subset Optional logical or numeric draw subset.
#' @param max_draws Optional cap on number of draws.
#' @param return_matrix Logical; if `TRUE` and `summary=FALSE`, return
#'   draw-by-point matrix instead of long-format data frame.
#' @param output_chunk_size Optional chunk size for long-format assembly.
#'
#' @return A prediction object as matrix, long-format data frame, or summary
#'   data frame depending on options.
#' @export
predict_surface <- function(draws,
                            factor_names,
                            newdata,
                            summary = FALSE,
                            probs = c(0.025, 0.5, 0.975),
                            draw_subset = NULL,
                            max_draws = NULL,
                            return_matrix = FALSE,
                            output_chunk_size = NULL) {
  UseMethod("predict_surface")
}

#' @rdname predict_surface
#' @export
predict_surface.brsm_fit <- function(draws,
                                     factor_names = NULL,
                                     newdata,
                                     summary = FALSE,
                                     probs = c(0.025, 0.5, 0.975),
                                     draw_subset = NULL,
                                     max_draws = NULL,
                                     return_matrix = FALSE,
                                     output_chunk_size = NULL) {
  if (is.null(factor_names)) {
    factor_names <- draws$factor_names
  }
  predict_surface.default(
    draws = as_brsm_draws(draws, factor_names = factor_names),
    factor_names = factor_names,
    newdata = newdata,
    summary = summary,
    probs = probs,
    draw_subset = draw_subset,
    max_draws = max_draws,
    return_matrix = return_matrix,
    output_chunk_size = output_chunk_size
  )
}

#' @rdname predict_surface
#' @export
predict_surface.brmsfit <- function(draws,
                                    factor_names,
                                    newdata,
                                    summary = FALSE,
                                    probs = c(0.025, 0.5, 0.975),
                                    draw_subset = NULL,
                                    max_draws = NULL,
                                    return_matrix = FALSE,
                                    output_chunk_size = NULL) {
  predict_surface.default(
    draws = as_brsm_draws(draws, factor_names = factor_names),
    factor_names = factor_names,
    newdata = newdata,
    summary = summary,
    probs = probs,
    draw_subset = draw_subset,
    max_draws = max_draws,
    return_matrix = return_matrix,
    output_chunk_size = output_chunk_size
  )
}

#' @rdname predict_surface
#' @export
predict_surface.default <- function(draws,
                                    factor_names,
                                    newdata,
                                    summary = FALSE,
                                    probs = c(0.025, 0.5, 0.975),
                                    draw_subset = NULL,
                                    max_draws = NULL,
                                    return_matrix = FALSE,
                                    output_chunk_size = NULL) {
  draws <- .brsm_validate_draws(draws)
  factor_names <- .brsm_validate_factor_names(factor_names)

  if (!is.data.frame(newdata)) {
    stop("newdata must be a data.frame.")
  }

  .brsm_check_columns(
    factor_names, newdata,
    "newdata must contain columns for all factor_names."
  )

  if (nrow(newdata) == 0) {
    stop("newdata must contain at least one row.")
  }

  probs <- .brsm_validate_probs(probs)

  draw_ids <- seq_len(nrow(draws))

  if (!is.null(draw_subset)) {
    if (is.logical(draw_subset)) {
      if (length(draw_subset) != nrow(draws)) {
        stop("Logical draw_subset must have length equal to nrow(draws).")
      }
      draw_ids <- draw_ids[draw_subset]
      draws <- draws[draw_subset, , drop = FALSE]
    } else if (is.numeric(draw_subset)) {
      if (length(draw_subset) == 0 || any(!is.finite(draw_subset))) {
        stop("Numeric draw_subset must contain one or more finite indices.")
      }
      idx <- as.integer(draw_subset)
      if (any(idx < 1L | idx > nrow(draws))) {
        stop("draw_subset indices are out of bounds.")
      }
      draw_ids <- draw_ids[idx]
      draws <- draws[idx, , drop = FALSE]
    } else {
      stop("draw_subset must be NULL, a logical vector, or numeric indices.")
    }
  }

  if (!is.null(max_draws)) {
    if (!is.numeric(max_draws) || length(max_draws) != 1 ||
          !is.finite(max_draws) || max_draws < 1) {
      stop("max_draws must be a finite numeric scalar >= 1.")
    }
    max_draws <- as.integer(max_draws)
    if (nrow(draws) > max_draws) {
      draw_ids <- draw_ids[seq_len(max_draws)]
      draws <- draws[seq_len(max_draws), , drop = FALSE]
    }
  }

  if (nrow(draws) == 0) {
    stop("No draws remain after applying draw_subset/max_draws.")
  }

  n_draws <- nrow(draws)
  n_points <- nrow(newdata)

  draw_cols <- colnames(draws)

  # Warn if newdata columns are not numeric before coercion
  if (!all(sapply(newdata[factor_names], is.numeric))) {
    warning("Some newdata columns were coerced to numeric.")
  }

  # Ensure numeric inputs
  X <- as.data.frame(lapply(newdata[factor_names], as.numeric))

  if (!"b_Intercept" %in% draw_cols) {
    stop("draws must contain 'b_Intercept'.")
  }

  combos <- if (length(factor_names) > 1) {
    utils::combn(factor_names, 2, simplify = FALSE)
  } else {
    NULL
  }

  # Build aligned design and coefficient matrices once, then multiply.
  design_terms <- list(rep(1, n_points))
  coef_terms <- list(draws$b_Intercept)

  # Linear terms
  for (f in factor_names) {
    coef_name <- paste0("b_", f)
    if (coef_name %in% draw_cols) {
      design_terms[[length(design_terms) + 1L]] <- X[[f]]
      coef_terms[[length(coef_terms) + 1L]] <- draws[[coef_name]]
    }
  }

  # Quadratic terms
  for (f in factor_names) {
    coef_name <- .brsm_find_quadratic_col(f, draw_cols)
    if (!is.na(coef_name)) {
      design_terms[[length(design_terms) + 1L]] <- X[[f]]^2
      coef_terms[[length(coef_terms) + 1L]] <- draws[[coef_name]]
    }
  }

  # Interaction terms
  if (!is.null(combos)) {
    for (pair in combos) {
      coef_name <- .brsm_find_interaction_col(pair[1], pair[2], draw_cols)
      if (!is.null(coef_name)) {
        design_terms[[length(design_terms) + 1L]] <- X[[pair[1]]] * X[[pair[2]]]
        coef_terms[[length(coef_terms) + 1L]] <- draws[[coef_name]]
      }
    }
  }

  design_mat <- do.call(cbind, design_terms)
  coef_mat <- do.call(cbind, coef_terms)
  pred_matrix <- coef_mat %*% t(design_mat)

  if (summary && isTRUE(return_matrix)) {
    warning("return_matrix is ignored when summary = TRUE.")
  }

  if (!summary && isTRUE(return_matrix)) {
    return(pred_matrix)
  }

  if (!summary) {
    if (is.null(output_chunk_size)) {
      # Target roughly 1e6 output rows per chunk to limit peak memory churn.
      # For Rcpp path: if total output <= 5M rows, use Rcpp directly.
      # Otherwise fall back to chunked assembly.
      output_chunk_size <- max(1L, as.integer(1e6 / n_points))
    }

    if (!is.numeric(output_chunk_size) || length(output_chunk_size) != 1 ||
          !is.finite(output_chunk_size) || output_chunk_size < 1) {
      stop("output_chunk_size must be a finite numeric scalar >= 1.")
    }

    output_chunk_size <- as.integer(output_chunk_size)

    point_ids <- seq_len(n_points)

    # Use Rcpp for matrix-to-long conversion
    # when output is moderate size.
    # Rcpp avoids rbind and make.unique operations.
    total_rows <- n_draws * n_points
    if (total_rows <= 5e6) {
      # Direct Rcpp conversion
      result <- matrix_to_long_format(pred_matrix, draw_ids, point_ids)
      # Bind with newdata columns (repeated for each draw)
      nd_repeated <- newdata[
        rep(seq_len(n_points), times = n_draws), ,
        drop = FALSE
      ]
      result <- cbind(
        result[, c("draw", "point_id"), drop = FALSE],
        nd_repeated,
        estimate = result$estimate
      )
      rownames(result) <- NULL
      return(result)
    }

    # Chunked assembly for very large outputs (> 5M rows)
    row_ids <- seq_len(n_draws)
    n_chunks <- ceiling(n_draws / output_chunk_size)
    chunks <- vector("list", n_chunks)

    for (k in seq_len(n_chunks)) {
      start_idx <- (k - 1L) * output_chunk_size + 1L
      end_idx <- min(k * output_chunk_size, n_draws)
      row_chunk <- row_ids[start_idx:end_idx]
      draw_chunk <- draw_ids[start_idx:end_idx]
      chunk_n <- length(row_chunk)

      nd_block <- newdata[rep(point_ids, times = chunk_n), , drop = FALSE]
      chunks[[k]] <- data.frame(
        draw = rep(draw_chunk, each = n_points),
        point_id = rep(point_ids, times = chunk_n),
        nd_block,
        estimate = as.vector(t(pred_matrix[row_chunk, , drop = FALSE]))
      )
    }

    result <- do.call(rbind, chunks)
    rownames(result) <- NULL
    return(result)
  }

  # Posterior summary
  summary_mat <- apply(pred_matrix, 2, function(x) {
    qs <- stats::quantile(x, probs = probs, names = FALSE, na.rm = TRUE)
    names(qs) <- paste0("q", formatC(probs * 100, format = "f", digits = 1))
    c(mean = mean(x, na.rm = TRUE), qs)
  })

  summary_mat <- t(summary_mat)

  result <- cbind(newdata, as.data.frame(summary_mat))
  rownames(result) <- NULL

  result
}