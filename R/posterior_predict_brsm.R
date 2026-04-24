#' Posterior Predictive Draws at New Points
#'
#' Generates posterior predictive values at arbitrary new predictor points.
#' Unlike [predict_surface()], this function can include residual uncertainty
#' (sigma) so outputs represent draws from the posterior predictive
#' distribution, not just posterior expectations.
#'
#' Supports S3 dispatch for \code{brsm_fit}, \code{brmsfit}, and data frames
#' of posterior draws.
#'
#' @param object A \code{brsm_fit} object, \code{brmsfit} object, or data
#'   frame of posterior draws.
#' @param factor_names Character vector of factor names. For
#'   \code{brsm_fit}, defaults to metadata when omitted.
#' @param newdata Data frame of predictor values at which to predict.
#' @param include_residual Logical; if \code{TRUE} (default), include
#'   residual noise using posterior \code{sigma} draws.
#' @param sigma Optional residual SD specification. If \code{NULL}, the
#'   function attempts to infer \code{sigma} from posterior draw columns.
#'   Can be a scalar or one value per posterior draw.
#' @param summary Logical; if \code{TRUE}, return posterior summaries per row.
#' @param probs Probabilities used for quantile summaries when
#'   \code{summary = TRUE}.
#' @param draw_subset Optional logical or numeric draw subset.
#' @param max_draws Optional cap on number of draws.
#' @param return_matrix Logical; if \code{TRUE} and \code{summary = FALSE},
#'   return a draw-by-point matrix instead of long-format data.
#' @param output_chunk_size Optional chunk size for long-format assembly when
#'   \code{summary = FALSE} and \code{return_matrix = FALSE}.
#' @param seed Optional random seed for predictive noise draws.
#'
#' @return A prediction object as matrix, long-format data frame, or summary
#'   data frame depending on options.
#' @export
posterior_predict_brsm <- function(object,
                                   factor_names,
                                   newdata,
                                   include_residual = TRUE,
                                   sigma = NULL,
                                   summary = FALSE,
                                   probs = c(0.025, 0.5, 0.975),
                                   draw_subset = NULL,
                                   max_draws = NULL,
                                   return_matrix = FALSE,
                                   output_chunk_size = NULL,
                                   seed = NULL) {
  UseMethod("posterior_predict_brsm")
}

#' @rdname posterior_predict_brsm
#' @export
posterior_predict_brsm.brsm_fit <- function(object,
                                            factor_names = NULL,
                                            newdata,
                                            include_residual = TRUE,
                                            sigma = NULL,
                                            summary = FALSE,
                                            probs = c(0.025, 0.5, 0.975),
                                            draw_subset = NULL,
                                            max_draws = NULL,
                                            return_matrix = FALSE,
                                            output_chunk_size = NULL,
                                            seed = NULL) {
  if (is.null(factor_names)) {
    factor_names <- object$factor_names
  }

  model_terms <- object$model_terms
  if (is.null(model_terms)) {
    model_terms <- "second_order"
  }

  require_quadratic <- model_terms %in% c("second_order", "pure_quadratic")
  require_interactions <- model_terms %in% c("second_order", "first_order_twi")

  draws_raw <- as.data.frame(object$fit)
  coef_draws <- as_brsm_draws.data.frame(
    draws_raw,
    factor_names = factor_names,
    require_quadratic = require_quadratic,
    require_interactions = require_interactions
  )

  posterior_predict_brsm.default(
    object = coef_draws,
    factor_names = factor_names,
    newdata = newdata,
    include_residual = include_residual,
    sigma = sigma,
    summary = summary,
    probs = probs,
    draw_subset = draw_subset,
    max_draws = max_draws,
    return_matrix = return_matrix,
    output_chunk_size = output_chunk_size,
    seed = seed,
    .sigma_source_draws = draws_raw
  )
}

#' @rdname posterior_predict_brsm
#' @export
posterior_predict_brsm.brmsfit <- function(object,
                                           factor_names,
                                           newdata,
                                           include_residual = TRUE,
                                           sigma = NULL,
                                           summary = FALSE,
                                           probs = c(0.025, 0.5, 0.975),
                                           draw_subset = NULL,
                                           max_draws = NULL,
                                           return_matrix = FALSE,
                                           output_chunk_size = NULL,
                                           seed = NULL) {
  if (is.null(factor_names)) {
    stop("factor_names must be supplied for brmsfit objects.")
  }

  draws_raw <- as.data.frame(object)
  coef_draws <- as_brsm_draws.data.frame(
    draws_raw,
    factor_names = factor_names,
    require_quadratic = FALSE,
    require_interactions = FALSE
  )

  posterior_predict_brsm.default(
    object = coef_draws,
    factor_names = factor_names,
    newdata = newdata,
    include_residual = include_residual,
    sigma = sigma,
    summary = summary,
    probs = probs,
    draw_subset = draw_subset,
    max_draws = max_draws,
    return_matrix = return_matrix,
    output_chunk_size = output_chunk_size,
    seed = seed,
    .sigma_source_draws = draws_raw
  )
}

#' @rdname posterior_predict_brsm
#' @export
posterior_predict_brsm.default <- function(object,
                                           factor_names,
                                           newdata,
                                           include_residual = TRUE,
                                           sigma = NULL,
                                           summary = FALSE,
                                           probs = c(0.025, 0.5, 0.975),
                                           draw_subset = NULL,
                                           max_draws = NULL,
                                           return_matrix = FALSE,
                                           output_chunk_size = NULL,
                                           seed = NULL,
                                           .sigma_source_draws = NULL) {
  coef_draws <- .brsm_validate_draws(object)
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

  if (!is.logical(include_residual) || length(include_residual) != 1L ||
      is.na(include_residual)) {
    stop("include_residual must be TRUE or FALSE.")
  }

  probs <- .brsm_validate_probs(probs)

  selected <- .brsm_select_prediction_draws(
    draws = coef_draws,
    draw_subset = draw_subset,
    max_draws = max_draws
  )

  coef_sel <- selected$draws
  draw_ids <- selected$draw_ids

  if (nrow(coef_sel) == 0) {
    stop("No draws remain after applying draw_subset/max_draws.")
  }

  sigma_source <- if (is.null(.sigma_source_draws)) object else .sigma_source_draws
  sigma_draws <- .brsm_resolve_sigma_draws(
    sigma_source = sigma_source,
    sigma = sigma,
    include_residual = include_residual,
    n_target = nrow(coef_sel),
    draw_ids = draw_ids
  )

  mean_matrix <- predict_surface.default(
    draws = coef_sel,
    factor_names = factor_names,
    newdata = newdata,
    summary = FALSE,
    draw_subset = NULL,
    max_draws = NULL,
    return_matrix = TRUE
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_draws <- nrow(mean_matrix)
  n_points <- ncol(mean_matrix)

  if (isTRUE(include_residual)) {
    noise <- matrix(stats::rnorm(n_draws * n_points),
      nrow = n_draws,
      ncol = n_points
    )
    pred_matrix <- mean_matrix + noise * sigma_draws
  } else {
    pred_matrix <- mean_matrix
  }

  if (summary && isTRUE(return_matrix)) {
    warning("return_matrix is ignored when summary = TRUE.")
  }

  if (!summary && isTRUE(return_matrix)) {
    return(pred_matrix)
  }

  if (!summary) {
    if (is.null(output_chunk_size)) {
      output_chunk_size <- max(1L, as.integer(1e6 / n_points))
    }

    if (!is.numeric(output_chunk_size) || length(output_chunk_size) != 1 ||
          !is.finite(output_chunk_size) || output_chunk_size < 1) {
      stop("output_chunk_size must be a finite numeric scalar >= 1.")
    }

    output_chunk_size <- as.integer(output_chunk_size)
    point_ids <- seq_len(n_points)

    total_rows <- n_draws * n_points
    if (total_rows <= 5e6) {
      result <- matrix_to_long_format(pred_matrix, draw_ids, point_ids)
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


.brsm_select_prediction_draws <- function(draws,
                                          draw_subset = NULL,
                                          max_draws = NULL) {
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

  list(draws = draws, draw_ids = draw_ids)
}


.brsm_resolve_sigma_draws <- function(sigma_source,
                                      sigma = NULL,
                                      include_residual = TRUE,
                                      n_target,
                                      draw_ids) {
  if (!isTRUE(include_residual)) {
    return(rep(0, n_target))
  }

  if (!is.null(sigma)) {
    if (!is.numeric(sigma) || any(!is.finite(sigma)) || any(sigma < 0)) {
      stop("sigma must be a finite non-negative numeric scalar or vector.")
    }
    if (length(sigma) == 1L) {
      return(rep(as.numeric(sigma), n_target))
    }
    if (length(sigma) == n_target) {
      return(as.numeric(sigma))
    }
    stop("sigma vector length must be 1 or match selected number of draws.")
  }

  source_df <- as.data.frame(sigma_source)
  sigma_candidates <- c("sigma", "b_sigma", "sigma_Intercept")
  sigma_col <- intersect(sigma_candidates, colnames(source_df))

  if (length(sigma_col) == 0L) {
    stop(
      "No posterior sigma column found. Provide sigma=... or set ",
      "include_residual = FALSE."
    )
  }

  sigma_draws <- as.numeric(source_df[[sigma_col[[1L]]]])
  if (length(sigma_draws) < max(draw_ids)) {
    stop("Sigma draw vector is shorter than selected draw indices.")
  }

  sigma_draws <- sigma_draws[draw_ids]

  if (any(!is.finite(sigma_draws)) || any(sigma_draws < 0)) {
    stop("Resolved sigma draws must be finite and non-negative.")
  }

  sigma_draws
}