#' Analyze Response Surface Along Ridges
#'
#' Performs ridge analysis by locating maximum and minimum points along
#' specified radii from the center of the experimental region. Computes
#' posterior distributions of ridge values.
#'
#' @param object A `brsm_fit` object, `brmsfit` object, or
#'   data frame of posterior draws with Bayesian coefficient columns.
#' @param factor_names Character vector of factor names
#'   (if object is data frame).
#' @param radii Sequence of radii at which to compute ridge values.
#' @param tol Numerical tolerance for ridge computations.
#' @param summary If TRUE, return summarized ridge values;
#'   if FALSE, return draws.
#'
#' @return Data frame with posterior ridge values or summary statistics.
#'
#' @export
posterior_ridge_analysis <- function(
    object,
    factor_names = NULL,
    radii = seq(0, 3, length.out = 10),
    tol = 1e-6,
    summary = TRUE) {
  UseMethod("posterior_ridge_analysis")
}

#' @rdname posterior_ridge_analysis
#' @export
posterior_ridge_analysis.brsm_fit <- function(
    object,
    factor_names = NULL,
    radii = seq(0, 3, length.out = 10),
    tol = 1e-6,
    summary = TRUE) {
  if (is.null(factor_names)) {
    factor_names <- object$factor_names
  }
  posterior_ridge_analysis.default(
    as_brsm_draws(object),
    factor_names = factor_names,
    radii = radii,
    tol = tol,
    summary = summary
  )
}

#' @rdname posterior_ridge_analysis
#' @export
posterior_ridge_analysis.default <- function(
    object,
    factor_names = NULL,
    radii = seq(0, 3, length.out = 10),
    tol = 1e-6,
    summary = TRUE) {
  draws <- .brsm_validate_draws(
    object,
    empty_message = "draws must contain posterior draws."
  )
  if (is.null(factor_names)) {
    stop("factor_names must be supplied if object does not contain metadata.")
  }
  factor_names <- .brsm_validate_factor_names(factor_names)
  draw_cols <- colnames(draws)
  if (!is.numeric(radii) || length(radii) == 0 ||
        any(!is.finite(radii)) || any(radii < 0)) {
    stop("radii must be a non-empty, finite, non-negative numeric vector.")
  }
  if (!is.finite(tol) || tol <= 0) {
    stop("tol must be a positive, finite number.")
  }

  p <- length(factor_names)
  combos <- if (p > 1) utils::combn(factor_names, 2, simplify = FALSE) else NULL
  factor_index <- stats::setNames(seq_len(p), factor_names)

  # Warn about missing coefficients; for interactions check all naming variants
  quad_cols <- vapply(
    factor_names,
    function(f) .brsm_find_quadratic_col(f, draw_cols),
    character(1)
  )
  base_cols <- c(paste0("b_", factor_names), stats::na.omit(quad_cols))
  missing_cols <- setdiff(base_cols, draw_cols)
  if (any(is.na(quad_cols))) {
    missing_cols <- c(
      missing_cols,
      paste0("b_I(", factor_names[is.na(quad_cols)], "^2)")
    )
  }
  if (!is.null(combos)) {
    for (pair in combos) {
      if (is.null(.brsm_find_interaction_col(pair[1], pair[2], draw_cols))) {
        missing_cols <- c(missing_cols, paste0("b_", pair[1], ":", pair[2]))
      }
    }
  }
  if (length(missing_cols) > 0) {
    warning("Missing coefficients: ", paste(missing_cols, collapse = ", "))
  }

  radii <- sort(radii)
  n_draws <- nrow(draws)
  n_r <- length(radii)
  I <- diag(p)

  linear_cols <- paste0("b_", factor_names)
  b_mat <- matrix(0, nrow = n_draws, ncol = p)
  colnames(b_mat) <- factor_names
  present_linear <- linear_cols %in% draw_cols
  if (any(present_linear)) {
    b_mat[, present_linear] <- as.matrix(
      draws[, linear_cols[present_linear], drop = FALSE]
    )
  }
  storage.mode(b_mat) <- "numeric"

  # Precompute radius-0 stationary points in batch with detailed status.
  h_array <- .brsm_hessian_array(draws, factor_names)
  sp_details <- stationary_points_batch_details(
    h_array = h_array,
    b_matrix = b_mat,
    kappa_thresh = 1e12
  )
  x0_mat <- sp_details$x_star
  x0_status_code <- sp_details$status_code
  x0_kappa_proxy <- sp_details$kappa_proxy

  # storage
  ridge_array <- array(
    NA_real_,
    dim = c(n_draws, n_r, p),
    dimnames = list(
      draw = seq_len(n_draws),
      radius = radii,
      factor = factor_names
    )
  )

  for (d in seq_len(n_draws)) {
    draw <- draws[d, ]
    b <- as.numeric(b_mat[d, ])

    B <- matrix(0, p, p)
    colnames(B) <- factor_names
    rownames(B) <- factor_names

    # quadratic terms
    for (i in seq_len(p)) {
      qname <- .brsm_find_quadratic_col(factor_names[i], draw_cols)
      if (!is.na(qname)) {
        B[i, i] <- draw[[qname]]
      }
    }

    # interactions
    if (p > 1 && !is.null(combos)) {
      for (pair in combos) {
        coef_name <- .brsm_find_interaction_col(pair[1], pair[2], draw_cols)
        if (!is.null(coef_name)) {
          val <- draw[[coef_name]] / 2
          i1 <- factor_index[pair[1]]
          i2 <- factor_index[pair[2]]
          B[i1, i2] <- val
          B[i2, i1] <- val
        }
      }
    }

    eig <- tryCatch(
      eigen(B, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) seq(-1, 1, length.out = p)
    )
    lower <- min(eig) - 100
    upper <- max(eig) + 100

    for (k in seq_along(radii)) {
      r <- radii[k]
      if (r == 0) {
        ridge_array[d, k, ] <- x0_mat[d, ]
        next
      }
      f_root <- function(lambda) {
        A <- B - lambda * I
        sol <- tryCatch(
          solve(A, -0.5 * b),
          error = function(e) rep(NA_real_, p)
        )
        if (any(is.na(sol))) {
          return(NA_real_)
        }
        sum(sol^2) - r^2
      }

      # Search for valid root brackets between eigenvalue poles.
      eig_sorted <- sort(unique(eig))
      eps <- max(
        1e-8,
        sqrt(.Machine$double.eps) * (1 + max(abs(eig_sorted)))
      )
      bounds <- c(min(eig_sorted) - 100, eig_sorted, max(eig_sorted) + 100)

      lambda_candidates <- numeric(0)

      for (ii in seq_len(length(bounds) - 1)) {
        a <- bounds[ii] + eps
        bnd <- bounds[ii + 1] - eps
        if (!is.finite(a) || !is.finite(bnd) || a >= bnd) {
          next
        }

        fa <- f_root(a)
        fb <- f_root(bnd)
        if (!is.finite(fa) || !is.finite(fb)) {
          next
        }

        if (fa == 0) {
          lambda_candidates <- c(lambda_candidates, a)
        }
        if (fb == 0) {
          lambda_candidates <- c(lambda_candidates, bnd)
        }
        if (fa * fb < 0) {
          rt <- tryCatch(
            stats::uniroot(f_root, interval = c(a, bnd), tol = tol)$root,
            error = function(e) NA_real_
          )
          if (!is.na(rt)) {
            lambda_candidates <- c(lambda_candidates, rt)
          }
        }
      }

      if (length(lambda_candidates) > 0) {
        lambda_candidates <- unique(round(lambda_candidates, 12))

        x_candidates <- lapply(lambda_candidates, function(lam) {
          A <- B - lam * I
          tryCatch(
            solve(A, -0.5 * b),
            error = function(e) rep(NA_real_, p)
          )
        })

        # Choose the candidate giving the highest fitted response.
        obj_vals <- vapply(x_candidates, function(xc) {
          if (any(is.na(xc))) {
            return(-Inf)
          }
          as.numeric(sum(b * xc) + t(xc) %*% B %*% xc)
        }, numeric(1))

        best <- which.max(obj_vals)
        if (length(best) == 1 && is.finite(obj_vals[best])) {
          ridge_array[d, k, ] <- x_candidates[[best]]
        }
      }
    }
  }

  # Posterior diagnostics
  ridge_valid <- apply(ridge_array, c(1, 2), function(x) !any(is.na(x)))
  n_valid_draws <- sum(ridge_valid)
  total_draws <- n_draws * n_r
  if (n_valid_draws == 0) {
    warning("No ridge solutions could be computed.")
  } else if (n_valid_draws < total_draws) {
    warning(
      total_draws - n_valid_draws,
      " ridge solutions could not be computed (returned NA)."
    )
  }

  # Per-radius diagnostics
  fail_counts <- colSums(!ridge_valid)
  if (any(fail_counts > 0)) {
    msg <- paste0(
      "Ridge solutions failed for ", fail_counts,
      " draws at radii: ",
      formatC(radii, digits = 6, format = "fg")
    )
    warning(paste(msg[fail_counts > 0], collapse = "; "))
  }

  # Attach per-draw solver diagnostics from radius-0 computation.
  x0_diag_info <- list(
    status_code = x0_status_code,
    status_label = .brsm_stationary_status_labels(x0_status_code),
    kappa_proxy = x0_kappa_proxy,
    status_counts = .brsm_stationary_status_counts(x0_status_code),
    n_draws = n_draws,
    n_excluded = sum(x0_status_code != 0L),
    pct_excluded = mean(x0_status_code != 0L),
    function_name = "posterior_ridge_analysis"
  )

  if (!summary) {
    ridge_df <- do.call(
      rbind,
      lapply(seq_len(n_draws), function(d) {
        draw_df <- data.frame(
          draw = d,
          radius = radii,
          matrix(ridge_array[d, , ], nrow = n_r, ncol = p, byrow = FALSE)
        )
        names(draw_df)[3:(2 + p)] <- factor_names
        draw_df
      })
    )
    rownames(ridge_df) <- NULL
    attr(ridge_df, "diagnostics") <- x0_diag_info
    return(ridge_df)
  }

  # summarize posterior ridge path
  results <- vector("list", n_r)
  for (k in seq_along(radii)) {
    mat <- ridge_array[, k, , drop = FALSE]
    # reshape to n_draws x p even when n_draws == 1
    mat <- matrix(mat, nrow = n_draws, ncol = p, byrow = FALSE)
    colnames(mat) <- factor_names

    summary_mat <- apply(mat, 2, function(x) {
      if (all(is.na(x))) {
        c(
          mean = NA_real_,
          sd = NA_real_,
          q2.5 = NA_real_,
          q50 = NA_real_,
          q97.5 = NA_real_
        )
      } else {
        qs <- stats::quantile(
          x, c(0.025, 0.5, 0.975),
          na.rm = TRUE, names = FALSE
        )
        c(
          mean = mean(x, na.rm = TRUE),
          sd = stats::sd(x, na.rm = TRUE),
          q2.5 = qs[1],
          q50 = qs[2],
          q97.5 = qs[3]
        )
      }
    })

    summary_df <- as.data.frame(t(summary_mat))
    summary_df$factor <- rownames(summary_df)
    summary_df$radius <- radii[k]
    rownames(summary_df) <- NULL
    results[[k]] <- summary_df[
      , c("radius", "factor", "mean", "sd", "q2.5", "q50", "q97.5")
    ]
  }

  result <- do.call(rbind, results)
  rownames(result) <- NULL
  attr(result, "diagnostics") <- x0_diag_info
  result
}