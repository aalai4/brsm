#' Canonical Analysis of Bayesian Response Surface
#'
#' Computes the eigendecomposition of the quadratic Hessian (\eqn{\mathbf{B}})
#' across posterior draws. Returns principal curvatures (eigenvalues),
#' canonical axes (eigenvectors), and optional canonical scores at the
#' stationary point.
#'
#' The Hessian \eqn{\mathbf{B}} of the quadratic surface
#' \eqn{y = \beta_0 + \mathbf{b}^\top\mathbf{x} + \mathbf{x}^\top\mathbf{B}\mathbf{x}}
#' has eigendecomposition \eqn{\mathbf{B} = \mathbf{M}\boldsymbol{\Lambda}\mathbf{M}^\top},
#' where the columns of \eqn{\mathbf{M}} are the canonical axes and the
#' diagonal of \eqn{\boldsymbol{\Lambda}} are the principal curvatures.
#'
#' @param object A \code{brsm_fit} object, \code{brmsfit} object, or data frame
#'   of posterior draws with Bayesian coefficient columns.
#' @param factor_names Character vector of factor names (required when
#'   \code{object} is a data frame).
#' @param include_scores Logical; if \code{TRUE} (default), compute canonical
#'   factor scores at the stationary point for each draw.
#' @param kappa_thresh Condition number threshold for Hessian singularity.
#'   Draws whose Hessian condition number exceeds this are excluded from score
#'   computation.
#' @param probs Probability levels for posterior summary intervals.
#' @param summary Logical; if \code{TRUE} (default), return a compact summary
#'   data frame with posterior means and credible intervals per canonical
#'   component. If \code{FALSE}, return per-draw eigenvalues, eigenvectors, and
#'   (optionally) canonical scores as a list.
#'
#' @return When \code{summary = TRUE}: a list with components
#'   \code{eigenvalues} (data frame of posterior summaries for each principal
#'   curvature), \code{eigenvectors} (data frame of posterior summaries for
#'   each axis loading), and optionally \code{scores} (data frame of posterior
#'   summaries for canonical factor scores at \eqn{x^*}).
#'
#'   When \code{summary = FALSE}: a list with components \code{eigenvalues}
#'   (draws × p matrix), \code{eigenvectors} (draws × p × p array), and
#'   optionally \code{scores} (draws × p matrix).
#'
#' @examples
#' \dontrun{
#' fit <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
#' ca <- canonical_analysis(fit)
#' ca$eigenvalues   # posterior summary of principal curvatures
#' ca$eigenvectors  # posterior summary of axis loadings
#' ca$scores        # canonical factor scores at stationary point
#' }
#'
#' @export
canonical_analysis <- function(object,
                               factor_names = NULL,
                               include_scores = TRUE,
                               kappa_thresh = 1e10,
                               probs = c(0.025, 0.5, 0.975),
                               summary = TRUE) {
  UseMethod("canonical_analysis")
}

#' @rdname canonical_analysis
#' @export
canonical_analysis.brsm_fit <- function(object,
                                        factor_names = NULL,
                                        include_scores = TRUE,
                                        kappa_thresh = 1e10,
                                        probs = c(0.025, 0.5, 0.975),
                                        summary = TRUE) {
  if (is.null(factor_names)) {
    factor_names <- object$factor_names
  }
  canonical_analysis.default(
    as_brsm_draws(object),
    factor_names = factor_names,
    include_scores = include_scores,
    kappa_thresh = kappa_thresh,
    probs = probs,
    summary = summary
  )
}

#' @rdname canonical_analysis
#' @export
canonical_analysis.brmsfit <- function(object,
                                       factor_names = NULL,
                                       include_scores = TRUE,
                                       kappa_thresh = 1e10,
                                       probs = c(0.025, 0.5, 0.975),
                                       summary = TRUE) {
  if (is.null(factor_names)) {
    stop("factor_names must be supplied for brmsfit objects.")
  }
  canonical_analysis.default(
    as_brsm_draws(object, factor_names = factor_names),
    factor_names = factor_names,
    include_scores = include_scores,
    kappa_thresh = kappa_thresh,
    probs = probs,
    summary = summary
  )
}

#' @rdname canonical_analysis
#' @export
canonical_analysis.default <- function(object,
                                       factor_names = NULL,
                                       include_scores = TRUE,
                                       kappa_thresh = 1e10,
                                       probs = c(0.025, 0.5, 0.975),
                                       summary = TRUE) {
  draws <- .brsm_validate_draws(object)
  if (is.null(factor_names)) {
    stop("factor_names must be supplied if object does not contain metadata.")
  }
  factor_names <- .brsm_validate_factor_names(factor_names)

  if (!is.logical(include_scores) || length(include_scores) != 1L ||
      is.na(include_scores)) {
    stop("include_scores must be TRUE or FALSE.")
  }
  if (!is.numeric(kappa_thresh) || length(kappa_thresh) != 1L ||
      !is.finite(kappa_thresh) || kappa_thresh <= 0) {
    stop("kappa_thresh must be a finite positive numeric scalar.")
  }
  probs <- .brsm_validate_probs(
    probs,
    require_finite = TRUE,
    message = "probs must be finite values between 0 and 1."
  )
  if (!is.logical(summary) || length(summary) != 1L || is.na(summary)) {
    stop("summary must be TRUE or FALSE.")
  }

  n_draws <- nrow(draws)
  p <- length(factor_names)

  H <- .brsm_hessian_array(draws, factor_names)

  eigenvalues <- matrix(NA_real_, nrow = n_draws, ncol = p)
  colnames(eigenvalues) <- paste0("lambda_", seq_len(p))

  eigenvectors <- array(NA_real_, dim = c(n_draws, p, p))
  dimnames(eigenvectors) <- list(
    draw = NULL,
    factor = factor_names,
    axis = paste0("W_", seq_len(p))
  )

  n_fail <- 0L

  for (d in seq_len(n_draws)) {
    H_d <- H[d, , ]
    eig <- tryCatch(
      eigen(H_d, symmetric = TRUE),
      error = function(e) NULL
    )
    if (is.null(eig) || any(is.na(eig$values)) || any(is.na(eig$vectors))) {
      n_fail <- n_fail + 1L
      next
    }
    eigenvalues[d, ] <- eig$values
    eigenvectors[d, , ] <- eig$vectors
  }

  if (n_fail > 0L) {
    warning(
      n_fail, " draw(s) could not be decomposed and were set to NA."
    )
  }

  # Canonical factor scores at x* per draw
  scores <- NULL
  if (isTRUE(include_scores)) {
    linear_cols <- paste0("b_", factor_names)
    if (!all(linear_cols %in% colnames(draws))) {
      warning(
        "Cannot compute canonical scores: missing linear coefficient ",
        "columns. Set include_scores = FALSE to suppress this warning."
      )
    } else {
      b_mat <- as.matrix(draws[, linear_cols, drop = FALSE])
      storage.mode(b_mat) <- "numeric"

      x_star_mat <- stationary_points_batch(
        h_array = H,
        b_matrix = b_mat,
        kappa_thresh = as.numeric(kappa_thresh)
      )

      scores <- matrix(NA_real_, nrow = n_draws, ncol = p)
      colnames(scores) <- paste0("z_", seq_len(p))

      for (d in seq_len(n_draws)) {
        x_star <- x_star_mat[d, ]
        if (any(!is.finite(x_star))) next
        M <- eigenvectors[d, , ]
        if (any(is.na(M))) next
        scores[d, ] <- as.numeric(t(M) %*% x_star)
      }

      n_score_fail <- sum(rowSums(is.na(scores)) > 0)
      if (n_score_fail > 0L) {
        warning(
          n_score_fail, " draw(s) produced NA canonical scores ",
          "(near-singular Hessian or failed decomposition)."
        )
      }
    }
  }

  if (!isTRUE(summary)) {
    out <- list(
      eigenvalues = eigenvalues,
      eigenvectors = eigenvectors
    )
    if (!is.null(scores)) out$scores <- scores
    return(out)
  }

  # Summarize eigenvalues
  ev_summary <- do.call(rbind, lapply(seq_len(p), function(k) {
    vals <- eigenvalues[, k]
    vals_ok <- vals[is.finite(vals)]
    if (length(vals_ok) == 0L) {
      row <- as.data.frame(matrix(NA_real_, nrow = 1,
        ncol = 2 + length(probs)))
    } else {
      q <- stats::quantile(vals_ok, probs = probs, na.rm = TRUE)
      row <- as.data.frame(t(c(mean = mean(vals_ok), sd = stats::sd(vals_ok),
        q)))
    }
    colnames(row) <- c(
      "mean", "sd",
      paste0("q", formatC(probs * 100, format = "f", digits = 1))
    )
    row$axis <- k
    row$label <- paste0("lambda_", k)
    row
  }))
  ev_summary <- ev_summary[, c("axis", "label",
    setdiff(names(ev_summary), c("axis", "label"))), drop = FALSE]
  rownames(ev_summary) <- NULL

  # Summarize eigenvectors
  evec_rows <- list()
  for (k in seq_len(p)) {
    for (i in seq_len(p)) {
      vals <- eigenvectors[, i, k]
      vals_ok <- vals[is.finite(vals)]
      if (length(vals_ok) == 0L) {
        q_row <- as.data.frame(matrix(NA_real_, nrow = 1,
          ncol = 2 + length(probs)))
        colnames(q_row) <- c(
          "mean", "sd",
          paste0("q", formatC(probs * 100, format = "f", digits = 1))
        )
      } else {
        q <- stats::quantile(vals_ok, probs = probs, na.rm = TRUE)
        q_row <- as.data.frame(t(c(
          mean = mean(vals_ok), sd = stats::sd(vals_ok), q
        )))
        colnames(q_row) <- c(
          "mean", "sd",
          paste0("q", formatC(probs * 100, format = "f", digits = 1))
        )
      }
      q_row$axis <- k
      q_row$factor <- factor_names[i]
      evec_rows[[length(evec_rows) + 1L]] <- q_row
    }
  }
  evec_summary <- do.call(rbind, evec_rows)
  evec_summary <- evec_summary[, c("axis", "factor",
    setdiff(names(evec_summary), c("axis", "factor"))), drop = FALSE]
  rownames(evec_summary) <- NULL

  out <- list(
    eigenvalues = ev_summary,
    eigenvectors = evec_summary
  )

  if (!is.null(scores)) {
    score_rows <- list()
    for (k in seq_len(p)) {
      vals <- scores[, k]
      vals_ok <- vals[is.finite(vals)]
      if (length(vals_ok) == 0L) {
        q_row <- as.data.frame(matrix(NA_real_, nrow = 1,
          ncol = 2 + length(probs)))
        colnames(q_row) <- c(
          "mean", "sd",
          paste0("q", formatC(probs * 100, format = "f", digits = 1))
        )
      } else {
        q <- stats::quantile(vals_ok, probs = probs, na.rm = TRUE)
        q_row <- as.data.frame(t(c(
          mean = mean(vals_ok), sd = stats::sd(vals_ok), q
        )))
        colnames(q_row) <- c(
          "mean", "sd",
          paste0("q", formatC(probs * 100, format = "f", digits = 1))
        )
      }
      q_row$axis <- k
      q_row$label <- paste0("z_", k)
      score_rows[[length(score_rows) + 1L]] <- q_row
    }
    score_summary <- do.call(rbind, score_rows)
    score_summary <- score_summary[, c("axis", "label",
      setdiff(names(score_summary), c("axis", "label"))), drop = FALSE]
    rownames(score_summary) <- NULL
    out$scores <- score_summary
  }

  out
}