#' Identify Stationary Points of Response Surface
#'
#' Computes stationary (critical) points of the response surface from posterior
#' draws, including computation of the eigenvalues of the Hessian matrix for
#' classification.
#'
#' @param object A `brsm_fit` object, `brmsfit` object, or
#'   data frame of posterior draws with Bayesian coefficient columns.
#' @param factor_names Character vector of factor names
#'   (if object is data frame).
#' @param kappa_thresh Threshold for computing condition numbers of Hessian.
#'
#' @return A data frame with stationary points and their properties.
#'
#' @export
stationary_point <- function(object, factor_names = NULL, kappa_thresh = 1e10) {
  UseMethod("stationary_point")
}

#' @rdname stationary_point
#' @export
stationary_point.brsm_fit <- function(
    object,
    factor_names = NULL,
    kappa_thresh = 1e10) {
  if (is.null(factor_names)) {
    factor_names <- object$factor_names
  }
  stationary_point.default(
    as_brsm_draws(object),
    factor_names = factor_names,
    kappa_thresh = kappa_thresh
  )
}

#' @rdname stationary_point
#' @export
stationary_point.default <- function(
    object,
    factor_names = NULL,
    kappa_thresh = 1e10) {
  draws <- .brsm_validate_draws(object)
  if (is.null(factor_names)) {
    stop("factor_names must be supplied if object does not contain metadata.")
  }
  factor_names <- .brsm_validate_factor_names(factor_names)

  n_draws <- nrow(draws)
  n_factors <- length(factor_names)

  # Store column names once
  draw_cols <- colnames(draws)

  # Linear coefficients
  linear_cols <- paste0("b_", factor_names)

  if (!all(linear_cols %in% draw_cols)) {
    stop(
      "missing linear columns: ",
      paste(setdiff(linear_cols, draw_cols), collapse = ", ")
    )
  }

  if (!all(vapply(
    draws[, linear_cols, drop = FALSE],
    is.numeric,
    logical(1)
  ))) {
    stop("draws must contain numeric columns for linear terms.")
  }

  b <- as.matrix(draws[, linear_cols, drop = FALSE])
  storage.mode(b) <- "numeric"

  # Hessian matrices
  H <- .brsm_hessian_array(draws, factor_names)

  # Output matrix
  x_star <- matrix(NA, nrow = n_draws, ncol = n_factors)
  colnames(x_star) <- factor_names

  for (d in seq_len(n_draws)) {
    H_d <- H[d, , ]

    if (kappa(H_d) > kappa_thresh) {
      x_star[d, ] <- NA
      next
    }

    x_star[d, ] <- tryCatch(
      -0.5 * solve(H_d, b[d, ]),
      error = function(e) rep(NA, n_factors)
    )
  }

  n_na <- sum(apply(x_star, 1, function(r) any(is.na(r))))
  if (n_na > 0) {
    warning(n_na, " draws have near-singular Hessians and were set to NA.")
  }

  x_star <- as.data.frame(x_star)
  rownames(x_star) <- NULL

  x_star
}