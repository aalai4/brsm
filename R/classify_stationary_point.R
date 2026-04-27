#' Classify Stationary Points as Maximum, Minimum, or Saddle
#'
#' Classifies stationary points based on the eigenvalues of the Hessian matrix.
#' Uses posterior distribution to compute probabilities of each classification.
#'
#' @param object A `brsm_fit` object, `brmsfit` object, or
#'   data frame of posterior draws with Bayesian coefficient columns.
#' @param factor_names Character vector of factor names
#'   (if object is data frame).
#' @param tol Numerical tolerance for eigenvalue comparisons.
#' @param return_kappa Logical; if \code{TRUE}, attach a \code{kappa} column
#'   to the output containing the LU-diagonal condition-number proxy of the
#'   Hessian for each draw. Useful for identifying near-singular draws without
#'   altering classification behaviour. Default \code{FALSE}.
#'
#' @return A data frame with columns \code{draw} and \code{classification}
#'   (a factor with levels \code{"maximum"}, \code{"minimum"},
#'   \code{"saddle"}, \code{"indeterminate"}). When \code{return_kappa = TRUE}
#'   an additional numeric column \code{kappa} is included.
#'
#' @export
classify_stationary_point <- function(
    object,
    factor_names = NULL,
    tol = 1e-8,
    return_kappa = FALSE) {
  UseMethod("classify_stationary_point")
}

#' @rdname classify_stationary_point
#' @export
classify_stationary_point.brsm_fit <- function(
    object,
    factor_names = NULL,
    tol = 1e-8,
    return_kappa = FALSE) {
  if (is.null(factor_names)) {
    factor_names <- object$factor_names
  }
  classify_stationary_point.default(
    as_brsm_draws(object),
    factor_names = factor_names,
    tol = tol,
    return_kappa = return_kappa
  )
}

#' @rdname classify_stationary_point
#' @export
classify_stationary_point.default <- function(
    object,
    factor_names = NULL,
    tol = 1e-8,
    return_kappa = FALSE) {
  draws <- .brsm_validate_draws(object)
  if (is.null(factor_names)) {
    stop("factor_names must be supplied if object does not contain metadata.")
  }
  factor_names <- .brsm_validate_factor_names(factor_names)

  if (!is.numeric(tol) || length(tol) != 1 || tol < 0) {
    stop("tol must be a non-negative numeric scalar.")
  }

  if (!is.logical(return_kappa) || length(return_kappa) != 1L ||
      is.na(return_kappa)) {
    stop("return_kappa must be TRUE or FALSE.")
  }

  n_draws <- nrow(draws)
  n_factors <- length(factor_names)

  H <- .brsm_hessian_array(draws, factor_names)

  class_vec <- rep(NA_character_, n_draws)
  kappa_vec <- rep(NA_real_, n_draws)

  for (d in seq_len(n_draws)) {
    H_d <- H[d, , ]

    eig <- tryCatch(
      eigen(H_d, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) rep(NA_real_, n_factors)
    )

    if (any(is.na(eig))) {
      next
    }

    # Condition-number proxy: ratio of largest to smallest absolute eigenvalue.
    abs_eig <- abs(eig)
    min_abs <- min(abs_eig)
    max_abs <- max(abs_eig)
    kappa_vec[d] <- if (min_abs > 0) max_abs / min_abs else Inf

    if (any(abs(eig) <= tol)) {
      class_vec[d] <- "indeterminate"
    } else if (all(eig > 0)) {
      class_vec[d] <- "minimum"
    } else if (all(eig < 0)) {
      class_vec[d] <- "maximum"
    } else {
      class_vec[d] <- "saddle"
    }
  }

  class_vec <- factor(
    class_vec,
    levels = c("maximum", "minimum", "saddle", "indeterminate")
  )

  n_na <- sum(is.na(class_vec))

  if (n_na > 0) {
    warning(
      n_na, " draws could not be classified due to ",
      "eigenvalue computation errors."
    )
  }

  out <- data.frame(
    draw = seq_len(n_draws),
    classification = class_vec,
    stringsAsFactors = FALSE
  )

  if (isTRUE(return_kappa)) {
    out$kappa <- kappa_vec
  }

  out
}