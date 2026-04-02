#' Compute Credible Region Around Response Surface Optimum
#'
#' Identifies the credible region around the stationary point of the response
#' surface at specified probability levels. Useful for uncertainty
#' quantification around the estimated optimum.
#'
#' @param object A `brsm_fit` object, `brmsfit` object, or
#'   data frame of posterior draws with Bayesian coefficient columns.
#' @param factor_names Character vector of factor names
#'   (if object is data frame).
#' @param probs Probability levels for credible intervals
#'   (default: 95% interval).
#' @param summary If TRUE, return region boundaries; if FALSE, return all draws.
#'
#' @return Data frame with credible region boundaries or posterior draws.
#'
#' @export
credible_optimum_region <- function(
    object,
    factor_names = NULL,
    probs = c(0.025, 0.975),
    summary = TRUE) {
  UseMethod("credible_optimum_region")
}

#' @rdname credible_optimum_region
#' @export
credible_optimum_region.brsm_fit <- function(
    object,
    factor_names = NULL,
    probs = c(0.025, 0.975),
    summary = TRUE) {
  if (is.null(factor_names)) {
    factor_names <- object$factor_names
  }
  credible_optimum_region.default(
    as_brsm_draws(object),
    factor_names = factor_names,
    probs = probs,
    summary = summary
  )
}

#' @rdname credible_optimum_region
#' @export
credible_optimum_region.default <- function(
    object,
    factor_names = NULL,
    probs = c(0.025, 0.975),
    summary = TRUE) {
  draws <- .brsm_validate_draws(object)
  if (is.null(factor_names)) {
    stop("factor_names must be supplied if object does not contain metadata.")
  }
  factor_names <- .brsm_validate_factor_names(factor_names)

  probs <- .brsm_validate_probs(
    probs,
    require_finite = TRUE,
    message = "probs must contain finite values between 0 and 1."
  )
  sp <- stationary_point(draws, factor_names)

  sp <- as.data.frame(sp)

  if (ncol(sp) != length(factor_names)) {
    stop("stationary_point() returned unexpected number of columns.")
  }

  if (nrow(sp) == 0) {
    stop("stationary_point() returned zero rows.")
  }

  colnames(sp) <- factor_names

  n_valid <- sum(stats::complete.cases(sp))

  if (n_valid == 0) {
    warning("All stationary points could not be computed.")
  } else if (n_valid < nrow(sp)) {
    warning(nrow(sp) - n_valid, " stationary points could not be computed.")
  }

  if (!summary) {
    return(sp)
  }

  summary_list <- lapply(factor_names, function(f) {
    x <- sp[[f]]

    qs <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
    names(qs) <- paste0("q", formatC(probs * 100, format = "f", digits = 1))

    if (all(is.na(x))) {
      stats <- c(mean = NA_real_, sd = NA_real_)
    } else {
      stats <- c(
        mean = mean(x, na.rm = TRUE),
        sd = stats::sd(x, na.rm = TRUE)
      )
    }

    c(stats, qs)
  })

  summary_df <- as.data.frame(do.call(rbind, summary_list))
  summary_df[] <- lapply(summary_df, as.numeric)

  rownames(summary_df) <- factor_names

  summary_df
}